"""Base class for resources that map reaction identifiers to EC numbers."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import ClassVar, Mapping

import requests

from ..ec import normalize_ec

LOGGER = logging.getLogger(__name__)

USER_AGENT = "ecmapper/0.1 (+https://github.com/cnapy-org/metabolic-pathway-maps)"


def default_cache_dir() -> Path:
    """Cache directory: ``$ECMAPPER_CACHE`` or ``$XDG_CACHE_HOME/ecmapper`` (``~/.cache/ecmapper``)."""
    if "ECMAPPER_CACHE" in os.environ:
        return Path(os.environ["ECMAPPER_CACHE"])
    xdg = os.environ.get("XDG_CACHE_HOME")
    return (Path(xdg) if xdg else Path.home() / ".cache") / "ecmapper"


#: values of :attr:`FileStatus.status`
FILE_STATUSES = ("current", "outdated", "missing", "updated", "downloaded", "local", "error")


@dataclass
class FileStatus:
    """Result of checking (or updating) one cached bulk file.

    Attributes
    ----------
    source : str
        Name of the source.
    file : str
        Name of the file in the cache.
    status : str
        ``current`` (the cached file is up to date), ``outdated`` (a newer
        version is available), ``missing`` (not downloaded yet), ``updated``
        (the newer version was downloaded), ``downloaded`` (a missing file was
        downloaded), ``local`` (a user supplied file that is not managed) or
        ``error`` (the check failed, see `detail`).
    cached : str
        When the cached file was downloaded (UTC, ISO format), empty if unknown.
    remote_modified : str
        ``Last-Modified`` of the file on the server, if the server sends it.
    detail : str
        Error message.
    """

    source: str
    file: str
    status: str
    cached: str = ""
    remote_modified: str = ""
    detail: str = ""

    @property
    def update_available(self) -> bool:
        return self.status in ("outdated", "missing")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ECSource(ABC):
    """A resource with a bulk table that maps identifiers of one namespace to EC numbers.

    The bulk files are downloaded on first use and cached in `cache_dir`. If a
    download fails, the source logs a warning and behaves as if it was empty.
    Next to every downloaded file, its ETag, Last-Modified date and SHA-256
    checksum are stored (``<file>.meta.json``); :meth:`check_updates` uses them
    to find out whether a newer version is available without downloading it
    again and :meth:`update` downloads only the files that changed.

    Parameters
    ----------
    cache_dir : path, optional
        Directory for the downloaded files, defaults to :func:`default_cache_dir`.
    refresh : bool
        Download the files again on the next load even if they are in the cache.
    local_files : mapping, optional
        Use these local files instead of downloading, keyed like :attr:`files`.
    session : requests.Session, optional
        Session used for the downloads.
    """

    #: short name of the resource, used in the evidence of a mapping
    name: ClassVar[str]
    #: canonical annotation namespace of the identifiers this source resolves
    namespace: ClassVar[str]
    #: bulk files {file name in cache: URL}
    files: ClassVar[dict[str, str]]

    def __init__(
        self,
        cache_dir: str | os.PathLike | None = None,
        refresh: bool = False,
        local_files: Mapping[str, str | os.PathLike] | None = None,
        session: requests.Session | None = None,
    ):
        self.cache_dir = Path(cache_dir) if cache_dir is not None else default_cache_dir()
        self.refresh = refresh
        self.local_files = {k: Path(v) for k, v in (local_files or {}).items()}
        self.session = session
        self.error: Exception | None = None
        self._table: dict[str, set[str]] | None = None
        self._xrefs: dict[str, dict[str, set[str]]] = {}

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.name}: {self.namespace}>"

    @property
    def table(self) -> dict[str, set[str]]:
        """The mapping identifier -> EC numbers, loaded on first access."""
        if self._table is None:
            self.load()
        return self._table

    @property
    def loaded(self) -> bool:
        """True once the bulk files were loaded (or loading failed)."""
        return self._table is not None

    @property
    def available(self) -> bool:
        """False if the bulk files could not be loaded."""
        if self._table is None:
            self.load()
        return self.error is None

    def load(self) -> None:
        """Download (if needed) and parse the bulk files."""
        try:
            paths = {name: self._fetch(name, url) for name, url in self.files.items()}
            table = self._parse(paths)
        except (OSError, ValueError, requests.RequestException) as e:
            LOGGER.warning("%s is not available and will be skipped: %s", self.name, e)
            self.error = e
            self._table, self._xrefs = {}, {}
            return
        self.refresh = False
        self._table = {k: v for k, v in table.items() if v}
        LOGGER.info("%s: loaded EC numbers for %d identifiers", self.name, len(self._table))

    def lookup(self, identifier: str) -> set[str]:
        """EC numbers for `identifier` (empty if unknown)."""
        return set(self.table.get(identifier, ()))

    def xrefs(self, identifier: str) -> dict[str, set[str]]:
        """Cross-references of `identifier` to other namespaces, if the resource provides them."""
        if self._table is None:
            self.load()
        return {ns: set(ids) for ns, ids in self._xrefs.get(identifier, {}).items()}

    @abstractmethod
    def _parse(self, paths: dict[str, Path]) -> dict[str, set[str]]:
        """Parse the bulk files into {identifier: {EC numbers}}."""

    @staticmethod
    def _add(table: dict[str, set[str]], identifier: str, raw_ec: str) -> None:
        """Add `raw_ec` to `table` if it is a valid EC number."""
        ec = normalize_ec(raw_ec)
        if ec is not None:
            table.setdefault(identifier, set()).add(ec)

    def check_updates(self) -> list[FileStatus]:
        """Check for every bulk file whether a newer version is available, without updating the cache.

        Uses conditional requests (``If-None-Match``/``If-Modified-Since``), so
        that nothing is downloaded if the file did not change. If the server
        sends neither ETag nor Last-Modified (KEGG), the file is downloaded and
        compared with the cached one by checksum.
        """
        return [self._sync(name, url, download=False) for name, url in self.files.items()]

    def update(self, force: bool = False) -> list[FileStatus]:
        """Download the bulk files for which a newer version is available (or that are missing).

        With `force` all files are downloaded again. If a file changed, the
        source is reloaded on its next use.
        """
        statuses = [self._sync(name, url, download=True, force=force) for name, url in self.files.items()]
        # the files were just synchronised, a pending refresh would download them again
        self.refresh = False
        if any(s.status in ("updated", "downloaded") for s in statuses):
            self.reset()
        return statuses

    def reset(self) -> None:
        """Forget the loaded table, it is loaded again on the next use."""
        self._table, self._xrefs, self.error = None, {}, None

    def _meta_path(self, name: str) -> Path:
        return self.cache_dir / f"{name}.meta.json"

    def _read_meta(self, name: str) -> dict:
        try:
            with open(self._meta_path(name), encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _write_meta(self, name: str, meta: dict) -> None:
        with open(self._meta_path(name), "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=1)

    def _request(self, url: str, headers: dict | None = None) -> requests.Response:
        session = self.session or requests.Session()
        return session.get(
            url, stream=True, timeout=(30, 600), headers={"User-Agent": USER_AGENT, **(headers or {})}
        )

    def _save_temp(self, name: str, response: requests.Response) -> tuple[str, str]:
        """Write the body of `response` to a temporary file in the cache, return its path and checksum."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        # write to a temporary file first so that an interrupted download is not cached
        fd, tmp = tempfile.mkstemp(dir=self.cache_dir, prefix=f".{name}.")
        digest = hashlib.sha256()
        try:
            with os.fdopen(fd, "wb") as f:
                for chunk in response.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
                    digest.update(chunk)
        except BaseException:
            os.unlink(tmp)
            raise
        return tmp, digest.hexdigest()

    @staticmethod
    def _new_meta(url: str, response: requests.Response, sha256: str) -> dict:
        return {
            "url": url,
            "etag": response.headers.get("ETag"),
            "last_modified": response.headers.get("Last-Modified"),
            "sha256": sha256,
            "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }

    def _download(self, name: str, url: str) -> Path:
        LOGGER.info("Downloading %s", url)
        path = self.cache_dir / name
        with self._request(url) as r:
            r.raise_for_status()
            tmp, sha256 = self._save_temp(name, r)
            os.replace(tmp, path)
            self._write_meta(name, self._new_meta(url, r, sha256))
        return path

    def _fetch(self, name: str, url: str) -> Path:
        if name in self.local_files:
            return self.local_files[name]
        path = self.cache_dir / name
        if path.exists() and not self.refresh:
            return path
        return self._download(name, url)

    def _sync(self, name: str, url: str, download: bool, force: bool = False) -> FileStatus:
        """Compare the cached file `name` with the server and download it if `download` and it changed."""
        status = FileStatus(self.name, name, "error")
        if name in self.local_files:
            status.status = "local"
            return status
        path = self.cache_dir / name
        try:
            if not path.exists():
                if download:
                    self._download(name, url)
                    status.status = "downloaded"
                    status.cached = self._read_meta(name).get("downloaded", "")
                else:
                    status.status = "missing"
                return status
            meta = self._read_meta(name)
            status.cached = meta.get("downloaded") or datetime.fromtimestamp(
                path.stat().st_mtime, timezone.utc
            ).isoformat(timespec="seconds")
            local_sha256 = _sha256(path)
            # validators describe the downloaded version, they are useless if the file was changed since
            intact = meta.get("sha256") == local_sha256
            headers = {}
            if not force and intact:
                if meta.get("etag"):
                    headers["If-None-Match"] = meta["etag"]
                if meta.get("last_modified"):
                    headers["If-Modified-Since"] = meta["last_modified"]
            with self._request(url, headers) as r:
                if r.status_code == 304:
                    status.status = "current"
                    status.remote_modified = r.headers.get("Last-Modified") or meta.get("last_modified") or ""
                    return status
                r.raise_for_status()
                status.remote_modified = r.headers.get("Last-Modified") or ""
                changed = self._changed(meta, r) if intact and not force else None
                if changed is not None and not download:
                    status.status = "outdated" if changed else "current"
                    return status
                if changed is False:
                    status.status = "current"
                    return status
                # no usable validators (or forced): compare the content
                tmp, sha256 = self._save_temp(name, r)
                new_meta = self._new_meta(url, r, sha256)
            if sha256 == local_sha256:
                os.unlink(tmp)
                # remember the validators so that the next check needs no download
                self._write_meta(name, {**new_meta, "downloaded": status.cached})
                status.status = "current"
            elif download:
                os.replace(tmp, path)
                self._write_meta(name, new_meta)
                status.status, status.cached = "updated", new_meta["downloaded"]
            else:
                os.unlink(tmp)
                status.status = "outdated"
        except (OSError, requests.RequestException) as e:
            LOGGER.warning("%s: checking %s failed: %s", self.name, name, e)
            status.status, status.detail = "error", str(e)
        return status

    @staticmethod
    def _changed(meta: dict, response: requests.Response) -> bool | None:
        """Whether the server's validators differ from the stored ones, None if they can not be compared."""
        etag = response.headers.get("ETag")
        if etag and meta.get("etag"):
            return etag != meta["etag"]
        last_modified = response.headers.get("Last-Modified")
        if last_modified and meta.get("last_modified"):
            return last_modified != meta["last_modified"]
        return None
