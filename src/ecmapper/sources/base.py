"""Base class for resources that map reaction identifiers to EC numbers."""

from __future__ import annotations

import logging
import os
import tempfile
from abc import ABC, abstractmethod
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


class ECSource(ABC):
    """A resource with a bulk table that maps identifiers of one namespace to EC numbers.

    The bulk files are downloaded on first use and cached in `cache_dir`. If a
    download fails, the source logs a warning and behaves as if it was empty.

    Parameters
    ----------
    cache_dir : path, optional
        Directory for the downloaded files, defaults to :func:`default_cache_dir`.
    refresh : bool
        Download the files again even if they are in the cache.
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

    def _fetch(self, name: str, url: str) -> Path:
        if name in self.local_files:
            return self.local_files[name]
        path = self.cache_dir / name
        if path.exists() and not self.refresh:
            return path
        LOGGER.info("Downloading %s", url)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        session = self.session or requests.Session()
        with session.get(url, stream=True, timeout=(30, 600), headers={"User-Agent": USER_AGENT}) as r:
            r.raise_for_status()
            # write to a temporary file first so that an interrupted download is not cached
            fd, tmp = tempfile.mkstemp(dir=self.cache_dir, prefix=f".{name}.")
            try:
                with os.fdopen(fd, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                os.replace(tmp, path)
            except BaseException:
                os.unlink(tmp)
                raise
        return path
