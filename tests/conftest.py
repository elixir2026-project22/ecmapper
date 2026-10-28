import shutil
import tarfile
from pathlib import Path

import pytest
import requests

from ecmapper.sources import SOURCES

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def cache_dir(tmp_path) -> Path:
    """A download cache that holds the fixture excerpts of all bulk files."""
    cache = tmp_path / "cache"
    cache.mkdir()
    for name in ("kegg_link_ec_rn.tsv", "Unique_ModelSEED_Reaction_ECs.txt", "rhea2ec.tsv",
                 "rhea-directions.tsv", "biggr_reactions.json"):
        shutil.copy(FIXTURES / name, cache / name)
    # the BKMS-react download is a tar.gz with a tab separated Reactions_BKMS.csv
    with tarfile.open(cache / "Reactions_BKMS.tar.gz", "w:gz") as tar:
        tar.add(FIXTURES / "Reactions_BKMS.csv", arcname="Reactions_BKMS.csv")
    return cache


@pytest.fixture
def sources(cache_dir):
    """All sources, reading the fixture files from the cache."""
    return {name: cls(cache_dir=cache_dir) for name, cls in SOURCES.items()}


@pytest.fixture(autouse=True)
def no_downloads(request, monkeypatch):
    """Fail if a test that is not marked as network test tries to download."""
    if request.node.get_closest_marker("network"):
        return

    def fail(*args, **kwargs):
        raise AssertionError("unexpected download in offline test")

    monkeypatch.setattr("requests.Session.get", fail)


class FakeResponse:
    """Minimal stand-in for a streamed requests.Response."""

    def __init__(self, content=b"", headers=None, status_code=200, error=None):
        self.content, self.headers, self.status_code, self.error = content, headers or {}, status_code, error

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def raise_for_status(self):
        if self.error:
            raise self.error

    def iter_content(self, chunk_size):
        yield self.content


class FakeServer:
    """Serves files {url: (content, headers)} and answers conditional requests like a web server."""

    def __init__(self, files=None):
        self.files = dict(files or {})
        self.requests = []

    def get(self, url, headers=None, **kwargs):
        headers = headers or {}
        self.requests.append((url, headers))
        if url not in self.files:
            return FakeResponse(status_code=404, error=requests.HTTPError(f"404 Not Found: {url}"))
        content, response_headers = self.files[url]
        etag, last_modified = response_headers.get("ETag"), response_headers.get("Last-Modified")
        if etag and headers.get("If-None-Match") == etag:
            return FakeResponse(status_code=304)
        if not etag and last_modified and headers.get("If-Modified-Since") == last_modified:
            return FakeResponse(status_code=304)
        return FakeResponse(content, response_headers)


@pytest.fixture
def server(monkeypatch):
    """A fake server that replaces requests.Session.get."""
    fake = FakeServer()
    monkeypatch.setattr("requests.Session.get", lambda session, url, **kwargs: fake.get(url, **kwargs))
    return fake
