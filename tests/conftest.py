import shutil
import tarfile
from pathlib import Path

import pytest

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
