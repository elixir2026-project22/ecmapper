"""Tests against the real resources, run with ``pytest -m network``."""

from pathlib import Path

import cobra
import pytest

from ecmapper import ECMapper, check_updates, summary, update_sources
from ecmapper.sources import SOURCES

pytestmark = pytest.mark.network

E_COLI_CORE = Path(__file__).parents[2] / "models" / "bigg2-ecoli" / "e_coli_core.xml"


@pytest.fixture(scope="module")
def cache_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("cache")


@pytest.mark.parametrize(
    "name, identifier, ec",
    [
        ("kegg", "R00200", "2.7.1.40"),
        ("modelseed", "rxn00001", "3.6.1.1"),
        ("rhea", "16110", "2.7.1.11"),
        ("bkms", "PGLUCISOM-RXN", "5.3.1.9"),
        ("biggr", "PFK", "2.7.1.11"),
    ],
)
def test_source(cache_dir, name, identifier, ec):
    source = SOURCES[name](cache_dir=cache_dir)
    assert ec in source.lookup(identifier)
    assert source.available
    assert len(source.table) > 1000


@pytest.mark.skipif(not E_COLI_CORE.exists(), reason="e_coli_core.xml not found")
def test_e_coli_core(cache_dir):
    model = cobra.io.read_sbml_model(str(E_COLI_CORE))
    mappings = ECMapper(cache_dir=cache_dir).map_model(model)
    stats = summary(mappings)
    assert stats["origin"]["boundary"] == 20
    # the model is only annotated with BiGG IDs; all enzymatic reactions except CYTBD and
    # the PTS reactions are found (53 in 10/2026), unmapped are mostly transport reactions
    assert stats["origin"]["resources"] >= 50
    assert "5.3.1.9" in mappings["PGI"].ec_numbers
    assert "2.7.1.11" in mappings["PFK"].ec_numbers


def test_check_updates(cache_dir):
    assert all(s.status in ("current", "downloaded", "updated") for s in update_sources(cache_dir=cache_dir))
    # all servers answer either with 304 Not Modified or (KEGG) with identical content
    assert [s.status for s in check_updates(cache_dir=cache_dir)] == ["current"] * 6
