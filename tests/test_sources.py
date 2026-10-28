import logging

import pytest
import requests

from ecmapper.sources import (
    BiGGrSource,
    KEGGSource,
    default_sources,
)


def test_kegg(sources):
    kegg = sources["kegg"]
    assert kegg.lookup("R00200") == {"2.7.1.40"}
    assert kegg.lookup("R00669") == {"3.5.1.14", "3.5.1.16"}
    assert kegg.lookup("R99999") == set()


def test_modelseed(sources):
    seed = sources["modelseed"]
    assert seed.lookup("rxn00469") == {"3.5.1.14", "3.5.1.16"}
    # entries with whitespace or EC- prefix are cleaned up
    assert seed.lookup("rxn02157") == {"1.14.11.13"}
    assert seed.lookup("rxn05015") == {"2.3.1.122"}
    # letter suffixes and transporter classification numbers are no EC numbers
    assert seed.lookup("rxn00062") == set()
    assert seed.lookup("rxn05145") == set()


def test_rhea_directions(sources):
    rhea = sources["rhea"]
    for rhea_id in ("16109", "16110", "16111", "16112"):
        assert rhea.lookup(rhea_id) == {"2.7.1.11"}
    # master reaction without EC number
    assert rhea.lookup("40843") == set()
    assert rhea.lookup("40844") == set()


def test_bkms(sources):
    bkms = sources["bkms"]
    assert bkms.lookup("PGLUCISOM-RXN") == {"5.3.1.9"}
    assert bkms.lookup("PEPSYNTH-RXN") == {"2.7.1.40", "2.7.9.2"}
    # rows with several MetaCyc IDs, BRENDA provisional numbers
    assert bkms.lookup("RXN-7738") == {"3.1.1.14"}
    assert bkms.lookup("CHLOROPHYLLASE-RXN") == {"3.1.1.14"}
    assert "3.5.1.B4" in bkms.lookup("ACETYLORNDEACET-RXN")
    # incomplete EC numbers are padded
    assert bkms.lookup("THYMIDYLATE-5-PHOSPHATASE-RXN") == {"3.1.3.-"}
    # SPONTANEOUS and empty EC numbers are skipped
    assert bkms.lookup("RXN-3521") == set()
    assert bkms.lookup("PRENYL-PYROPHOSPHATASE-RXN") == set()


def test_biggr(sources):
    biggr = sources["biggr"]
    assert biggr.lookup("PFK") == {"2.7.1.11"}
    assert biggr.lookup("PGI") == set()
    # legacy IDs are resolved
    assert biggr.lookup("ACODA_1") == {"3.5.1.16"}
    xrefs = biggr.xrefs("PGI")
    assert xrefs["kegg.reaction"] == {"R00771"}
    assert xrefs["metacyc.reaction"] == {"PGLUCISOM-RXN"}
    # cross-references from obsolete ModelSEED reactions are dropped by default
    assert xrefs["seed.reaction"] == {"rxn00558", "rxn33838"}
    assert biggr.xrefs("UNKNOWN") == {}


def test_biggr_include_obsolete(cache_dir):
    biggr = BiGGrSource(cache_dir=cache_dir, include_obsolete=True)
    assert biggr.xrefs("PGI")["seed.reaction"] == {"rxn00558", "rxn33838", "rxn34259", "rxn37760"}


def test_local_files(tmp_path, cache_dir):
    kegg = KEGGSource(cache_dir=tmp_path / "empty", local_files={"kegg_link_ec_rn.tsv": cache_dir / "kegg_link_ec_rn.tsv"})
    assert kegg.lookup("R00200") == {"2.7.1.40"}


class FakeResponse:
    def __init__(self, content=b"", error=None):
        self.content, self.error = content, error

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def raise_for_status(self):
        if self.error:
            raise self.error

    def iter_content(self, chunk_size):
        yield self.content


def test_download_and_cache(tmp_path, monkeypatch):
    calls = []

    def get(session, url, **kwargs):
        calls.append(url)
        return FakeResponse(b"rn:R00200\tec:2.7.1.40\n")

    monkeypatch.setattr(requests.Session, "get", get)
    assert KEGGSource(cache_dir=tmp_path).lookup("R00200") == {"2.7.1.40"}
    assert calls == [KEGGSource.files["kegg_link_ec_rn.tsv"]]
    assert (tmp_path / "kegg_link_ec_rn.tsv").exists()
    # the second instance uses the cache
    assert KEGGSource(cache_dir=tmp_path).lookup("R00200") == {"2.7.1.40"}
    assert len(calls) == 1
    # unless refresh is set
    KEGGSource(cache_dir=tmp_path, refresh=True).lookup("R00200")
    assert len(calls) == 2


def test_failed_download(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(
        requests.Session, "get", lambda *a, **k: FakeResponse(error=requests.HTTPError("403 Forbidden"))
    )
    kegg = KEGGSource(cache_dir=tmp_path)
    with caplog.at_level(logging.WARNING):
        assert kegg.lookup("R00200") == set()
    assert not kegg.available
    assert "kegg is not available" in caplog.text
    # nothing (and no partial file) is cached
    assert list(tmp_path.iterdir()) == []


def test_default_sources(tmp_path):
    assert [s.name for s in default_sources(cache_dir=tmp_path)] == ["kegg", "modelseed", "rhea", "bkms", "biggr"]
    assert [s.name for s in default_sources(["rhea"], cache_dir=tmp_path)] == ["rhea"]
    with pytest.raises(ValueError, match="unknown source"):
        default_sources(["chebi"])
