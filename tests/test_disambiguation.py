import pytest
from cobra import Reaction

from ecmapper import ECEvidence, ECMapper, ReactionECMapping
from ecmapper.disambiguation import covers, drop_covered, select


def reaction(annotation=None, notes=None, rid="R1"):
    r = Reaction(rid)
    r.annotation = annotation or {}
    r.notes = notes or {}
    return r


@pytest.fixture
def mapper(sources):
    return ECMapper(sources.values())


def test_covers():
    assert covers("1.1.1.-", "1.1.1.1")
    assert covers("1.1.-.-", "1.1.1.-")
    assert not covers("1.1.1.-", "1.1.2.1")
    assert not covers("1.1.1.1", "1.1.1.1")
    assert not covers("1.1.1.1", "1.1.1.-")


def test_drop_covered():
    assert drop_covered(["2.7.1.-", "2.7.1.165", "2.7.1.31"]) == ["2.7.1.165", "2.7.1.31"]
    assert drop_covered(["1.1.-.-", "1.1.1.-"]) == ["1.1.1.-"]
    assert drop_covered(["1.1.1.-", "1.2.1.1"]) == ["1.1.1.-", "1.2.1.1"]


def test_bkms_consensus(mapper):
    # glycerate 2-kinase: ModelSEED rxn08647 has 2.7.1.-, 2.7.1.165 and 2.7.1.31,
    # BKMS-react aligns the KEGG reaction R08572 only with 2.7.1.165
    m = mapper.map_reaction(reaction({"seed.reaction": "rxn08647", "kegg.reaction": "R08572"}))
    assert m.ec_numbers == ["2.7.1.-", "2.7.1.165", "2.7.1.31"]
    assert (m.selection, m.selected) == ("bkms", ["2.7.1.165"])
    assert not m.ambiguous
    # the same with the MetaCyc reaction
    m = mapper.map_reaction(reaction({"seed.reaction": "rxn08647", "biocyc": "META:GKI-RXN"}))
    assert (m.selection, m.selected) == ("bkms", ["2.7.1.165"])


def test_bkms_consensus_disabled(sources):
    m = ECMapper(sources.values(), bkms_consensus=False).map_reaction(
        reaction({"seed.reaction": "rxn08647", "kegg.reaction": "R08572"})
    )
    assert (m.selection, m.selected) == ("ambiguous", ["2.7.1.165", "2.7.1.31"])
    assert m.ambiguous


def test_bkms_without_consensus_stays_ambiguous(mapper):
    m = mapper.map_reaction(reaction({"kegg.reaction": "R00200", "biocyc": "META:PEPSYNTH-RXN"}))
    assert (m.selection, m.selected) == ("ambiguous", ["2.7.1.40", "2.7.9.2"])


def test_single(mapper):
    m = mapper.map_reaction(reaction({"kegg.reaction": "R00200"}))
    assert (m.selection, m.selected) == ("single", ["2.7.1.40"])


def test_model_ec_numbers_are_kept(mapper):
    m = mapper.map_reaction(reaction({"ec-code": ["2.7.1.11", "2.7.1.40"]}))
    assert (m.selection, m.selected) == ("annotation", ["2.7.1.11", "2.7.1.40"])
    m = mapper.map_reaction(reaction(notes={"EC Number": "3.6.1.1"}))
    assert (m.selection, m.selected) == ("notes", ["3.6.1.1"])
    # with always_query, the resources do not change the selection
    m = ECMapper(mapper.sources.values(), always_query=True).map_reaction(
        reaction({"ec-code": "2.7.1.11", "kegg.reaction": "R00200"})
    )
    assert m.ec_numbers == ["2.7.1.11", "2.7.1.40"]
    assert (m.selection, m.selected) == ("annotation", ["2.7.1.11"])


def test_none_and_boundary(mapper):
    m = mapper.map_reaction(reaction({"kegg.reaction": "R99999"}))
    assert (m.selection, m.selected) == ("none", [])
    m = mapper.map_reaction(reaction({"kegg.reaction": "R00200"}), boundary=True)
    assert (m.selection, m.selected) == ("boundary", [])


def test_identifiers(mapper):
    m = mapper.map_reaction(reaction({"bigg.reaction": "PGI"}))
    assert m.identifiers["bigg.reaction"] == ["PGI"]
    assert m.identifiers["kegg.reaction"] == ["R00771"]
    assert m.identifiers["metacyc.reaction"] == ["PGLUCISOM-RXN"]
    # namespaces without a source are recorded as well
    assert "metanetx.reaction" in m.identifiers
    assert ReactionECMapping("R", "none").identifiers == {}


def test_select_without_bkms():
    evidence = [ECEvidence("2.7.1.165", "kegg", "R08572"), ECEvidence("2.7.1.31", "modelseed", "rxn08647")]
    m = ReactionECMapping("R", "resources", evidence, identifiers={"kegg.reaction": ["R08572"]})
    assert select(m) == (["2.7.1.165", "2.7.1.31"], "ambiguous")
