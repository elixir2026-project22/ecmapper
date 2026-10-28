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


def test_metabolites(mapper):
    from test_references import pyk

    r = pyk()
    r.annotation = {"kegg.reaction": "R00200", "biocyc": "META:PEPSYNTH-RXN"}
    m = mapper.map_reaction(r)
    # BKMS-react does not agree (PEPSYNTH-RXN: 2.7.1.40 and 2.7.9.2), the metabolites decide
    assert m.ec_numbers == ["2.7.1.40", "2.7.9.2"]
    assert (m.selection, m.selected) == ("metabolites", ["2.7.1.40"])
    assert m.scores["2.7.1.40"] == 1.0 and m.scores["2.7.9.2"] < 0.6

    m = ECMapper(mapper.sources.values(), compare_metabolites=False).map_reaction(r)
    assert (m.selection, m.selected, m.scores) == ("ambiguous", ["2.7.1.40", "2.7.9.2"], {})
    # a too high minimum score or margin leaves it ambiguous
    m = ECMapper(mapper.sources.values(), min_margin=0.5).map_reaction(r)
    assert (m.selection, m.selected) == ("ambiguous", ["2.7.1.40", "2.7.9.2"])


def test_current_ec(mapper, sources):
    m = mapper.map_reaction(reaction({"ec-code": ["4.2.1.52", "2.7.1.40"]}))
    assert (m.selection, m.selected) == ("annotation", ["2.7.1.40", "4.3.3.7"])
    # only deleted EC numbers: the model's EC numbers are kept
    assert mapper.map_reaction(reaction({"ec-code": "1.1.1.74"})).selected == ["1.1.1.74"]
    m = ECMapper(sources.values(), current_ec=False).map_reaction(reaction({"ec-code": "4.2.1.52"}))
    assert m.selected == ["4.2.1.52"]

    evidence = [ECEvidence("4.2.1.52", "modelseed", "rxn1"), ECEvidence("4.3.3.7", "kegg", "R1")]
    selection = mapper.disambiguator.select(ReactionECMapping("R", "resources", evidence))
    assert (selection.selected, selection.selection) == (["4.3.3.7"], "single")


def test_tie_breakers(sources):
    # acetylornithine deacetylase: KEGG R00669 is linked to 3.5.1.14 and 3.5.1.16
    r = reaction({"kegg.reaction": "R00669"})
    assert ECMapper(sources.values()).map_reaction(r).selection == "ambiguous"
    m = ECMapper(sources.values(), tie_breakers=["organism"], organism="ECOLI").map_reaction(r)
    assert (m.selection, m.selected) == ("organism", ["3.5.1.16"])
    # without Swiss-Prot entries of the organism nothing changes
    m = ECMapper(sources.values(), tie_breakers=["organism"], organism="NOSUCH").map_reaction(r)
    assert (m.selection, m.selected) == ("ambiguous", ["3.5.1.14", "3.5.1.16"])

    # via BiGGr, 3.5.1.16 is supported by five sources, 3.5.1.14 by three
    r = reaction({"bigg.reaction": "ACODA"})
    assert ECMapper(sources.values()).map_reaction(r).selection == "ambiguous"
    m = ECMapper(sources.values(), tie_breakers=["sources"]).map_reaction(r)
    assert (m.selection, m.selected) == ("sources", ["3.5.1.16"])

    with pytest.raises(ValueError, match="unknown tie-breaker"):
        ECMapper(sources.values(), tie_breakers=["genes"])
    with pytest.raises(ValueError, match="requires an organism"):
        ECMapper(sources.values(), tie_breakers=["organism"])


def test_resources(sources):
    mapper = ECMapper(sources.values())
    names = [r.name for r in mapper.resources]
    assert names == ["kegg", "modelseed", "rhea", "bkms", "biggr", "enzyme", "kegg-compounds", "biggr-metabolites"]
    # the KEGG source of the mapper is shared with the comparison
    assert mapper.disambiguator.comparison.kegg is sources["kegg"]
    mapper = ECMapper(sources.values(), current_ec=False, compare_metabolites=False)
    assert [r.name for r in mapper.resources] == ["kegg", "modelseed", "rhea", "bkms", "biggr"]
