import pytest
from cobra import Metabolite, Model, Reaction

from ecmapper import ECEvidence, ECMapper, annotate_model, summary, to_dataframe


@pytest.fixture
def model():
    m = Model("test")
    a, b = Metabolite("a_c", compartment="c"), Metabolite("b_c", compartment="c")
    a_e = Metabolite("a_e", compartment="e")
    m.add_metabolites([a, b, a_e])

    def add(rid, annotation=None, notes=None):
        r = Reaction(rid)
        r.add_metabolites({a: -1, b: 1})
        r.annotation = annotation or {}
        r.notes = notes or {}
        m.add_reactions([r])

    add("ANNOTATED", {"ec-code": "2.7.1.11", "kegg.reaction": "R00200"})
    add("NOTES", notes={"EC Number": "3.6.1.1"})
    add("RESOURCES", {"kegg.reaction": "R00200", "biocyc": "META:PEPSYNTH-RXN"})
    add("BIGG_ONLY", {"bigg.reaction": "PGI"})
    add("RHEA_DIRECTED", {"rhea": "16110"})
    add("NOTHING", {"bigg.reaction": "UNKNOWN"})
    m.add_boundary(a_e, type="exchange")
    return m


@pytest.fixture
def mapper(sources):
    return ECMapper(sources.values())


def test_map_model(model, mapper):
    mappings = mapper.map_model(model)
    assert {rid: m.origin for rid, m in mappings.items()} == {
        "ANNOTATED": "annotation",
        "NOTES": "notes",
        "RESOURCES": "resources",
        "BIGG_ONLY": "resources",
        "RHEA_DIRECTED": "resources",
        "NOTHING": "none",
        "EX_a_e": "boundary",
    }
    # the model's EC number is used, the resources are not queried
    assert mappings["ANNOTATED"].evidence == [ECEvidence("2.7.1.11", "annotation")]
    assert mappings["NOTES"].ec_numbers == ["3.6.1.1"]
    assert mappings["NOTHING"].evidence == []


def test_union_with_provenance(model, mapper):
    m = mapper.map_reaction(model.reactions.RESOURCES)
    assert m.ec_numbers == ["2.7.1.40", "2.7.9.2"]
    assert m.ambiguous
    assert m.sources_of("2.7.1.40") == ["kegg", "bkms"]
    assert m.sources_of("2.7.9.2") == ["bkms"]
    assert ECEvidence("2.7.9.2", "bkms", "PEPSYNTH-RXN") in m.evidence
    assert m.xrefs == {"kegg.reaction": ["R00200"], "metacyc.reaction": ["PEPSYNTH-RXN"]}


def test_follow_xrefs(model, sources):
    m = ECMapper(sources.values()).map_reaction(model.reactions.BIGG_ONLY)
    # PGI has no EC number in BiGGr, but cross-references to KEGG, ModelSEED and MetaCyc
    assert m.ec_numbers == ["5.3.1.9"]
    assert set(m.evidence) == {
        ECEvidence("5.3.1.9", "kegg", "R00771", "bigg.reaction:PGI"),
        ECEvidence("5.3.1.9", "modelseed", "rxn00558", "bigg.reaction:PGI"),
        ECEvidence("5.3.1.9", "bkms", "PGLUCISOM-RXN", "bigg.reaction:PGI"),
    }
    m = ECMapper(sources.values(), follow_xrefs=False).map_reaction(model.reactions.BIGG_ONLY)
    assert m.origin == "none"


def test_xrefs_already_in_model_are_not_duplicated(sources):
    r = Reaction("PFK")
    r.annotation = {"bigg.reaction": "PFK", "kegg.reaction": "R00756"}
    m = ECMapper(sources.values()).map_reaction(r)
    kegg = [e for e in m.evidence if e.source == "kegg"]
    assert kegg == [ECEvidence("2.7.1.11", "kegg", "R00756")]
    assert ECEvidence("2.7.1.11", "biggr", "PFK") in m.evidence
    assert ECEvidence("2.7.1.11", "rhea", "16109", "bigg.reaction:PFK") in m.evidence


def test_options(model, sources):
    mappings = ECMapper(sources.values(), always_query=True, use_notes=False, skip_boundary=False).map_model(model)
    annotated = mappings["ANNOTATED"]
    assert annotated.origin == "annotation"
    assert annotated.ec_numbers == ["2.7.1.11", "2.7.1.40"]
    assert annotated.sources == ["annotation", "kegg"]
    assert mappings["NOTES"].origin == "none"
    assert mappings["EX_a_e"].origin == "none"


def test_subset_of_sources(model, sources):
    mapper = ECMapper([sources["kegg"]])
    m = mapper.map_reaction(model.reactions.RESOURCES)
    assert m.sources == ["kegg"]
    assert mapper.map_reaction(model.reactions.BIGG_ONLY).origin == "none"
    assert not sources["bkms"].loaded


def test_sources_by_name(cache_dir):
    mapper = ECMapper(["rhea", "kegg"], cache_dir=cache_dir)
    assert list(mapper.sources) == ["rhea", "kegg.reaction"]
    with pytest.raises(ValueError, match="more than one source"):
        ECMapper(["kegg", "kegg"], cache_dir=cache_dir)


def test_dataframes(model, mapper):
    mappings = mapper.map_model(model)
    wide = to_dataframe(mappings).set_index("reaction_id")
    assert wide.loc["RESOURCES", "ec_numbers"] == "2.7.1.40;2.7.9.2"
    assert wide.loc["RESOURCES", "selection"] == "ambiguous"
    assert wide.loc["RESOURCES", "selected"] == "2.7.1.40;2.7.9.2"
    assert wide.loc["BIGG_ONLY", "selection"] == "single"
    assert wide.loc["EX_a_e", "selection"] == "boundary"
    assert wide.loc["RESOURCES", "n_ec"] == 2
    assert wide.loc["RESOURCES", "sources"] == "kegg;bkms"
    assert wide.loc["RESOURCES", "xrefs"] == "kegg.reaction:R00200;metacyc.reaction:PEPSYNTH-RXN"
    assert wide.loc["NOTHING", "ec_numbers"] == ""
    long = to_dataframe(mappings, long=True)
    rows = long[long.reaction_id == "BIGG_ONLY"]
    assert set(rows.source) == {"kegg", "modelseed", "bkms"}
    assert set(rows.via) == {"bigg.reaction:PGI"}
    assert rows.complete.all()
    assert rows.selected.all()
    assert len(long[long.reaction_id == "NOTHING"]) == 1


def test_summary(model, mapper):
    stats = summary(mapper.map_model(model))
    assert stats["reactions"] == 7
    assert stats["origin"] == {"annotation": 1, "notes": 1, "resources": 3, "none": 1, "boundary": 1}
    assert stats["with_ec"] == 5
    assert stats["ambiguous"] == 1
    assert stats["selection"] == {
        "annotation": 1, "notes": 1, "single": 2, "bkms": 0, "ambiguous": 1, "none": 1, "boundary": 1
    }
    assert stats["reactions_per_source"]["kegg"] == 2


def test_annotate_model(model, mapper):
    mappings = mapper.map_model(model)
    assert annotate_model(model, mappings) == 3
    assert model.reactions.RESOURCES.annotation["ec-code"] == ["2.7.1.40", "2.7.9.2"]
    assert model.reactions.ANNOTATED.annotation["ec-code"] == "2.7.1.11"
    assert "ec-code" not in model.reactions.NOTES.annotation
