from cobra import Reaction

from ecmapper.annotation import annotated_ecs, normalize_xref, notes_ecs, xrefs


def reaction(annotation=None, notes=None):
    r = Reaction("R1")
    r.annotation = annotation or {}
    r.notes = notes or {}
    return r


def test_annotated_ecs_str_and_list():
    assert annotated_ecs(reaction({"ec-code": "2.7.1.11"})) == ["2.7.1.11"]
    assert annotated_ecs(reaction({"ec-code": ["2.7.1.11", "2.7.1.11", "x"], "brenda": "2.7.1.40"})) == [
        "2.7.1.11",
        "2.7.1.40",
    ]
    assert annotated_ecs(reaction({"kegg.reaction": "R00200"})) == []


def test_notes_ecs():
    assert notes_ecs(reaction(notes={"EC Number": "3.6.1.1"})) == ["3.6.1.1"]
    assert notes_ecs(reaction(notes={"EC_Number": "1.1.1.1 or 1.1.1.2, 1.1.1.-"})) == ["1.1.1.1", "1.1.1.2", "1.1.1.-"]
    assert notes_ecs(reaction(notes={"PROTEIN_CLASS": "2.7.1.40"})) == ["2.7.1.40"]
    assert notes_ecs(reaction(notes={"SUBSYSTEM": "1.1.1.1"})) == []


def test_xrefs():
    r = reaction(
        {
            "bigg.reaction": "PGI",
            "kegg.reaction": ["R00771", "rn:R00771", "C00001"],
            "seed": ["rxn00558", "cpd00001"],
            "seed.reaction": "rxn00559",
            "rhea": ["RHEA:11816", "11817"],
            "biocyc": "META:PGLUCISOM-RXN",
            "metacyc.reaction": "PGLUCISOM-RXN",
            "sbo": "SBO:0000176",
        }
    )
    assert xrefs(r) == {
        "bigg.reaction": ["PGI"],
        "kegg.reaction": ["R00771"],
        "seed.reaction": ["rxn00559", "rxn00558"],
        "rhea": ["11816", "11817"],
        "metacyc.reaction": ["PGLUCISOM-RXN"],
    }


def test_normalize_xref():
    assert normalize_xref("bigg.reaction", "bigg.reaction:PGI") == "PGI"
    assert normalize_xref("kegg.reaction", "kegg.reaction:R00771") == "R00771"
    assert normalize_xref("metacyc.reaction", "ECOLI:PGLUCISOM-RXN") == "PGLUCISOM-RXN"
    assert normalize_xref("rhea", "rhea:10000") == "10000"
    assert normalize_xref("rhea", "abc") is None
