from cobra import Metabolite, Model, Reaction

from ecmapper import ECMapper, evaluate, evaluation_summary
from test_references import pyk


def test_evaluate(sources):
    model = Model("m")
    r = pyk()
    r.annotation = {"ec-code": "2.7.1.40", "kegg.reaction": "R00200", "biocyc": "META:PEPSYNTH-RXN", "bigg.reaction": "PYK"}
    model.add_reactions([r])
    a, b = Metabolite("a_c", compartment="c"), Metabolite("b_c", compartment="c")
    for rid, annotation in [
        ("GLYCK2", {"ec-code": "2.7.1.31", "kegg.reaction": "R08572"}),  # KEGG: 2.7.1.165 only
        ("TRANSFERRED", {"ec-code": "4.2.1.52"}),  # nothing found
        ("NO_EC", {"kegg.reaction": "R00200"}),  # no reference, not evaluated
    ]:
        x = Reaction(rid)
        x.add_metabolites({a: -1, b: 1})
        x.annotation = annotation
        model.add_reactions([x])
    mapper = ECMapper(sources.values())

    table = evaluate(model, mapper).set_index("reaction_id")
    assert list(table.index) == ["PYK", "GLYCK2", "TRANSFERRED"]
    assert table.loc["PYK", ["selection", "selected", "outcome"]].tolist() == ["metabolites", "2.7.1.40", "consistent"]
    assert table.loc["GLYCK2", ["selection", "selected", "outcome"]].tolist() == ["single", "2.7.1.165", "disjoint"]
    # the reference is made current
    assert table.loc["TRANSFERRED", ["reference", "outcome"]].tolist() == ["4.3.3.7", "not found"]
    # the model is not changed
    assert model.reactions.PYK.annotation["ec-code"] == "2.7.1.40"

    summary = evaluation_summary(evaluate(model, mapper))
    assert summary.loc["metabolites", "consistent"] == 1
    assert summary.loc["total", "total"] == 3

    # BiGG IDs only: PYK has no BiGGr entry in the fixtures
    table = evaluate(model, mapper, bigg_only=True).set_index("reaction_id")
    assert table.loc["PYK", "outcome"] == "not found"
