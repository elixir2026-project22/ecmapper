import pytest
from cobra import Metabolite, Model, Reaction

from conftest import FIXTURES
from ecmapper import ECMapper
from ecmapper.brenda_map import (
    BrendaMap,
    MapMatcher,
    highlight_svg,
    map_matches_to_dataframe,
    map_summary,
    name_key,
    node_coverage,
)
from test_references import pyk


@pytest.fixture
def brenda_map():
    return BrendaMap.from_json(FIXTURES / "brenda_map.json")


@pytest.fixture
def matcher(sources, brenda_map):
    return MapMatcher(brenda_map, ECMapper(sources.values()))


def met(mid, name, kegg=None, compartment="c"):
    m = Metabolite(mid, name=name, compartment=compartment)
    if kegg:
        m.annotation = {"kegg.compound": kegg}
    return m


@pytest.fixture
def model():
    model = Model("m")
    r = pyk()
    r.annotation = {"kegg.reaction": "R00200", "biocyc": "META:PEPSYNTH-RXN"}
    pyr, pep, atp, amp, pi, h2o = (met(f"{m}_c", m) for m in ("pyr", "pep", "atp", "amp", "pi", "h2o"))
    succ, fum = met("succ_c", "Succinate", "C00042"), met("fum_c", "Fumarate", "C00122")
    q8, q8h2 = met("q8_c", "Ubiquinone-8", "C17569"), met("q8h2_c", "Ubiquinol-8")
    arg_c, agm_c, co2 = met("arg_c", "L-Arginine", "C00062"), met("agm_c", "Agmatine", "C00179"), met("co2_c", "CO2", "C00011")
    arg_p, agm_p = met("arg_p", "L-Arginine", "C00062", "p"), met("agm_p", "Agmatine", "C00179", "p")
    g6p, f6p = met("g6p_c", "D-Glucose 6-phosphate"), met("f6p_c", "D-Fructose 6-phosphate")
    reactions = [r]
    for rid, annotation, stoichiometry in [
        ("PPS", {"biocyc": "META:PEPSYNTH-RXN"}, {pyr: -1, atp: -1, h2o: -1, pep: 1, amp: 1, pi: 1}),
        ("SUCDi", {"ec-code": "1.3.5.1"}, {succ: -1, q8: -1, fum: 1, q8h2: 1}),
        ("PGI", {"ec-code": "5.3.1.9"}, {g6p: -1, f6p: 1}),
        ("ARGDC", {}, {arg_c: -1, agm_c: 1, co2: 1}),
        ("ARGAGMt", {}, {arg_p: -1, agm_c: -1, arg_c: 1, agm_p: 1}),
        ("NOTHING", {}, {succ: -1, pyr: 1}),
    ]:
        x = Reaction(rid)
        x.add_metabolites(stoichiometry)
        x.annotation = annotation
        reactions.append(x)
    model.add_reactions(reactions)
    return model


def test_name_key():
    assert name_key("fructose-1,6-bisphosphate") == name_key("D-Fructose 1,6-diphosphate")
    assert name_key("(S)-malate") == name_key("L-Malate") == name_key("malate")
    assert name_key("Fumaric acid") == name_key("fumarate")
    assert name_key("all-trans-retinol") == name_key("retinol")
    assert name_key("NAD+") != name_key("NADH")
    assert name_key("H+[side 1]") == "h+"


def test_expand_names(matcher):
    from ecmapper.brenda_map import expand_names

    assert expand_names(["NAD(P)H", "ATP"]) == {"NAD(P)H", "NADH", "NADPH", "ATP"}
    matcher.node_keys("node1")  # builds the indexes
    assert {"name:nad+", "name:nadp+"} <= matcher.ligand_keys(["NAD(P)+"])


def test_read_map(brenda_map):
    assert len(brenda_map) == 5
    node = brenda_map.nodes["node2"]
    assert node.brenda_reaction_ids == ["BR48097"]
    assert [(m.names, m.role) for m in node.metabolites][2:4] == [(("ATP",), "consumed"), (("H2O",), "consumed")]
    assert brenda_map.nodes["node1"].metabolites[0].names == ("phosphoenolpyruvate", "PEP")
    assert (brenda_map.nodes["node1"].x, brenda_map.nodes["node1"].y) == (100, 200)


def test_match_model(matcher, model):
    matches = matcher.match_model(model)
    status = {rid: m.status for rid, m in matches.items()}
    assert status == {
        "PYK": "on map",
        "PPS": "on map",
        "SUCDi": "on map",
        "PGI": "on map",
        "ARGDC": "on map",
        "ARGAGMt": "not on map",
        "NOTHING": "not on map",
    }
    # pyruvate kinase: by EC number; the PEP synthase node is a candidate (common EC via BKMS) but does not match
    pyk_match = matches["PYK"]
    assert [(c.node.node_id, c.score, c.via) for c in pyk_match.nodes] == [("node1", 1.0, ("ec", "metabolites"))]
    assert [c.node.node_id for c in pyk_match.candidates] == ["node1", "node2"]
    # PEP synthase: by the BRENDA reaction of the node, aligned in BKMS with PEPSYNTH-RXN
    assert "reaction-id" in matches["PPS"].nodes[0].via
    # generic quinone and quinol of the map
    assert matches["SUCDi"].nodes[0].node.node_id == "node3"
    # only by metabolites
    assert matches["ARGDC"].nodes[0].via == ("metabolites",)


def test_max_gap(matcher, model):
    # PPS matches node2 (1.0); the pyruvate kinase node (0.5) is no match anyway, but with low thresholds
    # it is only accepted if the gap to the best node is allowed
    matcher.min_score = 0.4
    m = matcher.match_reaction(model.reactions.PPS)
    assert [c.node.node_id for c in m.nodes] == ["node2"]
    matcher.max_gap = 1.0
    m = matcher.match_reaction(model.reactions.PPS)
    assert [c.node.node_id for c in m.nodes] == ["node2", "node1"]


def test_inferred_compounds(matcher, model):
    # BR240 is aligned with KEGG R02740 (C00668, C05345): "fructose-6P" is beta-D-fructose 6-phosphate
    matches = matcher.match_model(model, ["PGI"])
    assert matcher._index["inferred"] == {13: "C05345"}
    assert matches["PGI"].nodes[0].score == 1.0


def test_related(matcher, model):
    model.reactions.PGI.annotation = {"ec-code": "5.3.1.9"}
    model.reactions.PGI.add_metabolites({model.metabolites.co2_c: 1})
    m = matcher.match_reaction(model.reactions.PGI)
    assert m.status == "on map"  # 2 of 3 metabolites: 2 / 2.5
    model.reactions.PGI.add_metabolites({model.metabolites.arg_c: 1, model.metabolites.agm_c: 1})
    m = matcher.match_reaction(model.reactions.PGI)
    assert m.status == "related"
    assert not m.nodes
    assert [c.via for c in m.candidates if c.node.node_id == "node4"] == [("ec", "metabolites")]
    assert map_matches_to_dataframe({"PGI": m}).node_id.tolist() == ["node4"]


def test_outputs(matcher, model, brenda_map, tmp_path):
    matches = matcher.match_model(model)
    table = map_matches_to_dataframe(matches)
    assert list(table.columns) == [
        "reaction_id", "reaction_name", "status", "node_id", "subsystem", "ec_label", "enzyme_name",
        "score", "via", "accepted", "x", "y",
    ]
    assert table[table.reaction_id == "PYK"].node_id.tolist() == ["node1"]
    assert table[table.reaction_id == "NOTHING"].node_id.tolist() == [""]
    assert len(map_matches_to_dataframe(matches, all_candidates=True).query("reaction_id == 'PYK'")) == 2

    coverage = node_coverage(brenda_map, matches).set_index("node_id")
    assert coverage.loc["node1", "model_reactions"] == "PYK"
    assert coverage.loc["node2", "model_reactions"] == "PPS"

    summary = map_summary(brenda_map, matches)
    assert summary["status"] == {"on map": 5, "related": 0, "not on map": 2, "boundary": 0}
    assert summary["nodes_with_model_reactions"] == 5
    assert summary["nodes_per_subsystem"]["glycolysis"] == 2

    output = tmp_path / "map.svg"
    assert highlight_svg(FIXTURES / "brenda_map.svg", output, matches) == 5
    svg = output.read_text()
    assert 'class="node nodeEN pw_glycolysis ecm-on-map"' in svg
    assert "<title>pyruvate kinase | model: PYK</title>" in svg
    assert "g.nodeEN.ecm-on-map rect" in svg
    # metabolite nodes are not changed
    assert '<g id="node99" class="node nodeHM pw_glycolysis lid2">' in svg


def test_boundary(matcher, model):
    model.add_boundary(model.metabolites.co2_c, type="sink")
    assert matcher.match_model(model)["SK_co2_c"].status == "boundary"


def test_resources(matcher):
    names = [r.name for r in matcher.resources]
    assert names[-1] == "kegg-compound-names"
    assert "bkms" in names and "biggr-metabolites" in names


def test_cli(model, cache_dir, tmp_path, capsys):
    import cobra
    import pandas as pd

    from ecmapper.cli import main

    sbml = tmp_path / "model.xml"
    cobra.io.write_sbml_model(model, str(sbml))
    out, nodes, svg = tmp_path / "matches.tsv", tmp_path / "nodes.tsv", tmp_path / "map.svg"
    args = [str(sbml), "--brenda-map", str(FIXTURES / "brenda_map.json"), "--cache-dir", str(cache_dir),
            "-o", str(out), "--node-table", str(nodes), "--svg", str(FIXTURES / "brenda_map.svg"), "--svg-out", str(svg)]
    assert main(args) == 0
    printed = capsys.readouterr().out
    assert "on map      5" in printed and "map nodes with model reactions: 5 of 5" in printed
    table = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert set(table[table.status == "on map"].reaction_id) == {"PYK", "PPS", "SUCDi", "PGI", "ARGDC"}
    assert pd.read_csv(nodes, sep="\t", keep_default_na=False).model_reactions.tolist()[0] == "PYK"
    assert "ecm-on-map" in svg.read_text()
