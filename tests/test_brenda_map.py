import pytest
from cobra import Metabolite, Model, Reaction

from conftest import FIXTURES
from ecmapper import ECMapper
from ecmapper.brenda_map import (
    BrendaMap,
    MapMatcher,
    MapMetabolite,
    MapNode,
    highlight_svg,
    map_matches_to_dataframe,
    map_summary,
    name_key,
    name_keys,
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
    # alpha-/beta- of anomers are ignored, others are part of the name
    assert name_key("alpha-D-glucose 6-phosphate") == name_key("D-Glucose 6-phosphate")
    assert name_key("beta-alanine") != name_key("L-alanine")
    assert name_key("beta-carotene") != name_key("alpha-carotene")
    assert name_key("alpha-ketoglutarate") == "alphaketoglutarate"
    assert name_key("UDP-alpha-N-acetyl-D-glucosamine") == name_key("UDP-N-acetyl-alpha-D-glucosamine")


def test_name_keys():
    assert name_keys("alpha-D-Glucose") == {"glucose", "alphadglucose"}
    assert name_keys("L-Alanine") == {"alanine", "lalanine"}
    assert name_keys("beta-Alanine") == {"betaalanine"}
    assert name_keys("Pyruvate") == {"pyruvate"}


def test_score_one_to_one():
    g6p, g6p_alpha, inositol_p = frozenset({"a", "b"}), frozenset({"a", "c"}), frozenset({"d"})
    # an anomerase: both anomers match the one glucose 6-phosphate of the node, but only one can
    assert MapMatcher.score([g6p, g6p_alpha], [g6p, inositol_p]) == pytest.approx(1 / 3)
    # ... and a reaction with only glucose 6-phosphate matches only one anomer of an anomerase node
    assert MapMatcher.score([g6p], [g6p, g6p_alpha]) == 0.5
    assert MapMatcher.score([g6p, g6p_alpha], [g6p_alpha, g6p]) == 1.0
    # a maximum matching: {x, y} first takes {y}, which {y} needs, and moves to {x}
    assert MapMatcher.score([frozenset("xy"), frozenset("y")], [frozenset("y"), frozenset("x")]) == 1.0


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
    # the two pyruvate nodes are linked, the PEP nodes not
    assert brenda_map.links == [("node91", "node92")]
    assert brenda_map.neighbours("node1") == [("node2", "node92")]
    assert [m.names for m in brenda_map.place_metabolites("node92")] == [("pyruvate",), ("pyruvate",)]
    assert brenda_map.neighbours("node3") == []


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
        "score", "via", "vicinity_of", "accepted", "x", "y",
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

    assert main(args + ["--map-vicinity", "propagate", "--min-map-score-vicinity", "0.5"]) == 0
    assert "only found in the vicinity of other nodes: 0 reactions, 0 nodes" in capsys.readouterr().out


# --- vicinity -------------------------------------------------------------------------------


def drawn(name, place=None, role="main"):
    return MapMetabolite((name,), None, role, place)


def glycolysis_map(second_pgi=False):
    """PGI - PFK - FBA - TPI, connected by drawn metabolites (places node21 ... node25)."""
    nodes = [
        MapNode("node11", "glycolysis", "5.3.1.9", ["5.3.1.9"], "glucose-6-phosphate isomerase", [],
                [drawn("glucose 6-phosphate", "node21"), drawn("fructose 6-phosphate", "node22")]),
        # Mg2+ is not in the model: 3 / 4
        MapNode("node12", "glycolysis", "2.7.1.11", ["2.7.1.11"], "6-phosphofructokinase", [],
                [drawn("fructose 6-phosphate", "node22"), drawn("fructose 1,6-bisphosphate", "node23"),
                 drawn("ATP", role="consumed"), drawn("ADP", role="produced"), drawn("Mg2+", role="consumed")]),
        # GAP is not recognised: 2 / 4
        MapNode("node13", "glycolysis", "4.1.2.13", ["4.1.2.13"], "fructose-bisphosphate aldolase", [],
                [drawn("fructose 1,6-bisphosphate", "node23"), drawn("dihydroxyacetone phosphate", "node24"),
                 drawn("GAP", "node25")]),
        MapNode("node14", "glycolysis", "5.3.1.1", ["5.3.1.1"], "triose-phosphate isomerase", [],
                [drawn("dihydroxyacetone phosphate", "node24"), drawn("GAP", "node25")]),
    ]
    links = []
    if second_pgi:
        nodes.append(MapNode("node15", "gluconeogenesis", "5.3.1.9", ["5.3.1.9"], "glucose-6-phosphate isomerase",
                             [], [drawn("glucose 6-phosphate", "node26"), drawn("fructose 6-phosphate", "node27")]))
        links.append(("node27", "node22"))
    return BrendaMap(nodes, links)


@pytest.fixture
def glycolysis():
    model = Model("glycolysis")
    g6p, f6p = met("g6p_c", "D-Glucose 6-phosphate"), met("f6p_c", "D-Fructose 6-phosphate")
    fdp, dhap = met("fdp_c", "D-Fructose 1,6-bisphosphate"), met("dhap_c", "Dihydroxyacetone phosphate")
    g3p, atp = met("g3p_c", "Glyceraldehyde 3-phosphate"), met("atp_c", "ATP")
    adp, h = met("adp_c", "ADP"), met("h_c", "H+")
    reactions = []
    for rid, ec, stoichiometry in [
        ("PGI", "5.3.1.9", {g6p: -1, f6p: 1}),
        ("PFK", None, {f6p: -1, atp: -1, fdp: 1, adp: 1, h: 1}),
        ("PFKX", "3.1.3.11", {f6p: -1, atp: -1, fdp: 1, adp: 1, h: 1}),  # contradicting EC number
        ("FBA", "4.1.2.13", {fdp: -1, dhap: 1, g3p: 1}),
        ("TPI", None, {dhap: -1, g3p: 1}),
    ]:
        r = Reaction(rid)
        r.add_metabolites(stoichiometry)
        r.annotation = {"ec-code": ec} if ec else {}
        reactions.append(r)
    model.add_reactions(reactions)
    return model


def statuses(matches):
    return {rid: (m.status, [c.node.node_id for c in m.nodes]) for rid, m in matches.items()}


def test_neighbours(sources):
    brenda_map = glycolysis_map(second_pgi=True)
    assert sorted(brenda_map.neighbours("node12")) == [("node11", "node22"), ("node13", "node23"), ("node15", "node22")]
    assert sorted(brenda_map.neighbours("node15")) == [("node11", "node22"), ("node12", "node22")]


def test_compatible_ecs(sources):
    matcher = MapMatcher(glycolysis_map(), ECMapper(sources.values()))
    assert matcher.compatible_ecs({"2.7.1.1"}, "node12")
    assert not matcher.compatible_ecs({"3.1.3.11", "2.7.4.-"}, "node12")
    assert matcher.compatible_ecs({"2.7.-.-"}, "node12")  # only incomplete EC numbers
    assert matcher.compatible_ecs(set(), "node12")


def test_no_vicinity(sources, glycolysis):
    matches = MapMatcher(glycolysis_map(), ECMapper(sources.values())).match_model(glycolysis)
    assert statuses(matches) == {
        "PGI": ("on map", ["node11"]),
        "PFK": ("not on map", []),  # 0.75, found only by metabolites
        "PFKX": ("not on map", []),
        "FBA": ("related", []),  # 0.5
        "TPI": ("not on map", []),
    }


def test_vicinity_neighbours(sources, glycolysis):
    matcher = MapMatcher(glycolysis_map(), ECMapper(sources.values()), vicinity="neighbours")
    matches = matcher.match_model(glycolysis)
    # only the neighbours of the node of PGI; PFKX has an EC number of another enzyme class
    assert statuses(matches) == {
        "PGI": ("on map", ["node11"]),
        "PFK": ("on map", ["node12"]),
        "PFKX": ("not on map", []),
        "FBA": ("related", []),
        "TPI": ("not on map", []),
    }
    pfk = matches["PFK"].nodes[0]
    assert (pfk.score, pfk.via, pfk.vicinity_of) == (0.75, ("metabolites", "vicinity"), ("PGI", "node11"))

    matcher.min_score_vicinity_unlinked = 0.8
    assert matcher.match_model(glycolysis)["PFK"].status == "not on map"


def test_vicinity_propagate(sources, glycolysis):
    matcher = MapMatcher(glycolysis_map(), ECMapper(sources.values()))
    matches = matcher.match_model(glycolysis)
    new = matcher.search_vicinity(glycolysis.reactions, matches, "propagate")
    assert new == [("PFK", "node12"), ("FBA", "node13")]
    fba = matches["FBA"]
    assert fba.status == "on map"
    assert (fba.nodes[0].score, fba.nodes[0].via, fba.nodes[0].vicinity_of) == (
        0.5, ("ec", "metabolites", "vicinity"), ("PFK", "node12")
    )
    # TPI matches only 1 of 3 metabolites
    assert matches["TPI"].status == "not on map"

    matcher.min_score_vicinity = 0.6
    assert matcher.search_vicinity(glycolysis.reactions, matcher.match_model(glycolysis), "propagate") == [
        ("PFK", "node12")
    ]
    with pytest.raises(ValueError):
        matcher.search_vicinity(glycolysis.reactions, matches, "everywhere")


def test_vicinity_unambiguous(sources, glycolysis):
    # PGI is on two nodes; propagation does not start from it, the search of all neighbours does
    matcher = MapMatcher(glycolysis_map(second_pgi=True), ECMapper(sources.values()))
    matches = matcher.match_model(glycolysis)
    assert [c.node.node_id for c in matches["PGI"].nodes] == ["node11", "node15"]
    assert matcher.search_vicinity(glycolysis.reactions, matches, "propagate") == []
    assert matcher.search_vicinity(glycolysis.reactions, matches, "neighbours") == [("PFK", "node12")]


def test_vicinity_outputs(sources, glycolysis, tmp_path):
    brenda_map = glycolysis_map()
    matches = MapMatcher(brenda_map, ECMapper(sources.values()), vicinity="propagate").match_model(glycolysis)
    table = map_matches_to_dataframe(matches).set_index("reaction_id")
    assert table.loc["FBA", "vicinity_of"] == "PFK@node12"
    assert table.loc["PGI", "vicinity_of"] == ""
    summary = map_summary(brenda_map, matches)
    assert (summary["vicinity_reactions"], summary["vicinity_nodes"]) == (2, 2)

    svg = tmp_path / "in.svg"
    svg.write_text(
        "<svg><style></style>"
        + "".join(f'<g id="node{i}" class="node nodeEN"><rect/><title>n{i}</title></g>' for i in range(11, 15))
        + "</svg>"
    )
    assert highlight_svg(svg, tmp_path / "out.svg", matches) == 3
    out = (tmp_path / "out.svg").read_text()
    assert 'class="node nodeEN ecm-on-map"><rect/><title>n11 | model: PGI</title>' in out
    assert 'class="node nodeEN ecm-vicinity"><rect/><title>n13 | model (vicinity): FBA</title>' in out
