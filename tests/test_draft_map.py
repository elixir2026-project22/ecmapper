import json
import math
import pickle
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("numba")

from ecmapper import ECMapper
from ecmapper.brenda_map import BrendaMap, MapMatcher
from ecmapper.draft_map import core, fast, relayout, run, split
from ecmapper.draft_map.core import DraftMap, Edge, MetLink, Node, Subsystem, localise_edges
from ecmapper.draft_map.layout import Ports, SequencePair, compact_axis, count_crossings, evaluate_arrangement
from ecmapper.draft_map.prepare import prepare
from test_brenda_map import model  # noqa: F401 (fixture)

BRENDA_MAP = Path(__file__).parents[1] / "BRENDAmap"


def no_overlaps(x, y, w, h, gap=0.0):
    n = len(x)
    for i in range(n):
        for j in range(i + 1, n):
            apart = (
                x[i] + w[i] + gap <= x[j] + 1e-6
                or x[j] + w[j] + gap <= x[i] + 1e-6
                or y[i] + h[i] + gap <= y[j] + 1e-6
                or y[j] + h[j] + gap <= y[i] + 1e-6
            )
            if not apart:
                return False
    return True


def test_compact_axis():
    # three boxes in a row with large gaps, one box below them
    x, w = np.array([0.0, 300.0, 600.0, 300.0]), np.full(4, 100.0)
    y, h = np.array([0.0, 0.0, 0.0, 300.0]), np.full(4, 50.0)
    cx = compact_axis(x, w, y, h, gap=20.0, ogap=0.0)
    assert cx[0] < cx[1] < cx[2]
    assert np.all(np.diff(cx[:3]) >= 100.0 + 20.0 - 1e-6)
    assert (cx - w / 2).min() == pytest.approx(0.0)
    assert cx[2] - cx[0] < 600.0
    # overlapping boxes are separated, in their order
    cx = compact_axis([0.0, 10.0], [100.0, 100.0], [0.0, 0.0], [50.0, 50.0], gap=20.0, ogap=0.0)
    assert cx[1] - cx[0] == pytest.approx(120.0)


def test_sequence_pair_packing():
    rng = np.random.default_rng(1)
    n, gap = 12, 5.0
    w, h = rng.uniform(10, 100, n), rng.uniform(10, 100, n)
    packing = SequencePair(w, h, gap)
    for _ in range(5):
        gp, gn = rng.permutation(n), rng.permutation(n)
        x, y = packing.pack(gp, gn)
        assert no_overlaps(x, y, w, h, gap)
        # the numba kernel of the annealing packs in the same way
        fx, fy = fast.pack(gp.astype(np.int64), gn.astype(np.int64), packing.w, packing.h)
        assert np.allclose(fx, x) and np.allclose(fy, y)


def test_crossings():
    p1 = np.array([[0.0, 0.0], [0.0, 10.0], [20.0, 0.0]])
    p2 = np.array([[10.0, 10.0], [10.0, 0.0], [20.0, 10.0]])
    assert count_crossings(p1, p2) == 1
    # links that share an end point do not cross
    assert count_crossings(np.array([[0.0, 0.0], [0.0, 0.0]]), np.array([[10.0, 10.0], [10.0, 0.0]])) == 0


def test_metrics_kernel():
    """The numba metrics of the annealing agree with :func:`evaluate_arrangement`."""
    rng = np.random.default_rng(2)
    n, links = 10, 30
    w, h = rng.uniform(100, 400, n), rng.uniform(100, 400, n)
    a = rng.integers(0, n, links)
    b = (a + rng.integers(1, n, links)) % n
    pa = np.stack([rng.uniform(0, w[a]), rng.uniform(0, h[a])], axis=1)
    pb = np.stack([rng.uniform(0, w[b]), rng.uniform(0, h[b])], axis=1)
    ports = Ports([f"s{i}" for i in range(n)], w, h, a, b, pa, pb, [], h)
    x, y = SequencePair(w, h, 50.0).pack(rng.permutation(n), rng.permutation(n))
    expected = evaluate_arrangement(ports, x, y)
    flips = np.zeros(n, np.bool_)
    crossings, box_crossings, length = fast.metrics(x, y, w, h, a, b, pa[:, 0].copy(), pa[:, 1].copy(),
                                                    pb[:, 0].copy(), pb[:, 1].copy(), flips, flips, h, 10.0)
    assert crossings == expected["crossings"]
    assert box_crossings == expected["box_crossings"]
    assert length == pytest.approx(expected["link_length"])


def test_localise_edges():
    dm = DraftMap()
    for args in [
        ("e1", "enzyme", "A", 0, 0),
        ("gly_a", "main", "A", 100, 0, 72),
        ("gly_b", "main", "B", 1000, 0, 72),
        ("ser_b", "main", "B", 1000, 100, 65),
        ("thr_b", "main", "B", 1000, 200, 82),
        ("e2", "enzyme", "B", 1100, 100),
    ]:
        dm.nodes[args[0]] = Node(*args[:5], 60, 30, args[0], ligand_id=args[5] if len(args) > 5 else None)
    dm.edges = [Edge("e1", m, "main", False, True) for m in ("gly_a", "gly_b", "ser_b", "thr_b")]
    dm.edges += [Edge("e2", "gly_b", "main", False, True), Edge("e2", "ser_b", "main", False, True)]
    assert localise_edges(dm) == 3
    by_enzyme = {(e.enzyme, e.met) for e in dm.edges}
    # glycine is drawn in A already: the edge to B is removed; serine and threonine get copies in A
    assert by_enzyme == {("e1", "gly_a"), ("e1", "ser_b_in_A"), ("e1", "thr_b_in_A"), ("e2", "gly_b"), ("e2", "ser_b")}
    assert dm.nodes["ser_b_in_A"].subsystem == "A" and dm.nodes["ser_b_in_A"].ligand_id == 65
    # threonine is not used in B any more
    assert "thr_b" not in dm.nodes
    assert {(l.a, l.b) for l in dm.links} == {("gly_a", "gly_b"), ("ser_b_in_A", "ser_b")}


def prepared(sources, model, tmp_path):  # noqa: F811
    brenda_map = BrendaMap.from_json(BRENDA_MAP / "BRENDA_map_reactions.json")
    matcher = MapMatcher(brenda_map, ECMapper(sources.values()), vicinity="propagate")
    data = prepare(model, matcher, {"PGI": "Glycolysis/Gluconeogenesis"})
    prep = tmp_path / "prep.pkl"
    with open(prep, "wb") as f:
        pickle.dump(data, f)
    return data, prep


def run_draft_map(prep, out, *options):
    run.main([str(prep), str(BRENDA_MAP / "BRENDA_map_reactions.json"), str(BRENDA_MAP / "BRENDA_map.svg"),
              str(out), "--iterations", "300", *options])


def test_draft_map(sources, model, tmp_path):  # noqa: F811
    """The whole pipeline on the BRENDA map with a small model."""
    data, prep = prepared(sources, model, tmp_path)
    assert data["reactions"]["PGI"]["subsystem"] == "Glycolysis/Gluconeogenesis"
    assert data["reactions"]["ARGAGMt"]["transport"]
    on_map = {rid for rid, m in data["matches"].items() if m["status"] == "on map"}
    assert {"PYK", "PGI"} <= on_map
    out = tmp_path / "out"
    # all subsystems with model reactions (rules A-D in test_draft_map_redundant)
    run_draft_map(prep, out, "--keep-subsets", "--aside-coverage", "2")

    stats = json.loads((out / "stats.json").read_text())
    nodes = {
        c["node_id"] for m in data["matches"].values() for c in m["candidates"] if c["accepted"]
    }
    assert stats["step2_pruned"]["enzyme_nodes"] == len(nodes)
    assert 0 < stats["step1_subsystems"]["subsystems"] < stats["original"]["subsystems"]
    positions = json.loads((out / "reaction_positions.json").read_text())
    # NOTHING (succinate -> pyruvate) is not on the map, but both metabolites are drawn: added in step 3
    assert set(positions) == on_map | {"NOTHING"}
    assert stats["step3_completed"]["enzyme_nodes_added"] == 1
    # step 3b: the subsystem of the added reaction was laid out again (with what is installed)
    assert len(stats["relayout"]) == 1 and stats["relayout"][0]["layout"] in relayout.METHODS
    assert (out / "step3b_relayout.svg").is_file()
    assert set(json.loads((out / "draft_map.maps").read_text())) == set(positions)
    assert (out / "draft_map.svg").stat().st_size > 0

    # Escher map: all reactions on the map, segments between existing nodes
    header, body = json.loads((out / "draft_map.escher.json").read_text())
    assert header["schema"] == "https://escher.github.io/escher/jsonschema/1-0-0#"
    assert {r["bigg_id"] for r in body["reactions"].values()} == set(positions)
    for reaction in body["reactions"].values():
        for segment in reaction["segments"].values():
            assert segment["from_node_id"] in body["nodes"] and segment["to_node_id"] in body["nodes"]
        assert {m["bigg_id"] for m in reaction["metabolites"]} == set(data["reactions"][reaction["bigg_id"]]["metabolites"])
    metabolites = {n["bigg_id"] for n in body["nodes"].values() if n["node_type"] == "metabolite"}
    assert metabolites <= set(data["metabolites"])



def test_draft_map_packed(sources, model, tmp_path):  # noqa: F811
    """--pack: no inter-pathway links, the boxes packed without overlaps."""
    data, prep = prepared(sources, model, tmp_path)
    out = tmp_path / "out"
    run_draft_map(prep, out, "--pack")
    stats = json.loads((out / "stats.json").read_text())
    assert stats["step3_completed"]["inter_pathway_links"] > 0
    assert stats["final"]["inter_pathway_links"] == 0
    assert stats["arrangements"]["packed"]["area_ratio"] < stats["arrangements"]["original order (compacted)"]["area_ratio"]
    assert 'class="link inter-pw-link' not in (out / "draft_map.svg").read_text()
    with open(out / "draft_map.pkl", "rb") as f:
        dm = pickle.load(f)
    boxes = list(dm.subsystems.values())
    assert no_overlaps([b.x for b in boxes], [b.y for b in boxes], [b.w for b in boxes], [b.h for b in boxes], 80.0)
    assert all(dm.nodes[e.enzyme].subsystem == dm.nodes[e.met].subsystem for e in dm.edges)
    with pytest.raises(SystemExit):
        run_draft_map(prep, out, "--pack", "--hub-subsystems", "6")


def test_split_map():
    """Links between maps become references on both ends, links inside a map stay."""
    dm = DraftMap()
    for k, group in (("glycolysis", "pw1"), ("serine-metabolism", "pw3"), ("threonine-metabolism", "pw3"),
                     ("lipid-metabolism", "pw2")):
        dm.subsystems[k] = Subsystem(k, k, group, "#fff", 0, 0, 100, 100)
    for i, s in enumerate(dm.subsystems):
        dm.nodes[f"e{i}"] = Node(f"e{i}", "enzyme", s, 0, 0, 60, 30, f"e{i}", reactions=[f"R{i}"])
        dm.nodes[f"m{i}"] = Node(f"m{i}", "main", s, 50, 0, 60, 30, f"m{i}", ligand_id=1)
        dm.edges.append(Edge(f"e{i}", f"m{i}", "main", False, True))
    dm.links = [MetLink("m0", "m1"), MetLink("m1", "m2"), MetLink("m2", "m3")]
    dm.orig_centres = {k: (0, 0) for k in dm.subsystems}
    order, assignment = split.assign_maps(dm, ["glycolysis"], min_reactions=0)
    assert order == [split.CENTRAL, "Lipids", "Amino acids"]
    assert assignment["serine-metabolism"] == assignment["threonine-metabolism"] == "Amino acids"
    maps = dict(split.split_map(dm, order, assignment))
    assert [(l.a, l.b) for l in maps["Amino acids"].links] == [("m1", "m2")]
    assert maps[split.CENTRAL].nodes["m0"].stubs == ["Amino acids"]
    assert maps["Amino acids"].nodes["m1"].stubs == [split.CENTRAL]
    assert maps["Amino acids"].nodes["m2"].stubs == ["Lipids"]
    assert not hasattr(maps["Amino acids"].nodes["e1"], "stubs")
    # small colour group maps are merged into the map they have most links with
    order, assignment = split.assign_maps(dm, ["glycolysis"], min_reactions=2)
    assert order == [split.CENTRAL, "Amino acids"] and assignment["lipid-metabolism"] == "Amino acids"
    # explicit maps come first
    order, assignment = split.assign_maps(dm, ["glycolysis"], {"Serine": ["serine-metabolism"]}, min_reactions=0)
    assert order[0] == "Serine" and assignment["serine-metabolism"] == "Serine"


def test_draft_map_split(sources, model, tmp_path):  # noqa: F811
    """--split: every drawn reaction is on a map, the Escher maps are complete, the overview links the maps."""
    _, prep = prepared(sources, model, tmp_path)
    whole, out = tmp_path / "whole", tmp_path / "out"
    run_draft_map(prep, whole, "--pack")
    run_draft_map(prep, out, "--split", "--min-map-reactions", "0")
    index = json.loads((out / "maps.json").read_text())
    assert index["maps"][0]["name"] == split.CENTRAL and index["maps"][-1]["name"] == split.ASIDE
    assert len(index["maps"]) > 1
    drawn = set(json.loads((whole / "reaction_positions.json").read_text()))
    drawn |= set(json.loads((whole / "set_aside.reaction_positions.json").read_text()))
    assert set().union(*(m["reactions"] for m in index["maps"])) == drawn
    references = 0
    for entry in index["maps"]:
        prefix = entry["file"].removesuffix(".escher.json")
        _, body = json.loads((out / entry["file"]).read_text())
        assert {r["bigg_id"] for r in body["reactions"].values()} == set(entry["reactions"])
        assert set(json.loads((out / f"{prefix}.maps").read_text())) == set(entry["reactions"])
        references += sum(t["text"].startswith("→ ") for t in body["text_labels"].values())
        with open(out / f"{prefix}.pkl", "rb") as f:
            m = pickle.load(f)
        boxes = list(m.subsystems.values())
        assert no_overlaps([b.x for b in boxes], [b.y for b in boxes], [b.w for b in boxes], [b.h for b in boxes], 80.0)
    stats = json.loads((out / "stats.json").read_text())
    # one label per drawn node with references that is used in the Escher map
    references_total = sum(s["references"] for s in stats["maps"].values())
    assert 0 < references <= references_total
    overview = (out / "overview.svg").read_text()
    assert all(f'href="{m["file"].removesuffix(".escher.json")}.svg"' in overview for m in index["maps"])
    with pytest.raises(SystemExit):
        run_draft_map(prep, out, "--central", "glycolysis")
    # --split --pack: no links on the maps, the references stay
    packed = tmp_path / "packed"
    run_draft_map(prep, packed, "--split", "--pack", "--min-map-reactions", "0")
    stats = json.loads((packed / "stats.json").read_text())
    assert all(s["final"]["inter_pathway_links"] == 0 for s in stats["maps"].values())
    assert sum(s["references"] for s in stats["maps"].values()) == references_total


def test_drop_redundant_subsystems():
    """A: proper subsets are removed; B, C: (mostly) covered subsystems and D: few matched nodes are set aside."""
    dm = DraftMap()
    content = {
        "big": ["R1", "R2", "R3"],
        "subset": ["R1", "R2"],  # A: proper subset of big
        "equal1": ["R30"],  # equal sets: the first (by name) is set aside (B), the other stays
        "equal2": ["R30"],
        "mixed": ["R1", "R50", *(f"R{i}" for i in range(5, 14))],  # 10 of 11 in big and other (C)
        "other": [*(f"R{i}" for i in range(5, 14)), "R40", "R41"],
        "spread": ["R3", "R41", "R60"],  # 2 of 3 elsewhere: stays with coverage 0.9
        "sparse": ["R20"],  # one of five enzyme nodes matched (D)
    }
    for k, rs in content.items():
        dm.subsystems[k] = Subsystem(k, k, "pw1", "#fff", 0, 0, 100, 100)
        for r in rs:
            dm.nodes[f"{k}_{r}"] = Node(f"{k}_{r}", "enzyme", k, 0, 0, 60, 30, r, reactions=[r])
    for i in range(4):
        dm.nodes[f"sparse_{i}"] = Node(f"sparse_{i}", "enzyme", "sparse", 0, 0, 60, 30, "")
    result = core.drop_redundant_subsystems(dm, coverage=0.9, min_matched=0.25)
    assert [d["subsystem"] for d in result["dropped"]] == ["subset"]
    assert "subset" not in dm.subsystems and not dm.members("subset")
    assert {a["subsystem"]: a["rule"] for a in result["set_aside"]} == {"equal1": "B", "mixed": "C", "sparse": "D"}
    assert {k for k, s in dm.subsystems.items() if not s.aside} == {"big", "equal2", "spread", "other"}
    # no rules
    dm2 = DraftMap()
    dm2.subsystems = {k: Subsystem(k, k, "pw1", "#fff", 0, 0, 100, 100) for k in ("a", "b")}
    dm2.nodes = {k: Node(k, "enzyme", k, 0, 0, 60, 30, "", reactions=["R1"]) for k in ("a", "b")}
    assert core.drop_redundant_subsystems(dm2, drop_subsets=False, coverage=2.0) == dict(dropped=[], set_aside=[])


def test_draft_map_set_aside(sources, model, tmp_path):  # noqa: F811
    """Set-aside subsystems are drawn on a separate map, not on the draft map."""
    _, prep = prepared(sources, model, tmp_path)
    out = tmp_path / "out"
    run_draft_map(prep, out, "--aside-coverage", "1.0")
    stats = json.loads((out / "stats.json").read_text())
    aside = {a["subsystem"] for a in stats["redundant_subsystems"]["set_aside"]}
    assert aside
    with open(out / "draft_map.pkl", "rb") as f:
        main = pickle.load(f)
    with open(out / "set_aside.pkl", "rb") as f:
        side = pickle.load(f)
    assert set(side.subsystems) == aside and not aside & set(main.subsystems)
    assert not {d["subsystem"] for d in stats["redundant_subsystems"]["dropped"]} & (set(main.subsystems) | aside)
    # every reaction of a set-aside subsystem is also on the draft map (rule B)
    assert {r for n in side.nodes.values() for r in n.reactions} <= {r for n in main.nodes.values() for r in n.reactions}
    _, body = json.loads((out / "set_aside.escher.json").read_text())
    assert {r["bigg_id"] for r in body["reactions"].values()} == {r for n in side.nodes.values() for r in n.reactions}


def relayout_map():
    """A subsystem with two map enzymes and an added reaction whose edges cross them, and one without additions."""
    dm = DraftMap()
    for k in ("S", "T"):
        dm.subsystems[k] = Subsystem(k, k, "pw1", "#fff", 0, 0, 1000, 600)
    for i, x, y in (("m1", 0, 0), ("m2", 400, 0), ("m3", 0, 200), ("m4", 400, 200)):
        dm.nodes[i] = Node(i, "main", "S", x, y, 60, 30, i)
    for i, x, y in (("e1", 200, 0), ("e2", 200, 200)):
        dm.nodes[i] = Node(i, "enzyme", "S", x, y, 60, 30, i, reactions=[i.upper()], status="on-map")
    dm.nodes["a1"] = Node("a1", "enzyme", "S", 600, 100, 60, 30, "A1", reactions=["A1"], status="added")
    dm.nodes["n1"] = Node("n1", "main", "S", -300, 100, 60, 30, "n1", status="new")
    dm.nodes["c1"] = Node("c1", "cofactor", "S", 640, 140, 40, 20, "ATP", status="added")
    dm.nodes["t1"] = Node("t1", "enzyme", "T", 0, 0, 60, 30, "t1", reactions=["T1"], status="on-map")
    dm.nodes["t2"] = Node("t2", "main", "T", 100, 0, 60, 30, "t2")
    dm.edges = [Edge("e1", "m1", "main", True, False), Edge("e1", "m2", "main", False, True),
                Edge("e2", "m3", "main", True, False), Edge("e2", "m4", "main", False, True),
                Edge("a1", "n1", "main", True, False), Edge("a1", "m2", "main", False, True),
                Edge("a1", "m3", "main", False, True), Edge("a1", "c1", "cofactor", True, False),
                Edge("t1", "t2", "main", False, True)]
    return dm


def test_relayout_spring_keeps_map_nodes():
    pytest.importorskip("networkx")
    dm = relayout_map()
    s = relayout.subsystem_map(dm, "S")
    relayout.spring_layout(s)
    for i in ("m1", "m2", "m3", "m4", "e1", "e2"):
        assert (s.nodes[i].x, s.nodes[i].y) == (dm.nodes[i].x, dm.nodes[i].y)
    assert (s.nodes["a1"].x, s.nodes["a1"].y) != (dm.nodes["a1"].x, dm.nodes["a1"].y)
    # the cofactor moves with its enzyme node
    assert s.nodes["c1"].x - s.nodes["a1"].x == pytest.approx(40) and s.nodes["c1"].y - s.nodes["a1"].y == pytest.approx(40)


def test_relayout_graphviz():
    pytest.importorskip("pygraphviz")
    dm = relayout_map()
    pinned = ("m1", "m2", "m3", "m4", "e1", "e2")
    for prog in ("fdp", "neato", "dot"):
        s = relayout.subsystem_map(dm, "S")
        relayout.LAYOUTS[prog](s)
        assert all(math.isfinite(n.x) and math.isfinite(n.y) for n in s.nodes.values())
        if prog == "fdp":
            # pinned nodes keep their relative positions (Graphviz translates the whole drawing)
            dx, dy = s.nodes["m1"].x - dm.nodes["m1"].x, s.nodes["m1"].y - dm.nodes["m1"].y
            for i in pinned:
                assert s.nodes[i].x - dm.nodes[i].x == pytest.approx(dx, abs=0.2)
                assert s.nodes[i].y - dm.nodes[i].y == pytest.approx(dy, abs=0.2)
        if prog == "dot":
            # layered along the reaction directions: substrate above enzyme above product
            assert s.nodes["m1"].y < s.nodes["e1"].y < s.nodes["m2"].y


def test_relayout_subsystems():
    pytest.importorskip("networkx")
    pytest.importorskip("pygraphviz")
    dm = relayout_map()
    before = {i: (n.x, n.y) for i, n in dm.nodes.items()}
    rows = relayout.relayout_subsystems(dm)
    assert [r["subsystem"] for r in rows] == ["S"]  # T has no added reactions
    assert {i: (n.x, n.y) for i, n in dm.nodes.items() if n.subsystem == "T"} == {"t1": before["t1"], "t2": before["t2"]}
    row = rows[0]
    assert set(row) >= {"current", "spring", "fdp", "dot", "neato", "layout"}
    best_keeping = min(row[m]["crossings_main"] for m in relayout.KEEPING)
    chosen = row[row["layout"]]["crossings_main"]
    if row["layout"] in relayout.REDRAWING:
        assert chosen < 0.75 * best_keeping
    else:
        assert chosen == best_keeping
    assert chosen <= row["current"]["crossings_main"]
    # no layouts: nothing changes
    dm = relayout_map()
    rows = relayout.relayout_subsystems(dm, methods=[])
    assert rows[0]["layout"] == "current" and {i: (n.x, n.y) for i, n in dm.nodes.items()} == before
