import copy
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
from ecmapper.draft_map import core, escher, fast, layout, relayout, run, split
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
    # all subsystems with model reactions (rules A-E in test_drop_redundant_subsystems)
    run_draft_map(prep, out, "--keep-subsets", "--keep-borrowed", "--aside-coverage", "2")

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
    # (default: the conservative placement, NOTHING crosses nothing and stays)
    assert len(stats["relayout"]) == 1 and stats["relayout"][0]["layout"] == relayout.CONSERVATIVE
    assert stats["relayout"][0]["left_out"] == []
    assert (out / "step3b_relayout.svg").is_file()
    assert set(json.loads((out / "draft_map.maps").read_text())) == set(positions)
    assert (out / "draft_map.svg").stat().st_size > 0
    # all layouts of step 3b of that subsystem side by side, then stop
    k = stats["relayout"][0]["subsystem"]
    compared = tmp_path / "compared"
    run_draft_map(prep, compared, "--keep-subsets", "--keep-borrowed", "--aside-coverage", "2",
                  "--compare-layouts", f"{k},no-such-subsystem")
    index = (compared / "layouts" / "index.html").read_text(encoding="utf-8")
    assert f"{k}_current.svg" in index and f"{k}_{relayout.CONSERVATIVE}.svg" in index
    assert (compared / "layouts" / f"{k}_current.svg").is_file() and not (compared / "draft_map.svg").exists()
    # place_conservative is used alone
    with pytest.raises(SystemExit):
        run_draft_map(prep, compared, "--relayout", f"{relayout.CONSERVATIVE},spring")

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
    # curved branches to the secondary metabolites; main segments straight unless routed around other nodes
    curved_main = 0
    for reaction in body["reactions"].values():
        for segment in reaction["segments"].values():
            ends = [body["nodes"][segment[k]] for k in ("from_node_id", "to_node_id")]
            secondary = any(n["node_type"] == "metabolite" and not n["node_is_primary"] for n in ends)
            if secondary:
                assert segment["b1"] is not None
            elif "metabolite" not in {n["node_type"] for n in ends}:
                assert segment["b1"] is None
            else:
                curved_main += segment["b1"] is not None
    assert curved_main == 0 or stats["routed_edges"] > 0



def test_draft_map_packed(sources, model, tmp_path):  # noqa: F811
    """--pack: no inter-pathway links, the boxes packed without overlaps."""
    data, prep = prepared(sources, model, tmp_path)
    out = tmp_path / "out"
    # condensed subsystems: with --uniform-shrink (the default), glycolysis keeps the long column of the BRENDA map,
    # and the annealing gives up area for the target aspect ratio, so packing takes more area than the original order
    run_draft_map(prep, out, "--pack", "--no-uniform-shrink")
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


def test_draft_map_split_assignment(sources, model, tmp_path):  # noqa: F811
    """--write-assignment writes the maps of the subsystems and the subsystems of the reactions; --assignment reads
    the edited table: a subsystem without map is left out, a new map column is a map, reactions move."""
    pytest.importorskip("openpyxl")
    from openpyxl import load_workbook

    from ecmapper.draft_map.assignment import read_assignment

    _, prep = prepared(sources, model, tmp_path)
    reference, out = tmp_path / "reference", tmp_path / "out"
    run_draft_map(prep, reference, "--split", "--pack", "--min-map-reactions", "0")
    xlsx = tmp_path / "assignment.xlsx"
    run_draft_map(prep, out, "--split", "--min-map-reactions", "0", "--write-assignment", str(xlsx))
    assert not (out / "maps.json").exists()
    index = json.loads((reference / "maps.json").read_text())
    table = read_assignment(xlsx)
    assert table.maps == [m["name"] for m in index["maps"]]
    assert table.subsystems == {s: [m["name"]] for m in index["maps"] for s in m["subsystems"]}
    drawn = {r for m in index["maps"] for r in m["reactions"]}
    assert {r for r, subs in table.reactions.items() if subs} == drawn
    # reactions added from the model (purple), those of the BRENDA map not
    from ecmapper.draft_map.assignment import ADDED_COLOUR, REACTION_COLUMNS

    ws = load_workbook(xlsx)["Reactions"]
    i_drawn = REACTION_COLUMNS.index("Drawn as")
    fills = {row[0].value: (row[i_drawn].value, row[0].fill.fgColor.rgb) for row in ws.iter_rows(min_row=2)}
    assert any(drawn_as == "added" for drawn_as, _ in fills.values())
    for drawn_as, fill in fills.values():
        assert fill.endswith(ADDED_COLOUR) == (drawn_as == "added")

    # edit: the first subsystem of the last map without map, the first of the central map on a new map "Mine",
    # PGI only in a subsystem of the central map
    left_out, moved = index["maps"][-1]["subsystems"][0], index["maps"][0]["subsystems"][0]
    target = index["maps"][0]["subsystems"][-1]
    wb = load_workbook(xlsx)
    ws = wb["Subsystems"]
    new_col = ws.max_column + 1
    ws.cell(1, new_col, "Mine")
    for row in ws.iter_rows(min_row=2):
        if row[0].value == left_out:
            for c in row[5:]:
                c.value = None
        elif row[0].value == moved:
            for c in row[5:]:
                c.value = None
            ws.cell(row[0].row, new_col, "x")
    ws = wb["Reactions"]
    for row in ws.iter_rows(min_row=2):
        if row[0].value == "PGI":
            row[2].value = target
    wb.save(xlsx)
    edited = tmp_path / "edited"
    run_draft_map(prep, edited, "--split", "--pack", "--assignment", str(xlsx))
    index = json.loads((edited / "maps.json").read_text())
    names = [m["name"] for m in index["maps"]]
    assert names[-1] == "Mine" and {m["name"]: m for m in index["maps"]}["Mine"]["subsystems"] == [moved]
    assert all(left_out not in m["subsystems"] for m in index["maps"])
    assert [m["name"] for m in index["maps"] if "PGI" in m["reactions"]] == [split.CENTRAL]
    with open(edited / f"{split.slug(split.CENTRAL)}.pkl", "rb") as f:
        central = pickle.load(f)
    assert {n.subsystem for n in central.nodes.values() if "PGI" in n.reactions} == {target}
    with pytest.raises(SystemExit):
        run_draft_map(prep, edited, "--split", "--assignment", str(xlsx), "--central", "glycolysis")


def test_duplicate_subsystems():
    """A subsystem on two maps: a copy on the second map with the links of the original, inside the map or as
    references."""
    dm = DraftMap()
    for k, group in (("glycolysis", "pw1"), ("serine-metabolism", "pw3"), ("lipid-metabolism", "pw2")):
        dm.subsystems[k] = Subsystem(k, k, group, "#fff", 0, 0, 100, 100)
    for i, s in enumerate(dm.subsystems):
        dm.nodes[f"e{i}"] = Node(f"e{i}", "enzyme", s, 0, 0, 60, 30, f"e{i}", reactions=[f"R{i}"])
        dm.nodes[f"m{i}"] = Node(f"m{i}", "main", s, 50, 0, 60, 30, f"m{i}", ligand_id=1)
        dm.edges.append(Edge(f"e{i}", f"m{i}", "main", False, True))
    dm.links = [MetLink("m0", "m1"), MetLink("m1", "m2")]
    dm.orig_centres = {k: (0, 0) for k in dm.subsystems}
    maps_of = {"glycolysis": ["A", "B"], "serine-metabolism": ["B"], "lipid-metabolism": ["C"]}
    assignment, copies = split.duplicate_subsystems(dm, maps_of)
    assert copies == {("glycolysis", "B"): "glycolysis@b"}
    assert assignment == {"glycolysis": "A", "glycolysis@b": "B", "serine-metabolism": "B", "lipid-metabolism": "C"}
    assert dm.nodes["e0@b"].reactions == ["R0"] and dm.nodes["m0@b"].subsystem == "glycolysis@b"
    assert ("e0@b", "m0@b") in {(e.enzyme, e.met) for e in dm.edges}
    assert dm.orig_centres["glycolysis@b"] == (0, 0)
    maps = dict(split.split_map(dm, ["A", "B", "C"], assignment))
    # on B, the copy is linked to serine; the original on A refers to B, serine to A (original) and C
    assert {(l.a, l.b) for l in maps["B"].links} == {("m0@b", "m1")}
    assert maps["A"].nodes["m0"].stubs == ["B"]
    assert maps["B"].nodes["m1"].stubs == ["A", "C"]
    assert not getattr(maps["B"].nodes["m0@b"], "stubs", None)


def test_user_arrangement():
    """Given positions stay if the boxes are at least the gap apart; else boxes are moved apart as little as
    possible, in order."""
    w, h = np.array([100.0, 100.0, 50.0]), np.array([100.0, 60.0, 50.0])
    x, y = np.array([10.0, 200.0, 10.0]), np.array([20.0, 20.0, 300.0])
    ux, uy = layout.user_arrangement(w, h, x, y, gap=80)
    assert np.allclose(ux, x) and np.allclose(uy, y)
    x[1] = 60.0  # overlaps the first box
    assert len(layout.too_close(x, y, w, h, 80)) == 1
    ux, uy = layout.user_arrangement(w, h, x, y, gap=80)
    assert no_overlaps(ux, uy, w, h, 79.99)
    assert np.allclose(ux, [10, 190, 10], atol=1e-3) and np.allclose(uy, y, atol=1e-3)


def test_anneal_progress_and_stop():
    w = np.array([100.0, 50.0, 80.0])
    h = np.array([60.0, 50.0, 40.0])
    none = np.zeros(0, int)
    ports = Ports(["a", "b", "c"], w, h, none, none, np.zeros((0, 2)), np.zeros((0, 2)), [], h)
    seen = []
    x, y, *_ = layout.anneal_arrangement(ports, np.zeros(3), np.array([0.0, 100.0, 200.0]), iterations=1000,
                                         progress=lambda f, cost: seen.append(f), stop=lambda: len(seen) >= 3)
    assert seen == [0.0, 0.01, 0.02]
    assert no_overlaps(x, y, w, h, 80.0)


def test_draft_map_split_copies_and_arrangement(sources, model, tmp_path):  # noqa: F811
    """A subsystem marked for two maps is drawn on both; a map with the positions of all its subsystems in the
    sheet Arrangement is drawn so, boxes that are too close are moved apart, maps lacking positions are packed."""
    pytest.importorskip("openpyxl")
    from openpyxl import load_workbook

    from ecmapper.draft_map.assignment import read_assignment

    _, prep = prepared(sources, model, tmp_path)
    reference = tmp_path / "reference"
    run_draft_map(prep, reference, "--split", "--pack", "--min-map-reactions", "0")
    xlsx = tmp_path / "assignment.xlsx"
    run_draft_map(prep, tmp_path / "out", "--split", "--pack", "--min-map-reactions", "0", "--write-assignment",
                  str(xlsx))
    index = json.loads((reference / "maps.json").read_text())
    central, last = index["maps"][0], index["maps"][-1]
    copied = central["subsystems"][0]
    with open(reference / f"{split.slug(central['name'])}.pkl", "rb") as f:
        reference_map = pickle.load(f)
    reactions = {r for n in reference_map.nodes.values() if n.subsystem == copied for r in n.reactions}
    # the copy on the last map; the central map in a column, 100 apart; the last map: the boxes at one place,
    # (with the copy) or without the copy
    wb = load_workbook(xlsx)
    ws = wb["Subsystems"]
    column = [c.value for c in ws[1]].index(last["name"]) + 1
    for row in ws.iter_rows(min_row=2):
        if row[0].value == copied:
            ws.cell(row[0].row, column, "x")
    ws = wb.create_sheet("Arrangement")
    ws.append(["Map", "Subsystem ID", "x", "y"])
    given, top = {}, 0.0
    for k in central["subsystems"]:
        given[k] = (0.0, top)
        ws.append([central["name"], k, 0.0, top])
        top += reference_map.subsystems[k].h + 100
    for k in [*last["subsystems"], copied]:
        ws.append([last["name"], k, 0.0, 0.0])
    wb.save(xlsx)
    table = read_assignment(xlsx)
    assert table.subsystems[copied] == [central["name"], last["name"]]
    assert table.positions[central["name"]] == given

    edited = tmp_path / "edited"
    run_draft_map(prep, edited, "--split", "--pack", "--assignment", str(xlsx))
    maps = {m["name"]: m for m in json.loads((edited / "maps.json").read_text())["maps"]}
    assert [name for name, m in maps.items() if copied in m["subsystems"]] == [central["name"], last["name"]]
    assert reactions <= set(maps[central["name"]]["reactions"]) & set(maps[last["name"]]["reactions"])
    stats = json.loads((edited / "stats.json").read_text())
    assert stats["split"]["subsystem_copies"] == 1
    assert stats["maps"][central["name"]]["arranged_by"] == stats["maps"][last["name"]]["arranged_by"] == "table"
    with open(edited / f"{split.slug(central['name'])}.pkl", "rb") as f:
        m = pickle.load(f)
    for k, (x, y) in given.items():
        assert (m.subsystems[k].x, m.subsystems[k].y) == pytest.approx((x, y), abs=0.01)
    with open(edited / f"{split.slug(last['name'])}.pkl", "rb") as f:
        m = pickle.load(f)
    assert f"{copied}@{split.slug(last['name'])}" in m.subsystems
    boxes = list(m.subsystems.values())
    assert no_overlaps([b.x for b in boxes], [b.y for b in boxes], [b.w for b in boxes], [b.h for b in boxes], 79.99)

    # without the position of the copy, the last map is packed
    wb = load_workbook(xlsx)
    ws = wb["Arrangement"]
    ws.delete_rows(ws.max_row)
    wb.save(xlsx)
    run_draft_map(prep, edited, "--split", "--pack", "--assignment", str(xlsx))
    stats = json.loads((edited / "stats.json").read_text())
    assert stats["maps"][central["name"]]["arranged_by"] == "table"
    assert stats["maps"][last["name"]]["arranged_by"] == "packed"


def test_command_line():
    """The arguments of a namespace give it again (the GUI builds the maps with them)."""
    parser = run.build_parser()
    argv = ["p.pkl", "m.json", "m.svg", "out", "--split", "--pack", "--gap", "60", "--relayout", "spring",
            "--assignment", "a.xlsx"]
    args = parser.parse_args(argv)
    again = run.command_line(args, parser)
    assert parser.parse_args(again) == args
    assert again[:4] == argv[:4] and "--iterations" not in again


def test_borrowed_subsystems():
    """Rule B: subsystems whose reactions all have a home in subsystems with reactions of their own."""
    f = {k: frozenset(v) for k, v in {
        "glycolysis": ["FBA", "PGK", "TPI", "ENO"],  # ENO: its own
        "ppp": ["RPE", "TKT", "GND"],  # GND: its own
        "gluconeogenesis": ["PPC", "ME2", "PPS"],  # PPS: its own
        "photosynthesis": ["FBA", "PGK", "TPI", "RPE", "TKT", "PPC", "ME2"],  # all borrowed from the three above
        "c4": ["PPC", "ME2", "PYK"],  # PYK: its own
        "pdh1": ["PDH"],  # twins without own reactions: neither is borrowed
        "pdh2": ["PDH"],
        "plant": ["FBA", "RBPC"],  # RBPC: its own, stays although FBA is glycolysis'
    }.items()}
    assert core.borrowed_subsystems(f) == [("photosynthesis", ["c4", "gluconeogenesis", "glycolysis", "plant", "ppp"])]
    f["glycolysis"] = f["glycolysis"] | {"PYK"}
    # now c4 is the smallest borrowed subsystem; photosynthesis follows
    assert core.borrowed_subsystems(f) == [("c4", ["gluconeogenesis", "glycolysis"]),
                                           ("photosynthesis", ["gluconeogenesis", "glycolysis", "plant", "ppp"])]
    # a chain: x borrows from y only once z (a subset of y, borrowed from y and w) is gone
    g = {k: frozenset(v) for k, v in {"w": ["R4", "R5"], "x": ["R1", "R2"], "y": ["R1", "R2", "R3"],
                                      "z": ["R3", "R4"]}.items()}
    assert core.borrowed_subsystems(g) == []  # y has no own reaction while x and z are there
    del g["z"]
    assert core.borrowed_subsystems(g) == [("x", ["y"])]
    assert core.borrowed_subsystems({"a": frozenset(), "b": frozenset(["R1"])}) == []


def test_drop_redundant_subsystems():
    """A: proper subsets and B: borrowed subsystems are removed; C, D: (mostly) covered subsystems and E: few
    matched nodes are set aside."""
    dm = DraftMap()
    content = {
        "big": ["R1", "R2", "R3", "R4"],
        "subset": ["R1", "R2"],  # A: proper subset of big
        "equal1": ["R30"],  # equal sets: the first (by name) is set aside (C), the other stays
        "equal2": ["R30"],
        "mixed": ["R1", "R50", *(f"R{i}" for i in range(5, 14))],  # 10 of 11 in big and other (D)
        "other": [*(f"R{i}" for i in range(5, 14)), "R40", "R41"],
        "spread": ["R3", "R41", "R60"],  # 2 of 3 elsewhere: stays with coverage 0.9
        "borrowed": ["R2", "R41"],  # B: in big, other and spread, which have reactions of their own (R4, R40, R60)
        "sparse": ["R20"],  # one of five enzyme nodes matched (E)
    }
    for k, rs in content.items():
        dm.subsystems[k] = Subsystem(k, k, "pw1", "#fff", 0, 0, 100, 100)
        for r in rs:
            dm.nodes[f"{k}_{r}"] = Node(f"{k}_{r}", "enzyme", k, 0, 0, 60, 30, r, reactions=[r])
    for i in range(4):
        dm.nodes[f"sparse_{i}"] = Node(f"sparse_{i}", "enzyme", "sparse", 0, 0, 60, 30, "")
    result = core.drop_redundant_subsystems(dm, coverage=0.9, min_matched=0.25)
    assert [(d["subsystem"], d["rule"]) for d in result["dropped"]] == [("subset", "A"), ("borrowed", "B")]
    assert result["dropped"][0]["subset_of"] == "big"
    assert result["dropped"][1]["drawn_in"] == ["big", "other", "spread"]
    assert not {"subset", "borrowed"} & set(dm.subsystems) and not dm.members("subset") and not dm.members("borrowed")
    assert {a["subsystem"]: a["rule"] for a in result["set_aside"]} == {"equal1": "C", "mixed": "D", "sparse": "E"}
    assert {k for k, s in dm.subsystems.items() if not s.aside} == {"big", "equal2", "spread", "other"}
    # no rules
    dm2 = DraftMap()
    dm2.subsystems = {k: Subsystem(k, k, "pw1", "#fff", 0, 0, 100, 100) for k in ("a", "b", "c")}
    dm2.nodes = {k: Node(k, "enzyme", k, 0, 0, 60, 30, "", reactions=rs)
                 for k, rs in (("a", ["R1"]), ("b", ["R1", "R2"]), ("c", ["R2"]))}
    assert core.drop_redundant_subsystems(dm2, drop_subsets=False, coverage=2.0,
                                          drop_borrowed=False) == dict(dropped=[], set_aside=[])
    assert set(dm2.subsystems) == {"a", "b", "c"}


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
    # every reaction of a set-aside subsystem is also on the draft map (rule C)
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
        if prog == "neato":
            # pinned softly: the map nodes drift a little (back where Graphviz moved the whole drawing)
            for i in pinned:
                assert math.dist((s.nodes[i].x, s.nodes[i].y), (dm.nodes[i].x, dm.nodes[i].y)) < 150
    # the ghost nodes of soft pinning: pinned at the map nodes, one heavy edge each, not in the result
    s = relayout.subsystem_map(dm, "S")
    spec, box, owner = relayout.graphviz_spec(s, **relayout.GRAPHVIZ["neato"])
    ghosts = {name: attrs for name, attrs in spec["nodes"] if name.startswith(relayout.GHOST)}
    assert set(ghosts) == {relayout.GHOST + i for i in pinned} and all(a["pos"].endswith("!") for a in ghosts.values())
    assert sum(len(e) == 3 for e in spec["edges"]) == len(pinned)
    positions = relayout.graphviz_positions(spec)
    relayout.apply_graphviz(s, positions, box, owner)
    assert set(s.nodes) == set(relayout.subsystem_map(dm, "S").nodes)
    # hard pins and dot: no ghosts
    for args in (dict(prog="fdp", K="2.0"), dict(prog="dot")):
        assert not any(n.startswith(relayout.GHOST) for n, _ in relayout.graphviz_spec(s, **args)[0]["nodes"])


def test_relayout_subsystems():
    pytest.importorskip("networkx")
    pytest.importorskip("pygraphviz")
    dm = relayout_map()
    before = {i: (n.x, n.y) for i, n in dm.nodes.items()}
    rows = relayout.relayout_subsystems(dm)
    assert [r["subsystem"] for r in rows] == ["S"]  # T has no added reactions
    assert {i: (n.x, n.y) for i, n in dm.nodes.items() if n.subsystem == "T"} == {"t1": before["t1"], "t2": before["t2"]}
    row = rows[0]
    assert set(row) >= {"current", "place", "spring", "fdp", "dot", "neato", "layout"}
    best_keeping = min(r["crossings_main"] for m, r in row.items() if m != "layout" and relayout._base(m) in relayout.KEEPING)
    chosen = row[row["layout"]]["crossings_main"]
    if relayout._base(row["layout"]) in relayout.REDRAWING:
        assert chosen < 0.75 * best_keeping and best_keeping - chosen >= 3
    else:
        assert chosen == best_keeping
    assert chosen <= row["current"]["crossings_main"]
    assert "outside" not in row  # no cycle
    # the minimum gain switches the redrawing layouts off
    rows = relayout.relayout_subsystems(relayout_map(), min_gain=10**6)
    assert relayout._base(rows[0]["layout"]) in relayout.KEEPING
    # no layouts: nothing changes
    dm = relayout_map()
    rows = relayout.relayout_subsystems(dm, methods=[], untangle=False)
    assert rows[0]["layout"] == "current" and {i: (n.x, n.y) for i, n in dm.nodes.items()} == before


def ring_map():
    """A 5-metabolite cycle of the BRENDA map with an added reaction inside it (2 members) and a chord (3)."""
    dm = DraftMap()
    dm.subsystems["R"] = Subsystem("R", "R", "pw1", "#fff", 0, 0, 1000, 1000)
    r = 300.0
    for i in range(5):
        a = 2 * math.pi * i / 5
        dm.nodes[f"m{i}"] = Node(f"m{i}", "main", "R", r * math.cos(a), r * math.sin(a), 60, 30, f"m{i}")
        b = 2 * math.pi * (i + 0.5) / 5
        dm.nodes[f"e{i}"] = Node(f"e{i}", "enzyme", "R", r * math.cos(b), r * math.sin(b), 60, 30, f"e{i}",
                                 reactions=[f"E{i}"], status="on-map")
        dm.edges += [Edge(f"e{i}", f"m{i}", "main", True, False), Edge(f"e{i}", f"m{(i + 1) % 5}", "main", False, True)]
    dm.nodes["a1"] = Node("a1", "enzyme", "R", 0, 0, 60, 30, "A1", reactions=["A1"], status="added")
    dm.nodes["n1"] = Node("n1", "main", "R", 0, 60, 60, 30, "n1", status="new")
    dm.nodes["c1"] = Node("c1", "cofactor", "R", 40, 40, 40, 20, "ATP", status="added")
    dm.nodes["a2"] = Node("a2", "enzyme", "R", 20, -40, 60, 30, "A2", reactions=["A2"], status="added")
    dm.edges += [Edge("a1", "m0", "main", True, False), Edge("a1", "m1", "main", False, True),
                 Edge("a1", "n1", "main", False, True), Edge("a1", "c1", "cofactor", True, False),
                 Edge("a2", "m0", "main", True, False), Edge("a2", "m2", "main", False, True),
                 Edge("a2", "m3", "main", False, True)]
    return dm


def _inside(poly, x, y):
    inside = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def test_outside_cycle():
    pytest.importorskip("networkx")
    dm = ring_map()
    s = relayout.subsystem_map(dm, "R")
    assert relayout.outside_cycle(s)
    poly = [(s.nodes[f"m{i}"].x, s.nodes[f"m{i}"].y) for i in range(5)]
    for i in range(5):  # the ring is unmoved
        for j in (f"m{i}", f"e{i}"):
            assert (s.nodes[j].x, s.nodes[j].y) == (dm.nodes[j].x, dm.nodes[j].y)
    # the added reaction between m0 and m1 is outside, beyond the middle of its two members, with its own nodes
    assert not _inside(poly, s.nodes["a1"].x, s.nodes["a1"].y)
    mid = ((poly[0][0] + poly[1][0]) / 2, (poly[0][1] + poly[1][1]) / 2)
    assert math.hypot(s.nodes["a1"].x, s.nodes["a1"].y) > math.hypot(*mid)
    for j in ("n1", "c1"):
        assert not _inside(poly, s.nodes[j].x, s.nodes[j].y)
    # the chord (3 members) stays
    assert (s.nodes["a2"].x, s.nodes["a2"].y) == (20, -40)
    # without a cycle nothing happens
    s = relayout.subsystem_map(relayout_map(), "S")
    assert not relayout.outside_cycle(s)


def test_relayout_subsystems_outside():
    pytest.importorskip("networkx")
    dm = ring_map()
    rows = relayout.relayout_subsystems(dm, methods=["spring"], untangle=False)
    row = rows[0]
    assert "outside" in row and relayout._base(row["layout"]) in relayout.KEEPING
    assert row["outside"]["crossings_main"] <= row["current"]["crossings_main"]
    rows = relayout.relayout_subsystems(ring_map(), methods=["spring"], untangle=False, outside=False)
    assert "outside" not in rows[0]


MAP_NODES = ("m1", "m2", "m3", "m4", "e1", "e2")


def crossed_map():
    """relayout_map with the added reaction across both map reactions (a1 above them, its new metabolite below)."""
    dm = relayout_map()
    dm.nodes["a1"].x, dm.nodes["a1"].y, dm.nodes["c1"].x, dm.nodes["c1"].y = 300, -100, 340, -60
    dm.nodes["n1"].x, dm.nodes["n1"].y = 300, 400
    return dm


def main_crossings(s):
    from ecmapper.pwl.layout import crossing_nodes

    return crossing_nodes(s, {i: (n.x, n.y) for i, n in s.nodes.items() if n.kind != "cofactor"})[0]


def assert_no_overlaps(s):
    boxes = [(i, n.box) for i, n in sorted(s.nodes.items())]
    for k, (i, a) in enumerate(boxes):
        for j, b in boxes[k + 1:]:
            assert not (a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]), (i, j)


def test_place_added_keeps_map_nodes():
    from ecmapper.pwl.layout import STEP

    dm = crossed_map()
    s = relayout.subsystem_map(dm, "S")
    assert main_crossings(s) > 0
    assert relayout.place_added(s)
    for i in MAP_NODES:
        assert (s.nodes[i].x, s.nodes[i].y) == (dm.nodes[i].x, dm.nodes[i].y)
    assert main_crossings(s) == 0
    assert_no_overlaps(s)
    assert all(math.isfinite(n.x) and math.isfinite(n.y) for n in s.nodes.values())
    # next to its metabolites, the new metabolite next to it, the cofactor moves with its enzyme node
    a1 = (s.nodes["a1"].x, s.nodes["a1"].y)
    for j in ("m2", "m3", "n1"):
        assert math.dist(a1, (s.nodes[j].x, s.nodes[j].y)) < STEP
    assert s.nodes["c1"].x - a1[0] == pytest.approx(40) and s.nodes["c1"].y - a1[1] == pytest.approx(40)
    # deterministic
    again = relayout.subsystem_map(dm, "S")
    relayout.place_added(again)
    assert {i: (n.x, n.y) for i, n in again.nodes.items()} == {i: (n.x, n.y) for i, n in s.nodes.items()}


def test_place_added_outside_cycle():
    pytest.importorskip("networkx")
    dm = ring_map()
    s = relayout.subsystem_map(dm, "R")
    assert relayout.place_added(s)
    poly = [(s.nodes[f"m{i}"].x, s.nodes[f"m{i}"].y) for i in range(5)]
    for i in range(5):  # the ring is unmoved
        for j in (f"m{i}", f"e{i}"):
            assert (s.nodes[j].x, s.nodes[j].y) == (dm.nodes[j].x, dm.nodes[j].y)
    # the added reaction between neighbours m0 and m1 is outside with its own nodes; the chord a2 may stay inside
    for j in ("a1", "n1", "c1"):
        assert not _inside(poly, s.nodes[j].x, s.nodes[j].y)
    assert main_crossings(s) == 0
    assert_no_overlaps(s)
    s = relayout.subsystem_map(dm, "R")
    assert relayout.place_added(s, outside=False)  # without the rule of the cycle: the ring is still unmoved
    for i in range(5):
        for j in (f"m{i}", f"e{i}"):
            assert (s.nodes[j].x, s.nodes[j].y) == (dm.nodes[j].x, dm.nodes[j].y)


def test_place_added_chain_and_linked():
    """A chain of added reactions from a map metabolite, and a linked metabolite of a map enzyme, all far away."""
    from ecmapper.pwl.layout import STEP

    dm = DraftMap()
    dm.subsystems["C"] = Subsystem("C", "C", "pw1", "#fff", 0, 0, 1000, 600)
    dm.nodes["m1"] = Node("m1", "main", "C", 0, 0, 60, 30, "m1")
    dm.nodes["m2"] = Node("m2", "main", "C", 400, 0, 60, 30, "m2")
    dm.nodes["e1"] = Node("e1", "enzyme", "C", 200, 0, 60, 30, "e1", reactions=["E1"], status="on-map")
    for i in ("a1", "a2"):
        dm.nodes[i] = Node(i, "enzyme", "C", 1000, 1000, 60, 30, i, reactions=[i.upper()], status="added")
    for i in ("x", "y"):
        dm.nodes[i] = Node(i, "main", "C", 1000, 1000, 60, 30, i, status="new")
    dm.nodes["l1"] = Node("l1", "main", "C", -2000, 0, 60, 30, "l1", status="linked")
    dm.edges = [Edge("e1", "m1", "main", True, False), Edge("e1", "m2", "main", False, True),
                Edge("e1", "l1", "main", False, True),
                Edge("a1", "m2", "main", True, False), Edge("a1", "x", "main", False, True),
                Edge("a2", "x", "main", True, False), Edge("a2", "y", "main", False, True)]
    s = relayout.subsystem_map(dm, "C")
    assert relayout.place_added(s)
    for i in ("m1", "m2", "e1"):
        assert (s.nodes[i].x, s.nodes[i].y) == (dm.nodes[i].x, dm.nodes[i].y)

    def near(i, j):
        return math.dist((s.nodes[i].x, s.nodes[i].y), (s.nodes[j].x, s.nodes[j].y)) < STEP

    assert near("a1", "m2") and near("x", "a1") and near("a2", "x") and near("y", "a2") and near("l1", "e1")
    assert main_crossings(s) == 0
    assert_no_overlaps(s)


def test_place_added_only_added_nodes():
    from ecmapper.pwl.layout import STEP

    dm = DraftMap()
    dm.subsystems["A"] = Subsystem("A", "A", "pw1", "#fff", 0, 0, 1000, 600)
    dm.nodes["a1"] = Node("a1", "enzyme", "A", 100, 100, 60, 30, "A1", reactions=["A1"], status="added")
    dm.nodes["c1"] = Node("c1", "cofactor", "A", 140, 140, 40, 20, "ATP", status="added")
    for i in ("n1", "n2"):
        dm.nodes[i] = Node(i, "main", "A", 100, 100, 60, 30, i, status="new")
    dm.edges = [Edge("a1", "n1", "main", True, False), Edge("a1", "n2", "main", False, True),
                Edge("a1", "c1", "cofactor", True, False)]
    s = relayout.subsystem_map(dm, "A")
    assert relayout.place_added(s)
    assert (s.nodes["a1"].x, s.nodes["a1"].y) == (100, 100)  # no placed partner: it keeps its place
    assert_no_overlaps(s)
    for i in ("n1", "n2"):
        assert math.dist((s.nodes[i].x, s.nodes[i].y), (100, 100)) < STEP
    # nothing to place
    assert not relayout.place_added(relayout.subsystem_map(relayout_map(), "T"))


def test_relayout_subsystems_place():
    assert "place" in relayout.METHODS and "place" in relayout.KEEPING
    assert relayout.available_methods(["place"]) == ["place"]
    dm = crossed_map()
    before = {i: (n.x, n.y) for i, n in dm.nodes.items()}
    rows = relayout.relayout_subsystems(dm, methods=["place"], untangle=False)
    row = rows[0]
    assert "place" in row and "outside" not in row and row["layout"] in relayout.KEEPING
    assert row["place"]["crossings_main"] < row["current"]["crossings_main"]
    assert row["layout"] == "place"
    for i in MAP_NODES:
        assert (dm.nodes[i].x, dm.nodes[i].y) == before[i]
    assert "placed again next to their metabolites in 1" in dm.log[-1]
    # without place, as before
    rows = relayout.relayout_subsystems(crossed_map(), methods=[], untangle=False)
    assert "place" not in rows[0] and rows[0]["layout"] == "current"
    pytest.importorskip("networkx")
    rows = relayout.relayout_subsystems(ring_map(), methods=["place", "spring"], untangle=False)
    assert "place" in rows[0] and "outside" not in rows[0] and "spring" in rows[0]
    rows = relayout.relayout_subsystems(ring_map(), methods=["place", "spring"], untangle=False, outside=False)
    assert "place" in rows[0]
    rows = relayout.relayout_subsystems(ring_map(), methods=["spring"], untangle=False)
    assert "outside" in rows[0] and "place" not in rows[0]


def walled_map():
    """ring_map with a map metabolite z inside the ring and w outside, and an added reaction z -> w (a3, with a
    cofactor and a new metabolite): it cannot be drawn without crossing the ring."""
    dm = ring_map()
    dm.nodes["z"] = Node("z", "main", "R", 0, 0, 60, 30, "z")
    dm.nodes["w"] = Node("w", "main", "R", 800, 0, 60, 30, "w")
    dm.nodes["a3"] = Node("a3", "enzyme", "R", 400, 0, 60, 30, "A3", reactions=["A3"], status="added")
    dm.nodes["c3"] = Node("c3", "cofactor", "R", 440, 40, 40, 20, "ATP", status="added")
    dm.nodes["n3"] = Node("n3", "main", "R", 400, 80, 60, 30, "n3", status="new")
    dm.edges += [Edge("a3", "z", "main", True, False), Edge("a3", "w", "main", False, True),
                 Edge("a3", "n3", "main", False, True), Edge("a3", "c3", "cofactor", True, False)]
    dm.links.append(MetLink("n3", "w", "added"))
    return dm


def all_crossings(s):
    return layout.internal_metrics(s)["crossings_all"]


def test_place_conservative_adds_no_crossings():
    dm = crossed_map()
    s = relayout.subsystem_map(dm, "S")
    assert main_crossings(s) > 0
    assert relayout.place_conservative(s) == set()  # the added reaction fits between the two map reactions
    for i in MAP_NODES:
        assert (s.nodes[i].x, s.nodes[i].y) == (dm.nodes[i].x, dm.nodes[i].y)
    assert all_crossings(s) == 0
    assert_no_overlaps(s)
    # deterministic
    again = relayout.subsystem_map(dm, "S")
    relayout.place_conservative(again)
    assert {i: (n.x, n.y) for i, n in again.nodes.items()} == {i: (n.x, n.y) for i, n in s.nodes.items()}


def test_place_conservative_leaves_out():
    pytest.importorskip("networkx")
    dm = walled_map()
    s = relayout.subsystem_map(dm, "R")
    s.links = list(dm.links)
    before = all_crossings(s)
    removed = relayout.place_conservative(s)
    # a3 crosses the ring wherever it goes: left out with its cofactor and new metabolite; a1 and a2 fit
    assert removed == {"a3", "c3", "n3"}
    assert not removed & set(s.nodes) and not any(e.enzyme == "a3" for e in s.edges) and not s.links
    assert {"a1", "a2", "n1", "c1", "z", "w"} <= set(s.nodes)
    assert all_crossings(s) <= before and all_crossings(s) == 0
    assert_no_overlaps(s)
    for i in range(5):  # the ring is unmoved
        for j in (f"m{i}", f"e{i}"):
            assert (s.nodes[j].x, s.nodes[j].y) == (dm.nodes[j].x, dm.nodes[j].y)
    # put into the subsystem by hand: drawn anyway
    dm = walled_map()
    dm.nodes["a3"].assigned = True
    s = relayout.subsystem_map(dm, "R")
    assert relayout.place_conservative(s) == set()
    assert (s.nodes["a3"].x, s.nodes["a3"].y) == (400, 0)


def test_relayout_subsystems_conservative():
    pytest.importorskip("networkx")
    assert relayout.CONSERVATIVE not in relayout.METHODS
    dm = walled_map()
    rows = relayout.relayout_subsystems(dm, methods=[relayout.CONSERVATIVE])
    row = rows[0]
    # no other layout, not even untangled; the left-out reaction is removed from the draft map
    assert set(row) == {"subsystem", "layout", "current", relayout.CONSERVATIVE, "left_out"}
    assert row["layout"] == relayout.CONSERVATIVE and row["left_out"] == ["A3"]
    assert not {"a3", "c3", "n3"} & set(dm.nodes) and not dm.links
    assert row[relayout.CONSERVATIVE]["crossings_all"] == 0
    assert "1 model reactions left out" in dm.log[-1]
    with pytest.raises(ValueError):
        relayout.relayout_subsystems(walled_map(), methods=[relayout.CONSERVATIVE, "spring"])


def test_relayout_by_hand():
    pytest.importorskip("networkx")
    s = relayout.subsystem_map(walled_map(), "R")
    results, _ = relayout.compare_layouts(s, untangle=False)
    for name in ("place", "spring", "circle", relayout.CONSERVATIVE):
        dm = walled_map()
        rows = relayout.relayout_subsystems(dm, [relayout.CONSERVATIVE], untangle=False, by_hand={"R": name})
        assert [r["layout"] for r in rows] == [name] and rows[0]["by_hand"]
        assert ("a3" in dm.nodes) == (name != relayout.CONSERVATIVE)  # left out only by the conservative placement
        # the same drawing as in the comparison
        c = relayout.subsystem_map(dm, "R")
        layout.condense_subsystems(c, keep_crossings=name == relayout.CONSERVATIVE)
        shown = results[name]["map"]
        assert {i: (round(n.x, 3), round(n.y, 3)) for i, n in c.nodes.items()} == \
            {i: (round(n.x, 3), round(n.y, 3)) for i, n in shown.nodes.items()}
    assert relayout.layout_methods("spring+untangle") == (["place", "spring"], False, True)
    # the other subsystems as step 3b chooses
    dm = walled_map()
    dm.subsystems["T"] = Subsystem("T", "T", "pw1", "#fff", 0, 0, 1000, 600)
    rows = relayout.relayout_subsystems(dm, ["place"], untangle=False, by_hand={"R": relayout.CONSERVATIVE, "X": "place"})
    assert [(r["subsystem"], r["layout"]) for r in rows] == [("R", relayout.CONSERVATIVE)]
    # --layout-of
    assert run.parse_layout_of(" a=place, b=spring+untangle ") == {"a": "place", "b": "spring+untangle"}
    assert run.format_layout_of({"b": "spring", "a": "place"}) == "a=place,b=spring"
    parser = run.build_parser()
    for bad in ("a", "a=nonsense"):
        args = parser.parse_args(["p", "j", "s", "o", "--layout-of", bad])
        with pytest.raises(SystemExit):
            run.check_args(parser, args)
    run.check_args(parser, parser.parse_args(["p", "j", "s", "o", "--layout-of", "a=neato+untangle,b=circle"]))


def template_map():
    """The BRENDA map with t1 (m1 + ATP -> m2 + ADP) and a2's template t2 (m2 -> m3), and the draft map in which
    t1 and t2 were pruned: a1 (m1 + GTP -> m2 + GDP) at t1's place, a2 with its new metabolite n3 at m3's place,
    a3 without template node."""
    from types import SimpleNamespace

    template = DraftMap()
    template.subsystems["S"] = Subsystem("S", "S", "pw1", "#fff", 0, 0, 1000, 600)
    for i, x, y, kind in (("m1", 0, 0, "main"), ("m2", 400, 0, "main"), ("m3", 800, 0, "main"),
                          ("atp", 160, 60, "cofactor"), ("adp", 240, 60, "cofactor"), ("t1", 200, 0, "enzyme"),
                          ("t2", 600, 0, "enzyme")):
        template.nodes[i] = Node(i, kind, "S", x, y, 60, 30, i)
    template.edges = [Edge("t1", "m1", "main", True, False), Edge("t1", "m2", "main", False, True),
                      Edge("t1", "atp", "cofactor", True, False), Edge("t1", "adp", "cofactor", False, True),
                      Edge("t2", "m2", "main", True, False), Edge("t2", "m3", "main", False, True)]
    keys = {"m1": {"K1"}, "m2": {"K2"}, "m3": {"K3"}, "atp": {"KATP"}, "adp": {"KADP"}}
    md = SimpleNamespace(
        metabolites={"m1_c": dict(keys={"K1"}), "gtp_c": dict(keys={"KGTP"}), "gdp_c": dict(keys={"KGDP"}),
                     "m3_c": dict(keys={"K3"}), "pi_c": dict(keys={"KPI"})},
        node_met_keys={(e.enzyme, e.met): keys[e.met] for e in template.edges},
        matches={"A1": dict(candidates=[dict(node_id="t1", via=("ec",))]),
                 "A2": dict(candidates=[dict(node_id="t2", via=("reaction-id",))])})
    completion = core.Completion.__new__(core.Completion)
    completion.md, completion.template, completion.template_edges = md, template, template.edges_of()
    dm = DraftMap()
    dm.subsystems["S"] = Subsystem("S", "S", "pw1", "#fff", 0, 0, 1000, 600)
    for i in ("m1", "m2"):
        dm.nodes[i] = copy.deepcopy(template.nodes[i])
    dm.nodes["e1"] = Node("e1", "enzyme", "S", 200, 200, 60, 30, "e1", reactions=["E1"], status="on-map")
    dm.nodes["a1"] = Node("a1", "enzyme", "S", 230, 90, 60, 30, "A1", reactions=["A1"], status="added")
    dm.nodes["gtp"] = Node("gtp", "cofactor", "S", 200, 135, 40, 20, "GTP", status="added")
    dm.nodes["gdp"] = Node("gdp", "cofactor", "S", 300, 135, 40, 20, "GDP", status="added")
    dm.nodes["a2"] = Node("a2", "enzyme", "S", 500, 300, 60, 30, "A2", reactions=["A2"], status="added")
    dm.nodes["n3"] = Node("n3", "main", "S", 700, 300, 60, 30, "n3", status="new")
    dm.nodes["a3"] = Node("a3", "enzyme", "S", 0, 300, 60, 30, "A3", reactions=["A3"], status="added")
    dm.nodes["n4"] = Node("n4", "main", "S", 0, 400, 60, 30, "n4", status="new")
    dm.edges = [Edge("e1", "m1", "main", True, False), Edge("e1", "m2", "main", False, True),
                Edge("a1", "m1", "main", True, False), Edge("a1", "m2", "main", False, True),
                Edge("a1", "gtp", "cofactor", True, False), Edge("a1", "gdp", "cofactor", False, True),
                Edge("a2", "m2", "main", True, False), Edge("a2", "n3", "main", False, True),
                Edge("a3", "m1", "main", True, False), Edge("a3", "n4", "main", False, True)]
    completion.dm = dm
    return completion, dm


def test_template_places():
    completion, dm = template_map()
    # step 3: the place of the pruned node with the reaction's EC number, those of its metabolites and cofactors
    assert completion.template_spot("A1", "S") == "t1" and completion.template_spot("A1", "T") is None
    assert completion.template_spot("A3", "S") is None
    used = set()
    assert completion.template_met("t1", "m1_c", "main", -1, used) == ("m1", 0, 0)
    others = frozenset({"K1", "K2", "KGDP"})
    assert completion.template_met("t1", "gtp_c", "cofactor", -1, used, others) == ("atp", 160, 60)  # its side
    assert completion.template_met("t1", "gdp_c", "cofactor", 1, used, others) == ("adp", 240, 60)
    assert completion.template_met("t1", "pi_c", "cofactor", 1, used, others) is None  # no template cofactor left
    dm.nodes["a1"].template, dm.nodes["a2"].template = ("t1", 200, 0), ("t2", 600, 0)
    dm.nodes["gtp"].template, dm.nodes["gdp"].template = ("atp", 160, 60), ("adp", 240, 60)
    dm.nodes["n3"].template = ("m3", 800, 0)
    s = relayout.subsystem_map(dm, "S")
    removed = relayout.place_template(s)
    # a1 and a2 at the template's places, a3 (no template node) left out with its new metabolite
    assert removed == {"a3", "n4"}
    assert {i: (s.nodes[i].x, s.nodes[i].y) for i in ("a1", "gtp", "gdp", "a2", "n3")} == {
        "a1": (200, 0), "gtp": (160, 60), "gdp": (240, 60), "a2": (600, 0), "n3": (800, 0)}
    assert (s.nodes["e1"].x, s.nodes["e1"].y) == (200, 200)  # the map is unchanged
    # a template metabolite that is still drawn (no free place), or a new metabolite without one: left out
    dm.nodes["n3"].template = ("m2", 400, 0)
    assert "a2" in relayout.place_template(relayout.subsystem_map(dm, "S"))
    del dm.nodes["n3"].template
    assert "a2" in relayout.place_template(relayout.subsystem_map(dm, "S"))
    # a place from which a link would pass through a node (m1 on the way to m2): left out
    dm.nodes["a1"].template = ("t1", -200, 0)
    assert "a1" in relayout.place_template(relayout.subsystem_map(dm, "S"))
    dm.nodes["a1"].template = ("t1", 200, 0)
    # step 3b
    rows = relayout.relayout_subsystems(dm, [relayout.TEMPLATE])
    assert rows[0]["layout"] == relayout.TEMPLATE and rows[0]["left_out"] == ["A2", "A3"]
    assert set(dm.nodes) == {"m1", "m2", "e1", "a1", "gtp", "gdp"} and "without such a place" in dm.log[-1]
    with pytest.raises(ValueError):
        relayout.relayout_subsystems(dm, [relayout.TEMPLATE, "place"])
    assert relayout.layout_methods(relayout.TEMPLATE) == ([relayout.TEMPLATE], False, False)


def test_condense_keeps_crossings():
    dm = crossed_map()
    s = relayout.subsystem_map(dm, "S")
    relayout.place_conservative(s)
    for args in (dict(compactness=0.5), dict(compactness=5.0)):
        c = relayout.subsystem_map(s, "S")
        layout.condense_subsystems(c, keep_crossings=True, **args)
        assert all_crossings(c) == 0


def test_shrink_uniform():
    dm = crossed_map()
    s = relayout.subsystem_map(dm, "S")
    relayout.place_conservative(s)
    c = relayout.subsystem_map(s, "S")
    fx, fy = layout.shrink_uniform(c, "S")
    assert 0.1 <= fx < 1 and 0.1 <= fy < 1
    # enzyme and main metabolite nodes are scaled, the cofactor keeps its offset to its enzyme node
    for i, n in s.nodes.items():
        if n.kind != "cofactor":
            assert (c.nodes[i].x, c.nodes[i].y) == pytest.approx((fx * n.x, fy * n.y))
    assert (c.nodes["c1"].x - c.nodes["a1"].x, c.nodes["c1"].y - c.nodes["a1"].y) == pytest.approx(
        (s.nodes["c1"].x - s.nodes["a1"].x, s.nodes["c1"].y - s.nodes["a1"].y))
    assert all_crossings(c) == 0
    assert_no_overlaps(c)
    # two nodes that are closer than the gaps already do not come closer: no shrink along x
    tight = relayout.subsystem_map(s, "S")
    tight.nodes["m2"].x = tight.nodes["e1"].x + 70
    assert layout.shrink_uniform(tight, "S")[0] == 1.0
    # step 4: the box is fitted around the shrunk drawing, smaller than without condensing
    sizes = {}
    for uniform in (True, False):
        c = relayout.subsystem_map(s, "S")
        layout.condense_subsystems(c, keep_order=None, uniform=uniform) if uniform else layout.fit_boxes_original(c)
        sizes[uniform] = c.subsystems["S"].w * c.subsystems["S"].h
        if uniform:
            dx = c.nodes["m1"].x - fx * s.nodes["m1"].x
            assert c.nodes["m4"].x == pytest.approx(fx * s.nodes["m4"].x + dx)
    assert sizes[True] < sizes[False]
    rows = relayout.relayout_subsystems(dm, [relayout.CONSERVATIVE], uniform=True)
    assert rows[0]["layout"] == relayout.CONSERVATIVE


def test_compare_layouts():
    pytest.importorskip("networkx")
    s = relayout.subsystem_map(walled_map(), "R")
    results, chosen = relayout.compare_layouts(s, untangle=False)
    assert {"current", "place", relayout.CONSERVATIVE, relayout.TEMPLATE, "spring", "circle"} <= set(results)
    assert chosen in results and chosen != relayout.CONSERVATIVE and relayout._base(chosen) != "circle"
    assert results[relayout.CONSERVATIVE]["left_out"] == ["A3"] and results["current"]["left_out"] == []
    for r in results.values():
        assert list(r["map"].subsystems) == ["R"] and r["map"].nodes
    assert "a3" in results["place"]["map"].nodes and "a3" not in results[relayout.CONSERVATIVE]["map"].nodes
    assert sorted(results, key=relayout.layout_key)[:4] == ["current", relayout.CONSERVATIVE, relayout.TEMPLATE,
                                                            "place"]
    # the input is unchanged; without added reactions only the current drawing (and the layout chosen by hand only)
    assert "a3" in s.nodes
    results, chosen = relayout.compare_layouts(relayout.subsystem_map(relayout_map(), "T"))
    assert set(results) - {relayout.AUTOLAYOUT} == {"current"} and chosen is None


def test_autolayout(monkeypatch):
    pytest.importorskip("sbml_autolayout")
    dm = ring_map()
    drawings = []
    for _ in range(2):
        s = relayout.subsystem_map(dm, "R")
        assert relayout.LAYOUTS[relayout.AUTOLAYOUT](s) is None
        drawings.append({i: (n.x, n.y) for i, n in s.nodes.items()})
    assert drawings[0] == drawings[1]  # reproducible
    # all enzyme and main nodes drawn again, at places of their own; the cofactor moves with its enzyme node
    main = [i for i, n in s.nodes.items() if n.kind != "cofactor"]
    assert all(drawings[0][i] != (dm.nodes[i].x, dm.nodes[i].y) for i in main)
    assert len({drawings[0][i] for i in main}) == len(main)
    assert s.nodes["c1"].x - s.nodes["a1"].x == pytest.approx(40) and s.nodes["c1"].y - s.nodes["a1"].y == pytest.approx(40)
    # a known cycle (found by the names of the BRENDA map, the titles) is drawn on a circle, which condensing keeps
    tca = DraftMap()
    tca.subsystems["C"] = Subsystem("C", "C", "pw1", "#fff", 0, 0, 1000, 1000)
    names = ["oxaloacetate", "citrate", "isocitrate", "2-oxoglutarate", "succinyl-CoA", "succinate", "fumarate", "L-malate"]
    for i, name in enumerate(names):
        tca.nodes[f"m{i}"] = Node(f"m{i}", "main", "C", 100 * i, 0, 60, 30, f"m{i}", title=name)
        tca.nodes[f"e{i}"] = Node(f"e{i}", "enzyme", "C", 100 * i + 50, 60, 60, 30, f"e{i}", reactions=[f"E{i}"],
                                  status="on-map")
        tca.edges += [Edge(f"e{i}", f"m{i}", "main", True, False), Edge(f"e{i}", f"m{(i + 1) % 8}", "main", False, True)]
    tca.nodes["c1"] = Node("c1", "cofactor", "C", 60, 90, 40, 20, "NADH")
    tca.edges.append(Edge("e2", "c1", "cofactor", False, True))
    relayout.autolayout_layout(tca)
    assert [sorted(group) for group in tca.subsystems["C"].rigid] == [sorted([*tca.nodes])]
    layout.condense_subsystems(tca)
    centre = [sum(getattr(tca.nodes[f"m{i}"], a) for i in range(8)) / 8 for a in "xy"]
    radii = [math.dist((tca.nodes[f"m{i}"].x, tca.nodes[f"m{i}"].y), centre) for i in range(8)]
    assert max(radii) - min(radii) < 1 and min(radii) > 100
    # no main edges, too many reactions: no layout
    lone = DraftMap()
    lone.subsystems["L"] = Subsystem("L", "L", "pw1", "#fff", 0, 0, 100, 100)
    lone.nodes["e1"] = Node("e1", "enzyme", "L", 0, 0, 60, 30, "e1", reactions=["E1"], status="on-map")
    assert relayout.autolayout_layout(lone) is False
    assert relayout.autolayout_layout(relayout.subsystem_map(dm, "R"), max_reactions=3) is False
    # chosen by hand only, also in a subsystem without added reactions
    assert relayout.AUTOLAYOUT not in relayout.METHODS
    assert relayout.layout_methods("autolayout+untangle") == ([relayout.AUTOLAYOUT], False, True)
    dm = relayout_map()
    rows = relayout.relayout_subsystems(dm, [relayout.CONSERVATIVE], by_hand={"T": relayout.AUTOLAYOUT})
    row = next(r for r in rows if r["subsystem"] == "T")
    assert row["layout"] == relayout.AUTOLAYOUT and row["by_hand"] and dm.subsystems["T"].layout == relayout.AUTOLAYOUT
    for k in ("S", "T"):
        results, chosen = relayout.compare_layouts(relayout.subsystem_map(relayout_map(), k))
        assert relayout.AUTOLAYOUT in results and chosen != relayout.AUTOLAYOUT
    assert sorted(results, key=relayout.layout_key)[-1] == relayout.AUTOLAYOUT
    # --layout-of takes it, --relayout does not
    parser = run.build_parser()
    run.check_args(parser, parser.parse_args(["p", "j", "s", "o", "--layout-of", "a=autolayout,b=autolayout+untangle"]))
    with pytest.raises(SystemExit):
        run.check_args(parser, parser.parse_args(["p", "j", "s", "o", "--relayout", "place,autolayout"]))
    # without the package: not among the layouts, the drawing stays
    monkeypatch.setattr(relayout, "available_methods",
                        lambda methods=relayout.METHODS, old=relayout.available_methods:
                        [m for m in old(methods) if m != relayout.AUTOLAYOUT])
    dm = relayout_map()
    before = {i: (n.x, n.y) for i, n in dm.nodes.items()}
    rows = relayout.relayout_subsystems(dm, [], by_hand={"T": relayout.AUTOLAYOUT})
    assert {r["subsystem"]: r["layout"] for r in rows} == {"S": "current", "T": "current"}
    assert {i: (n.x, n.y) for i, n in dm.nodes.items()} == before
    assert relayout.AUTOLAYOUT not in relayout.compare_layouts(relayout.subsystem_map(relayout_map(), "T"))[0]


def test_short_labels():
    assert core.short_label("D-glyceraldehyde 3-phosphate", "g3p") == "g3p"
    assert core.short_label("acetyl-CoA", "accoa") == "acetyl-CoA"  # short enough
    assert core.short_label("D-glyceraldehyde 3-phosphate") == "D-glycerald…"
    assert core.short_label("D-glyceraldehyde 3-phosphate", length=0) == "D-glyceraldehyde 3-phosphate"
    dm = DraftMap()
    dm.subsystems["A"] = Subsystem("A", "A", "pw1", "#fff", 0, 0, 1000, 600)
    dm.nodes["m1"] = Node("m1", "main", "A", 0, 0, 200, 30, "D-glyceraldehyde 3-phosphate", model_mets={"g3p"})
    dm.nodes["m2"] = Node("m2", "main", "A", 0, 0, 200, 30, "2-dehydro-3-deoxy-D-gluconate", ligand_id=None)
    dm.nodes["e1"] = Node("e1", "enzyme", "A", 0, 0, 200, 30, "a very long enzyme name")
    assert core.abbreviate_labels(dm) == 2
    assert dm.nodes["m1"].label == "g3p" and dm.nodes["m1"].title == "D-glyceraldehyde 3-phosphate"
    assert dm.nodes["m1"].w == core.text_width("g3p")
    assert dm.nodes["m2"].label == "2-dehydro-3…"
    # the identity of the metabolite is still the full name
    assert core.metabolite_key(dm.nodes["m2"]) == "2-dehydro-3-deoxy-D-gluconate"
    assert dm.nodes["e1"].label == "a very long enzyme name"


def test_untangle_layout_moves_only_added_nodes():
    pytest.importorskip("pygraphviz")
    from ecmapper.pwl.layout import crossing_nodes

    dm = relayout_map()
    # the added reaction across both map reactions: a1 above them, n1 below
    dm.nodes["a1"].x, dm.nodes["a1"].y, dm.nodes["c1"].x, dm.nodes["c1"].y = 300, -100, 340, -60
    dm.nodes["n1"].x, dm.nodes["n1"].y = 300, 400
    s = relayout.subsystem_map(dm, "S")

    def crossings(s):
        return crossing_nodes(s, {i: (n.x, n.y) for i, n in s.nodes.items() if n.kind != "cofactor"})[0]

    before = crossings(s)
    assert before >= 2
    assert relayout.untangle_layout(s)
    assert crossings(s) < before
    for i in ("m1", "m2", "m3", "m4", "e1", "e2"):
        assert (s.nodes[i].x, s.nodes[i].y) == (dm.nodes[i].x, dm.nodes[i].y)
    # the cofactor moves with its enzyme node
    assert s.nodes["c1"].x - s.nodes["a1"].x == pytest.approx(40) and s.nodes["c1"].y - s.nodes["a1"].y == pytest.approx(40)


def escher_reaction(nodes, rid, mid, sub, prod, metabolites):
    """An Escher reaction drawn as by Escher's builder: midmarker `mid`, multimarkers `sub`, `prod`."""
    segments = {f"{rid}a": dict(from_node_id=mid, to_node_id=sub, b1=None, b2=None),
                f"{rid}b": dict(from_node_id=mid, to_node_id=prod, b1=None, b2=None)}
    for met, c in metabolites.items():
        marker = prod if c > 0 else sub
        a, b = (marker, met) if c > 0 else (met, marker)
        segments[f"{rid}{met}"] = dict(from_node_id=a, to_node_id=b, b1=None, b2=None)
    return dict(bigg_id=rid, label_x=0, label_y=0, segments=segments,
                metabolites=[dict(bigg_id=nodes[m]["bigg_id"], coefficient=c) for m, c in metabolites.items()])


def test_curve_reactions():
    """Straight reaction axes and main lines, curved cofactor branches, reactions at one place side by side."""
    def met(i, x, y, primary=True):
        return dict(node_type="metabolite", x=x, y=y, bigg_id=i, label_x=x, label_y=y, node_is_primary=primary)

    nodes = {"A": met("a", 0, 0), "B": met("b", 0, 400), "C": met("c", 150, 120, False), "D": met("d", 0, 0, False),
             "E": met("e", 400, 0)}
    for i in ("m1", "m2"):
        nodes[i] = dict(node_type="midmarker", x=0, y=200)
    for i in ("s1", "p1", "s2", "p2"):
        nodes[i] = dict(node_type="multimarker", x=0, y=200)
    body = dict(nodes=nodes, text_labels={}, reactions={
        "1": escher_reaction(nodes, "R1", "m1", "s1", "p1", {"A": -1, "C": -1, "B": 1, "D": 1}),
        "2": escher_reaction(nodes, "R2", "m2", "s2", "p2", {"A": -1, "B": 1})})
    before = {i: (n["x"], n["y"]) for i, n in nodes.items() if n["node_type"] == "metabolite"}
    result = escher.curve_reactions(body, movable={"D"}, lane=25.0)
    assert result == dict(drawn_reactions=2, moved_reactions=0, curved_segments=2)
    # metabolite nodes stay, except the movable one; D goes next to the product multimarker of R1 (below it)
    assert all((nodes[i]["x"], nodes[i]["y"]) == xy for i, xy in before.items() if i != "D")
    assert nodes["D"]["y"] > nodes["p1"]["y"] and abs(nodes["D"]["x"]) > 25
    # vertical axis (A above, B below): the two reactions side by side, 25 apart; multimarkers on the axes
    xs = sorted(nodes[m]["x"] for m in ("m1", "m2"))
    assert xs[1] - xs[0] == pytest.approx(25.0) and nodes["m1"]["y"] == nodes["m2"]["y"] == 200
    for mid, sub, prod in (("m1", "s1", "p1"), ("m2", "s2", "p2")):
        assert nodes[sub]["x"] == pytest.approx(nodes[mid]["x"]) == pytest.approx(nodes[prod]["x"])
        assert nodes[sub]["y"] < nodes[mid]["y"] < nodes[prod]["y"]
    # main metabolites: straight; cofactors: curves that leave the multimarker along the axis (b1 of a product
    # segment, b2 of a substrate segment)
    segments = body["reactions"]["1"]["segments"]
    assert segments["R1a"]["b1"] is None and segments["R1b"]["b2"] is None
    assert segments["R1A"]["b1"] is None and segments["R1B"]["b1"] is None
    to_d, from_c = segments["R1D"], segments["R1C"]
    assert to_d["b1"]["x"] == pytest.approx(nodes["p1"]["x"]) and to_d["b1"]["y"] > nodes["p1"]["y"]
    assert from_c["b2"]["x"] == pytest.approx(nodes["s1"]["x"]) and from_c["b2"]["y"] < nodes["s1"]["y"]
    # labels stacked right of the lanes
    labels = [body["reactions"][r] for r in ("1", "2")]
    assert all(r["label_x"] > max(xs) for r in labels) and labels[0]["label_y"] != labels[1]["label_y"]
    # a reaction drawn in another way is left alone
    nodes["x"] = dict(node_type="midmarker", x=5, y=5)
    body["reactions"]["3"] = dict(bigg_id="R3", metabolites=[dict(bigg_id="a", coefficient=-1)],
                                  segments={"s": dict(from_node_id="A", to_node_id="x", b1=None, b2=None)})
    assert escher.curve_reactions(body)["drawn_reactions"] == 2
    assert body["reactions"]["3"]["segments"]["s"]["b1"] is None


def test_curve_reactions_smooth():
    """With smooth, main lines are curves, too; through a metabolite between two reactions without a kink."""
    def met(i, x, y):
        return dict(node_type="metabolite", x=x, y=y, bigg_id=i, label_x=x, label_y=y, node_is_primary=True)

    def chain():  # A -R1-> B -R2-> C with a corner at B
        nodes = {"A": met("a", 0, 0), "B": met("b", 300, 0), "C": met("c", 300, 300),
                 "m1": dict(node_type="midmarker", x=150, y=0), "m2": dict(node_type="midmarker", x=300, y=150)}
        for i, (x, y) in (("s1", (150, 0)), ("p1", (150, 0)), ("s2", (300, 150)), ("p2", (300, 150))):
            nodes[i] = dict(node_type="multimarker", x=x, y=y)
        return dict(nodes=nodes, text_labels={}, reactions={
            "1": escher_reaction(nodes, "R1", "m1", "s1", "p1", {"A": -1, "B": 1}),
            "2": escher_reaction(nodes, "R2", "m2", "s2", "p2", {"B": -1, "C": 1})})

    body = chain()
    assert escher.curve_reactions(body, pull=False)["curved_segments"] == 0
    assert all(s["b1"] is None for r in body["reactions"].values() for s in r["segments"].values())
    body = chain()
    assert escher.curve_reactions(body, pull=False, smooth=True)["curved_segments"] == 4
    nodes, segments = body["nodes"], {**body["reactions"]["1"]["segments"], **body["reactions"]["2"]["segments"]}
    assert all(segments[s]["b1"] and segments[s]["b2"] for s in ("R1A", "R1B", "R2B", "R2C"))
    assert segments["R1a"]["b1"] is None  # the axis stays straight
    # at B (end of R1B, start of R2B) the control points are on one line through B, on either side
    bx, by = nodes["B"]["x"], nodes["B"]["y"]
    u = (segments["R1B"]["b2"]["x"] - bx, segments["R1B"]["b2"]["y"] - by)
    v = (segments["R2B"]["b1"]["x"] - bx, segments["R2B"]["b1"]["y"] - by)
    assert u[0] * v[1] - u[1] * v[0] == pytest.approx(0.0, abs=1e-6) and u[0] * v[0] + u[1] * v[1] < 0
    # A has one segment: its control point points towards the multimarker
    a, s1 = nodes["A"], nodes["s1"]
    b = segments["R1A"]["b1"]
    assert (b["x"] - a["x"]) * (s1["x"] - a["x"]) + (b["y"] - a["y"]) * (s1["y"] - a["y"]) > 0


def test_pull_reactions():
    """A U-shaped reaction moves towards the line between its metabolites, into free space, its cofactor with it."""
    def met(i, x, y, primary=True):
        return dict(node_type="metabolite", x=x, y=y, bigg_id=i, label_x=x, label_y=y, node_is_primary=primary)

    def u_shaped(blocker=None):
        nodes = {"S": met("s", 0, 0), "P": met("p", 0, 400), "C": met("c", 450, 120, False),
                 "m": dict(node_type="midmarker", x=400, y=200), "a": dict(node_type="multimarker", x=400, y=200),
                 "b": dict(node_type="multimarker", x=400, y=200)}
        if blocker:
            nodes["X"] = met("x", *blocker)
        body = dict(nodes=nodes, text_labels={},
                    reactions={"1": escher_reaction(nodes, "R1", "m", "a", "b", {"S": -1, "C": -1, "P": 1})})
        return body, escher.curve_reactions(body, lane=60.0)

    body, result = u_shaped()
    nodes = body["nodes"]
    assert result["moved_reactions"] == 1
    # to 90 from the line x = 0, the reaction axis parallel to it; the cofactor and the label on the far side
    assert nodes["m"]["x"] == pytest.approx(90.0) and nodes["m"]["y"] == pytest.approx(200.0)
    assert nodes["a"]["x"] == pytest.approx(90.0) and nodes["a"]["y"] < 200 < nodes["b"]["y"]
    assert nodes["C"]["x"] > 90 and nodes["C"]["y"] < nodes["a"]["y"]
    assert body["reactions"]["1"]["label_x"] > 90
    # a node in the way: the reaction stops before it
    body, result = u_shaped(blocker=(90, 200))
    assert result["moved_reactions"] == 1 and body["nodes"]["m"]["x"] > 120
    body, result = u_shaped(blocker=(400, 200))  # where the reaction is, and in each place it could go
    for x in (90, 152, 214, 276, 338):
        body["nodes"][f"X{x}"] = met(f"x{x}", x, 200)
    assert escher.curve_reactions(body, lane=60.0)["moved_reactions"] == 0


def test_lane_room():
    """Enzyme nodes with several model reactions get room for them side by side, across their reactions."""
    dm = DraftMap()
    dm.subsystems["s"] = Subsystem("s", "s", "pw1", "#fff", 0, 0, 1000, 1000)
    dm.nodes["e"] = Node("e", "enzyme", "s", 100, 100, 60, 30, "e", reactions=["R1", "R2", "R3"], status="added")
    dm.nodes["f"] = Node("f", "enzyme", "s", 300, 100, 60, 30, "f", reactions=["R4"], status="added")
    dm.nodes["g"] = Node("g", "enzyme", "s", 100, 300, 60, 30, "g", reactions=["R5", "R6"], status="added")
    for i, x, y in (("a", 100, 0), ("b", 100, 200), ("c", 0, 300), ("d", 200, 300)):
        dm.nodes[i] = Node(i, "main", "s", x, y, 60, 30, i, ligand_id=1)
    dm.edges += [Edge("e", "a", "main", True, False), Edge("e", "b", "main", False, True),
                 Edge("g", "c", "main", True, False), Edge("g", "d", "main", False, True)]
    # e: metabolites above and below -> reactions next to each other; g: left and right -> above each other
    assert core.lane_room(dm) == {"e": (2 * core.LANE, 0.0), "g": (0.0, core.LANE)}
    # step 4 keeps the room free: e and f at least half their widths, the room and the gap apart
    layout.condense_subsystems(dm)
    assert dm.nodes["f"].x - dm.nodes["e"].x >= (60 + 2 * core.LANE + 60) / 2 + 22.0 - 1e-6


def test_place_node_labels():
    """A metabolite label avoids the lines at its node; short map names for references."""
    def met(i, x, y):
        return dict(node_type="metabolite", x=x, y=y, bigg_id=i, label_x=x + 15, label_y=y - 22, node_is_primary=True)

    # lines from M to the upper right and to the right: the label goes to the lower right
    nodes = {"M": met("m", 0, 0), "A": met("a", 400, -400), "B": met("b", 400, 0)}
    reactions = {"1": dict(bigg_id="R1", label_x=None, segments={
        "a": dict(from_node_id="M", to_node_id="A", b1=None, b2=None),
        "b": dict(from_node_id="M", to_node_id="B", b1=None, b2=None)})}
    body = dict(nodes=nodes, reactions=reactions, text_labels={})
    escher.place_node_labels(body)
    assert nodes["M"]["label_x"] > 0 and nodes["M"]["label_y"] > 0
    assert escher.short_name("Nucleotides and cofactors") == "Nucleotides"
    assert escher.short_name("Energy and C1 metabolism") == "Energy" and escher.short_name("Lipids") == "Lipids"


def test_place_text_labels():
    """A text label goes where it covers the fewest segments and nodes."""
    def met(i, x, y):
        return dict(node_type="metabolite", x=x, y=y, bigg_id=i, label_x=x + 15, label_y=y - 22, node_is_primary=True)

    # a horizontal line through the node at (0, 0): the label goes below or above it, not right or left
    nodes = {"L": met("l", -2000, 0), "M": met("m", 0, 0), "R": met("r", 2000, 0)}
    reactions = {"1": dict(bigg_id="R1", label_x=None, segments={
        "a": dict(from_node_id="L", to_node_id="M", b1=None, b2=None),
        "b": dict(from_node_id="M", to_node_id="R", b1=None, b2=None)})}
    body = dict(nodes=nodes, reactions=reactions, text_labels={"t": dict(text="→ Lipids", x=0, y=0)})
    escher.place_text_labels(body, {"t": (0.0, 0.0, 0.0)})
    label = body["text_labels"]["t"]
    assert (label["x"], label["y"]) in ((-20, 66), (-20, -50))


def cycle_map():
    """A subsystem with a cycle of 6 metabolites (an enzyme node between each two), a parallel enzyme node, a branch."""
    dm = DraftMap()
    dm.subsystems["c"] = Subsystem("c", "cycle", "pw1", "#fff", 0, 0, 1000, 1000)
    mets = [f"m{i}" for i in range(6)]
    for i, m in enumerate(mets):  # an irregular loop
        angle = 2 * math.pi * i / 6
        dm.nodes[m] = Node(m, "main", "c", 400 + 250 * math.cos(angle) * (1 + 0.3 * (i % 2)),
                           400 + 180 * math.sin(angle), 80, 30, m, ligand_id=i)
    for i, (a, b) in enumerate(zip(mets, mets[1:] + mets[:1])):
        e = f"e{i}"
        na, nb = dm.nodes[a], dm.nodes[b]
        dm.nodes[e] = Node(e, "enzyme", "c", (na.x + nb.x) / 2, (na.y + nb.y) / 2, 60, 30, e, reactions=[f"R{i}"],
                           status="on-map")
        dm.edges += [Edge(e, a, "main", True, False), Edge(e, b, "main", False, True)]
    dm.nodes["p"] = Node("p", "enzyme", "c", 520, 420, 60, 30, "p", reactions=["P"], status="added")  # parallel to e0
    dm.edges += [Edge("p", "m0", "main", True, False), Edge("p", "m1", "main", False, True)]
    dm.nodes["x"] = Node("x", "main", "c", 900, 400, 80, 30, "x", ligand_id=9)
    dm.nodes["ex"] = Node("ex", "enzyme", "c", 800, 400, 60, 30, "ex", reactions=["X"], status="on-map")
    dm.edges += [Edge("ex", "m0", "main", True, False), Edge("ex", "x", "main", False, True)]
    for n in dm.nodes.values():
        if n.kind == "main":
            n.status = "map"
    return dm


def test_circle_layout():
    """The cycle of main metabolites on a circle, a parallel enzyme outside it; condensing keeps the circle."""
    pytest.importorskip("networkx")
    dm = cycle_map()
    cycle, graph = relayout.main_cycle(dm)
    assert len(cycle) == 6 and set(cycle) == {f"m{i}" for i in range(6)}
    assert graph["m0"]["m1"]["enzymes"] == {"e0", "p"}  # parallel enzyme nodes give one edge
    assert relayout.circle_layout(dm)
    ring = [f"m{i}" for i in range(6)] + [f"e{i}" for i in range(6)]
    cx = sum(dm.nodes[i].x for i in ring) / len(ring)
    cy = sum(dm.nodes[i].y for i in ring) / len(ring)
    radii = [math.dist((dm.nodes[i].x, dm.nodes[i].y), (cx, cy)) for i in ring]
    assert max(radii) - min(radii) < 1e-6 * max(radii) + 1.0
    assert math.dist((dm.nodes["p"].x, dm.nodes["p"].y), (cx, cy)) > radii[0]  # outside
    assert set(dm.subsystems["c"].rigid[0]) >= set(ring) | {"p"}
    layout.condense_subsystems(dm)
    cx = sum(dm.nodes[i].x for i in ring) / len(ring)
    cy = sum(dm.nodes[i].y for i in ring) / len(ring)
    radii2 = [math.dist((dm.nodes[i].x, dm.nodes[i].y), (cx, cy)) for i in ring]
    assert radii2 == pytest.approx(radii)
    # no cycle of the BRENDA map: no circle
    dm = cycle_map()
    dm.nodes["e3"].status = "added"
    assert relayout.main_cycle(dm) is None and not relayout.circle_layout(dm)


def test_route_edges():
    """An edge through another node is routed around it, the others stay straight; the routes as one cubic each."""
    pytest.importorskip("pygraphviz")
    from ecmapper.draft_map import routes

    dm = DraftMap()
    dm.subsystems["s"] = Subsystem("s", "s", "pw1", "#fff", 0, 0, 1000, 1000)
    dm.nodes["e"] = Node("e", "enzyme", "s", 300, 100, 60, 30, "e", reactions=["R1"])
    dm.nodes["a"] = Node("a", "main", "s", 300, 500, 80, 30, "a", ligand_id=1)
    dm.nodes["b"] = Node("b", "main", "s", 300, 300, 120, 30, "b", ligand_id=2)  # in the way of e -> a
    dm.nodes["c"] = Node("c", "main", "s", 500, 100, 80, 30, "c", ligand_id=3)
    dm.edges += [Edge("e", "a", "main", False, True), Edge("e", "c", "main", True, False)]
    assert routes.route_edges(dm) == 1
    route = dm.edges[0].route
    assert dm.edges[1].route is None
    assert math.dist(route[0], (300, 100)) < math.dist(route[-1], (300, 100))  # from the enzyme node
    b = dm.nodes["b"]
    assert not any(b.x - b.w / 2 < x < b.x + b.w / 2 and b.y - b.h / 2 < y < b.y + b.h / 2
                   for x, y in routes.sample(route))
    assert routes.bezier_path(route).startswith("M ")
    # fit: points of a cubic give back its control points (nearly: the fit uses chord-length parameters)
    p0, p1, p2, p3 = (0, 0), (100, 0), (200, 100), (200, 200)
    b1, b2 = routes.fit_cubic(p0, p3, routes.sample([p0, p1, p2, p3], 30)[1:-1])
    assert b1 == pytest.approx(p1, abs=6.0) and b2 == pytest.approx(p2, abs=6.0)


def test_read_compressed_model(model, tmp_path):  # noqa: F811
    """Models and subsystems files compressed with gzip, bzip2, xz or in a zip archive (with one model file)."""
    import bz2
    import gzip
    import lzma
    import zipfile

    import cobra

    from ecmapper.draft_map.prepare import read_model, read_subsystems

    model.reactions.get_by_id("PGI").subsystem = "Glycolysis"
    xml, js = tmp_path / "m.xml", tmp_path / "m.json"
    cobra.io.write_sbml_model(model, str(xml))
    cobra.io.save_json_model(model, str(js))
    reactions = {r.id for r in model.reactions}
    for name, compress in (("m.xml.gz", gzip.compress), ("m.xml.bz2", bz2.compress), ("m.json.xz", lzma.compress)):
        path = tmp_path / name
        path.write_bytes(compress((xml if ".xml" in name else js).read_bytes()))
        assert {r.id for r in read_model(path).reactions} == reactions
    with zipfile.ZipFile(tmp_path / "m.zip", "w") as z:
        z.write(js, "model/m.json")
        z.writestr("README.txt", "not a model")
    assert {r.id for r in read_model(tmp_path / "m.zip").reactions} == reactions
    assert read_subsystems(tmp_path / "m.zip")["PGI"] == "Glycolysis"
    (tmp_path / "s.json.gz").write_bytes(gzip.compress(js.read_bytes()))
    assert read_subsystems(tmp_path / "s.json.gz") == read_subsystems(js)
    with zipfile.ZipFile(tmp_path / "two.zip", "w") as z:
        z.write(js, "a.json")
        z.write(xml, "b.xml")
    with pytest.raises(ValueError, match="found 2"):
        read_model(tmp_path / "two.zip")
