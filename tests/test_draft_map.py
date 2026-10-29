import json
import pickle
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("numba")

from ecmapper import ECMapper
from ecmapper.brenda_map import BrendaMap, MapMatcher
from ecmapper.draft_map import fast, run
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


def test_draft_map(sources, model, tmp_path):  # noqa: F811
    """The whole pipeline on the BRENDA map with a small model."""
    brenda_map = BrendaMap.from_json(BRENDA_MAP / "BRENDA_map_reactions.json")
    matcher = MapMatcher(brenda_map, ECMapper(sources.values()), vicinity="propagate")
    data = prepare(model, matcher, {"PGI": "Glycolysis/Gluconeogenesis"})
    assert data["reactions"]["PGI"]["subsystem"] == "Glycolysis/Gluconeogenesis"
    assert data["reactions"]["ARGAGMt"]["transport"]
    on_map = {rid for rid, m in data["matches"].items() if m["status"] == "on map"}
    assert {"PYK", "PGI"} <= on_map
    prep = tmp_path / "prep.pkl"
    with open(prep, "wb") as f:
        pickle.dump(data, f)
    out = tmp_path / "out"
    run.main([str(prep), str(BRENDA_MAP / "BRENDA_map_reactions.json"), str(BRENDA_MAP / "BRENDA_map.svg"),
              str(out), "--iterations", "300"])

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
