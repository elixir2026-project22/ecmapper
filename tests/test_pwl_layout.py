import math

import numpy as np
import pytest

pytest.importorskip("networkx")

from ecmapper.pwl import parse  # noqa: E402
from ecmapper.pwl.compare import procrustes, similarity  # noqa: E402
from ecmapper.pwl.layout import (  # noqa: E402
    STEP,
    choose,
    hint_positions,
    layout_candidates,
    place_cofactors,
    remove_overlaps,
    straighten,
    to_draftmap,
)

DOC = """pwl 2.0;
subsystem demo 'demo' {
  flow down;
  compound A 'compound A';
  A -[r1: 1.1.1.1 | NAD>NADH ^R]-> B
    -L[r2: 2.2.2.2]-> C
    -R[r3: 3.3.3.3 || r4: 4.4.4.4]- D
    -B[r5: ~?]- E;
  D -[r6: 5.5.5.5]-> (F + G);
  react [r7: 6.6.6.6] {G, H};
}
"""


def _pos():
    doc = parse(DOC)
    dm = to_draftmap(doc)
    return doc, dm, hint_positions(doc, dm)


def test_to_draftmap():
    doc, dm, _ = _pos()
    kinds = [n.kind for n in dm.nodes.values()]
    assert kinds.count("enzyme") == 7 and kinds.count("main") == 8 and kinds.count("cofactor") == 2
    assert dm.nodes["A"].label == "compound A" and dm.nodes["r1"].label == "1.1.1.1"
    e = next(e for e in dm.edges if e.enzyme == "r1" and e.met == "A")
    assert (e.arrow_enzyme, e.arrow_met) == (True, False)  # substrate
    e = next(e for e in dm.edges if e.enzyme == "r3" and e.met == "D")
    assert (e.arrow_enzyme, e.arrow_met) == (True, True)  # direction unknown: double-headed
    cof = [e for e in dm.edges if e.kind == "cofactor"]
    assert [(e.arrow_enzyme, e.arrow_met) for e in cof] == [(True, False), (False, True)]


def test_turtle_follows_turn_letters():
    _, _, pos = _pos()
    d = lambda a, b: (round((pos[b][0] - pos[a][0]) / STEP), round((pos[b][1] - pos[a][1]) / STEP))  # noqa: E731
    assert d("A", "B") == (0, 1)  # flow down
    assert d("B", "C") == (1, 0)  # L: counter-clockwise on screen, from down to right
    assert d("C", "D") == (0, 1)  # R: clockwise, from right to down
    # B: back (up) from D; the place is taken by C, so E goes beside it
    assert pos["E"][1] < pos["D"][1] and pos["E"] != pos["C"]
    # the group (F + G) one step below D, side by side
    assert pos["F"][1] == pos["G"][1] and pos["F"][1] > pos["D"][1] and pos["F"][0] != pos["G"][0]
    # enzymes between their metabolites; the two reactions of the fan side by side
    assert pos["r1"] == pytest.approx(((pos["A"][0] + pos["B"][0]) / 2, (pos["A"][1] + pos["B"][1]) / 2))
    assert pos["r3"][1] == pytest.approx(pos["r4"][1]) and abs(pos["r3"][0] - pos["r4"][0]) > 60


def test_cofactors_on_the_given_side():
    doc, dm, pos = _pos()
    place_cofactors(dm, pos)
    nad, nadh = pos["r1~0"], pos["r1~1"]
    # reaction A -> B points down; ^R is to the right of the axis on screen: smaller x
    assert nad[0] < pos["r1"][0] and nadh[0] < pos["r1"][0]
    assert nad[1] < pos["r1"][1] < nadh[1]  # consumed upstream, produced downstream


def test_remove_overlaps_and_straighten():
    doc, dm, pos = _pos()
    pos = dict(pos, B=(pos["A"][0] + 5, pos["A"][1] + 3))
    pos = remove_overlaps(dm, pos)
    a, b = dm.nodes["A"], dm.nodes["B"]
    assert abs(pos["A"][0] - pos["B"][0]) >= (a.w + b.w) / 2 or abs(pos["A"][1] - pos["B"][1]) >= (a.h + b.h) / 2
    tilted = {"r1": (0.0, 0.0), "A": (20.0, -150.0)}
    assert straighten(dm, tilted)["A"][0] == straighten(dm, tilted)["r1"][0]


def test_procrustes_and_similarity_are_invariant():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(12, 2))
    rot = np.array([[0, -1], [1, 0]]) @ np.array([[-1, 0], [0, 1]])  # rotation and reflection
    assert procrustes(a, 3 * a @ rot.T + 5) == pytest.approx(0, abs=1e-9)
    assert procrustes(a, rng.normal(size=(12, 2))) > 0.3
    doc, dm, pos = _pos()
    place_cofactors(dm, pos)
    for i, (x, y) in pos.items():
        dm.nodes[i].x, dm.nodes[i].y = x, y
    import copy

    other = copy.deepcopy(dm)
    for n in other.nodes.values():
        n.x, n.y = -2 * n.y, 2 * n.x  # rotated by 90 degrees and scaled
    s = similarity(other, dm)
    assert s["procrustes"] == pytest.approx(0, abs=1e-6) and s["directions"] == 1 and s["neighbours"] == 1


def test_candidates_and_choice():
    doc = parse(DOC)
    candidates = layout_candidates(doc, methods=("hints", "hints+spring", "spring", "kamada_kawai"))
    assert set(candidates) == {"hints", "hints+spring", "spring", "kamada_kawai"}
    for dm in candidates.values():
        s = dm.subsystems["demo"]
        assert all(s.x <= n.x - n.w / 2 and n.x + n.w / 2 <= s.x + s.w for n in dm.nodes.values())
        lengths = sorted(math.dist((dm.nodes[e.enzyme].x, dm.nodes[e.enzyme].y), (dm.nodes[e.met].x, dm.nodes[e.met].y))
                         for e in dm.edges if e.kind == "main")
        assert 0.3 * STEP < lengths[len(lengths) // 2] < STEP
    assert choose(candidates) in candidates


def test_graphviz_candidates():
    pytest.importorskip("pygraphviz")
    doc = parse(DOC)
    candidates = layout_candidates(doc, methods=("hints+neato", "dot", "neato"))
    assert set(candidates) == {"hints+neato", "dot", "neato"}


def test_layout_directory(tmp_path):
    pytest.importorskip("pygraphviz")
    pytest.importorskip("scipy")
    from ecmapper.pwl.__main__ import main

    src = tmp_path / "pwl" / "subsystems"
    src.mkdir(parents=True)
    (src / "demo.pwl").write_text(DOC)
    out = tmp_path / "out"
    assert main(["layout", str(tmp_path / "pwl"), str(out), "--methods", "hints,neato"]) == 0
    svg = (out / "layout" / "demo.svg").read_text()
    assert "nodeEN" in svg and "compound A" in svg


def test_untangle_does_not_add_crossings():
    pytest.importorskip("pygraphviz")
    from ecmapper.pwl.layout import crossing_nodes, untangle

    doc, dm, pos = _pos()
    tangled = dict(pos, A=pos["C"], C=pos["A"])  # swap two metabolites: their links cross others
    before = crossing_nodes(dm, tangled)[0]
    assert before > 0
    after = crossing_nodes(dm, untangle(dm, tangled))[0]
    assert after <= before


RING = """pwl 2.0;
subsystem cyc 'cycle' {
  flow right;
  ring cw {
    A -[r1: 1.1.1.1]-> B -R[r2: 1.1.1.2]-> C -[r3: 1.1.1.3]-> D -R[r4: 1.1.1.4]-> E -[r5: 1.1.1.5]-> F
      -R[r6: 1.1.1.6]-> A;
    chord (X + C) -[r7: 2.2.2.1]-> F;
    chord (X + A) -[r8: 2.2.2.2]-> B2;
    chord E -[r9: 2.2.2.3]-> M;
    chord D -[r10: 2.2.2.4]-> M;
  }
}
"""


def test_side_metabolites_and_ring_on_a_circle():
    from ecmapper.pwl.layout import place_side_metabolites, side_metabolites

    doc = parse(RING)
    dm = to_draftmap(doc)
    side = side_metabolites(doc, dm)
    assert side == {"X"}  # only in groups, linked to two reactions
    dm.side = side
    pos = hint_positions(doc, dm)
    assert "X" not in pos
    cycle = ["A", "B", "C", "D", "E", "F"]
    cx = sum(pos[i][0] for i in cycle) / 6
    cy = sum(pos[i][1] for i in cycle) / 6
    radii = [math.dist(pos[i], (cx, cy)) for i in cycle]
    assert max(radii) - min(radii) < 1  # on a circle
    # M reacts with the neighbours D and E: outside the circle; X is placed near its enzymes afterwards
    assert math.dist(pos["M"], (cx, cy)) > max(radii)
    pos = place_side_metabolites(dm, pos, side)
    assert math.dist(pos["X"], ((pos["r7"][0] + pos["r8"][0]) / 2, (pos["r7"][1] + pos["r8"][1]) / 2)) < 1.5 * STEP


def test_shared_cofactors():
    from ecmapper.pwl.layout import share_cofactors

    doc = parse("""pwl 2.0;
subsystem s 's' {
  A -[r1: 1.1.1.1 | NAD>NADH || r2: 1.1.1.2 | NAD>NADH]-> B -[r3: 1.1.1.3 | ATP>ADP]-> C;
}
""")
    dm = to_draftmap(doc)
    pos = hint_positions(doc, dm)
    place_cofactors(dm, pos)
    before = sum(n.kind == "cofactor" for n in dm.nodes.values())
    pos = share_cofactors(dm, pos)
    labels = sorted(n.label for n in dm.nodes.values() if n.kind == "cofactor")
    assert before == 6 and labels == ["ADP", "ATP", "NAD", "NADH"]  # the fan shares NAD and NADH
    nad = next(i for i, n in dm.nodes.items() if n.label == "NAD")
    assert {e.enzyme for e in dm.edges if e.met == nad} == {"r1", "r2"}
    assert set(pos) >= {i for i in dm.nodes}
