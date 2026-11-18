import json
import math

import pytest

from ecmapper.brenda_escher import (Met, Rxn, align, build, circle_rings, drawn_sides, load, main_segment, place_axes,
                                    quality, tangents)
from ecmapper.brenda_escher import main as escher_main


def export():
    """A small map export: a -> b at a corner (with a cofactor couple), b -> c + d."""
    def node(i, x, y, cls="nodeHM"):
        return dict(node_id=i, display_name=i.upper(), brenda_name=i, brenda_ligand_id=ord(i[0]), node_class=cls,
                    subsystem="s", x=x, y=y, width=40, height=30)

    def entry(i):
        return dict(node_id=i, display_name=i.upper(), brenda_name=i, brenda_ligand_id=ord(i[0]))

    def reaction(i, x, y, mains, consumed=(), produced=()):
        return dict(node_id=i, subsystem="s", ec_label="1.1.1." + i[-1], enzyme_name=i,
                    layout=dict(x=x, y=y, width=60, height=30), main_metabolites=[entry(m) for m in mains],
                    cofactors=dict(consumed=[entry(m) for m in consumed], produced=[entry(m) for m in produced],
                                   undirected=[]))

    nodes = [node("a", 0, 0), node("b", 200, 200), node("c", 500, 150), node("d", 500, 250),
             node("x", 150, -60, "nodeCM"), node("y", 260, 60, "nodeCM")]
    return dict(subsystems=[dict(id="s", label_text="test", bbox=dict(x=0, y=-80, width=500, height=350),
                                 title_node=dict(x=250, y=-80))],
                reactions=[reaction("r1", 200, 0, ["a", "b"], ["x"], ["y"]), reaction("r2", 350, 200, ["b", "c", "d"])],
                metabolite_nodes=nodes)


def met(node, x, y):
    return Met(node, node, node, x, y)


def rxn(rid, x, y, subs, prods, consumed=(), produced=()):
    return Rxn(rid, "s", rid, rid, x, y, {-1: list(subs), 1: list(prods)}, {-1: list(consumed), 1: list(produced)})


def test_straight_reaction_stays_straight():
    r = rxn("r", 300, 0, [met("a", 0, 0)], [met("b", 600, 0)])
    assert place_axes([r]) == 0
    assert r.mid == (300, 0) and r.axis == (1.0, 0.0)
    assert r.markers == {-1: (270.0, 0.0), 1: (330.0, 0.0)}
    t = tangents([r])
    assert main_segment((0, 0), t["r", "a"], r.markers[-1], (-1.0, 0.0)) is None


def test_corner_is_rounded():
    # the enzyme node at the corner: substrate to its left, product below it
    r = rxn("r", 300, 0, [met("a", 0, 0)], [met("b", 300, 300)])
    assert place_axes([r]) == 1
    assert r.mid == pytest.approx((225, 75))  # a quarter of the leg along each leg
    assert r.axis == pytest.approx((math.sqrt(0.5), math.sqrt(0.5)))
    t = tangents([r])
    c_met, c_marker = main_segment((0, 0), t["r", "a"], r.markers[-1], (-r.axis[0], -r.axis[1]))
    assert c_met[1] == pytest.approx(0) and c_met[0] > 0  # leaves the substrate along the drawn link
    kx, ky = r.markers[-1]
    assert (c_marker[0] - kx) * r.axis[1] - (c_marker[1] - ky) * r.axis[0] == pytest.approx(0)  # arrives along the axis


def test_midmarker_keeps_clear_of_other_reactions():
    straight = rxn("r1", 225, 75, [met("c", 225, -300)], [met("d", 225, 450)])
    corner = rxn("r2", 300, 0, [met("a", 0, 0)], [met("b", 300, 300)])
    place_axes([straight, corner])
    assert corner.mid == (300, 0)


def test_pathway_runs_through_a_metabolite():
    # b is the corner of the pathway a -> b -> c, on enzyme nodes beside the lines between the metabolites
    def pathway():
        return [rxn("r1", 250, 150, [met("a", 0, 0)], [met("b", 400, 400)]),
                rxn("r2", 250, 650, [met("b", 400, 400)], [met("c", 0, 800)])]

    reactions = pathway()
    place_axes(reactions)
    t = tangents(reactions)
    assert t["r1", "b"] == pytest.approx((0.0, -1.0)) and t["r2", "b"] == pytest.approx((0.0, 1.0))
    assert t["r1", "a"] == pytest.approx((250 / math.hypot(250, 150), 150 / math.hypot(250, 150)))
    # a tip keeps its directions
    reactions = [rxn("r1", 300, 0, [met("a", 0, 0)], [met("b", 600, 0)]),
                 rxn("r3", 300, 60, [met("b", 600, 0)], [met("c", 0, 120)])]
    place_axes(reactions)
    assert tangents(reactions)["r1", "b"] == (-1.0, 0.0)


def test_straight_and_upright_lines_are_kept():
    # two straight reactions at a right angle: the corner stays at b
    reactions = [rxn("r1", 300, 0, [met("a", 0, 0)], [met("b", 600, 0)]),
                 rxn("r2", 600, 300, [met("b", 600, 0)], [met("c", 600, 600)])]
    place_axes(reactions)
    t = tangents(reactions)
    assert t["r1", "b"] == (-1.0, 0.0) and t["r2", "b"] == (0.0, 1.0)
    # a straight reaction, then one with its enzyme node beside the line: it leaves b along the straight one
    reactions = [rxn("r1", 300, 0, [met("a", 0, 0)], [met("b", 600, 0)]),
                 rxn("r2", 850, 100, [met("b", 600, 0)], [met("c", 1200, 0)])]
    place_axes(reactions)
    assert reactions[0].straight == {"a", "b"} and not reactions[1].straight
    t = tangents(reactions)
    assert t["r1", "b"] == (-1.0, 0.0) and t["r2", "b"] == (1.0, 0.0)


def test_nearly_straight_reaction_is_straightened():
    r = rxn("r", 300, 12, [met("a", 0, 0)], [met("b", 600, 0)])
    assert place_axes([r]) == 0
    assert r.mid == pytest.approx((300, 0)) and r.axis == (1.0, 0.0)
    # the trunk is the pair of metabolites that are most nearly opposite: c joins it from above
    r = rxn("r", 300, 0, [met("a", 0, 0), met("c", 300, -300)], [met("b", 600, 0)])
    place_axes([r])
    assert r.straight == {"a", "b"} and r.axis == (1.0, 0.0)


def test_align():
    reactions = [rxn("r1", 300, 4, [met("a", 0, 0)], [met("b", 600, 8)]),
                 rxn("r2", 604, 300, [met("b", 600, 8)], [met("c", 596, 600)])]
    assert align(reactions) == 5
    r1, r2 = reactions
    assert len({r1.y, r1.main[-1][0].y, r1.main[1][0].y, r2.main[-1][0].y}) == 1  # a row
    assert len({r2.x, r2.main[-1][0].x, r2.main[1][0].x, r1.main[1][0].x}) == 1  # a column
    assert r2.main[1][0].y == 600 and r1.main[-1][0].x == 0


def test_cycle_on_a_circle():
    # a cycle of six metabolites, drawn as a flat hexagon, with two reactions between n0 and n1
    corners = [(0, -350), (500, -150), (520, 150), (0, 350), (-520, 150), (-500, -150)]
    names = [f"n{k}" for k in range(6)]
    reactions = []
    for k in range(6):
        (ax, ay), (bx, by) = corners[k], corners[(k + 1) % 6]
        reactions.append(rxn(f"r{k}", (ax + bx) / 2, (ay + by) / 2, [met(names[k], ax, ay)],
                             [met(names[(k + 1) % 6], bx, by)]))
    reactions.append(rxn("r6", 300, -350, [met("n0", 0, -350)], [met("n1", 500, -150)]))
    circles = circle_rings(reactions, [names], least=6)
    assert len(circles) == 1
    radius = circles[0]["radius"]
    for r in reactions:
        for m in r.main[-1] + r.main[1]:
            assert math.hypot(m.x, m.y) == pytest.approx(radius)
    place_axes(reactions, circles=circles)
    assert math.hypot(*reactions[2].mid) == pytest.approx(radius)  # an arc of the circle
    inner, outer = sorted(math.hypot(*r.mid) for r in (reactions[0], reactions[6]))
    assert outer - inner == pytest.approx(180) and (inner + outer) / 2 == pytest.approx(radius)
    assert math.hypot(*reactions[6].mid) == outer  # its enzyme node is outside
    t = tangents(reactions)
    x, y = reactions[2].main[-1][0].x, reactions[2].main[-1][0].y
    assert t["r2", "n2"][0] * x + t["r2", "n2"][1] * y == pytest.approx(0)  # along the circle
    assert t["r1", "n2"] == pytest.approx((-t["r2", "n2"][0], -t["r2", "n2"][1]))


def test_drawn_sides():
    mets = [met("a", -100, -10), met("b", -100, 10), met("c", 100, 0)]
    one, other = drawn_sides(0, 0, mets)
    assert sorted([sorted(m.node for m in one), sorted(m.node for m in other)]) == [["a", "b"], ["c"]]


def test_build_map(tmp_path):
    path = tmp_path / "map.json"
    path.write_text(json.dumps(export()), encoding="utf-8")
    reactions, data = load(path)
    assert [[m.node for m in reactions[1].main[side]] for side in (-1, 1)] == [["b"], ["c", "d"]]
    assert [m.node for m in reactions[0].main[-1]] == ["a"]  # the consumed cofactor is drawn near a
    result = build(reactions, data, "test")
    body = result[1]
    assert len(body["reactions"]) == len(data["reactions"])
    nodes = body["nodes"]
    for r in body["reactions"].values():
        types = sorted(nodes[i]["node_type"] for s in r["segments"].values()
                       for i in (s["from_node_id"], s["to_node_id"]))
        assert types.count("midmarker") == 2 and types.count("multimarker") == 2 + len(r["segments"]) - 2
        assert all((s["b1"] is None) == (s["b2"] is None) for s in r["segments"].values())
    q = quality(result)
    assert q["close_nodes"] == 0
    assert escher_main([str(path), str(tmp_path), "--whole"]) == 0
    report = json.loads((tmp_path / "quality.json").read_text())
    assert "BRENDA_map" in report and (tmp_path / "BRENDA_map.escher.json").exists()



def test_svg_and_ids(tmp_path):
    """svg_text draws an Escher map in Escher's look; build records the objects of the Escher IDs with `new_id`."""
    import xml.etree.ElementTree as ET

    from ecmapper.brenda_escher import svg_text
    from ecmapper.draft_map.escher import _Ids

    path = tmp_path / "map.json"
    path.write_text(json.dumps(export()), encoding="utf-8")
    reactions, data = load(path)
    for r in reactions:  # as from_draft sets them
        r.key = ("enz_" + r.id, r.label)
        for side in (-1, 1):
            for m in r.main[side] + r.cof[side]:
                m.key = ("node", m.node, m.token)
    record = {}
    result = build(reactions, data, "test", new_id=_Ids(record=record))
    body = result[1]
    for i, r in body["reactions"].items():
        assert record[i][0] == "reaction"
    mids = [i for i, n in body["nodes"].items() if n["node_type"] == "midmarker"]
    assert sorted(record[i][0] for i in mids) == ["mid"] * len(mids)
    assert all(record[i][0] == "node" for i, n in body["nodes"].items()
               if n["node_type"] == "metabolite" and n["node_is_primary"])
    assert [record[i] for i in body["text_labels"]] == [["title", "s"]]
    assert len({tuple(key) for key in record.values()}) == len(record)  # each object once
    svg, (x, y, w, h) = svg_text(result)
    root = ET.fromstring(svg)
    ns = "{http://www.w3.org/2000/svg}"
    segments = sum(len(r["segments"]) for r in body["reactions"].values())
    metabolites = sum(n["node_type"] == "metabolite" for n in body["nodes"].values())
    assert len(root.findall(f".//{ns}path")) == segments
    assert len(root.findall(f".//{ns}g[@class='met']/{ns}circle")) == metabolites
    assert root.findall(f".//{ns}polygon")  # arrowheads
    for n in body["nodes"].values():
        assert x <= n["x"] <= x + w and y <= n["y"] <= y + h

PWL = """pwl 2.0;
subsystem demo 'demo' {
  flow down;
  compound A 'compound A';
  A -[r1: 1.1.1.1 | NAD>NADH ^R]-> B
    -L[r2: 2.2.2.2]-> C
    -R[r3: 3.3.3.3 || r4: 4.4.4.4]- D;
  D -[r6: 5.5.5.5]-> (F + G);
  react [r7: 6.6.6.6] {G, H};
  ring cw {
    P -[q1: 1.2.3.4]-> Q -R[q2: 1.2.3.5]-> R -R[q3: 1.2.3.6]-> S -R[q4: 1.2.3.7]-> T -R[q5: 1.2.3.8]-> P;
  }
}
"""


def test_maps_from_pwl(tmp_path):
    pytest.importorskip("networkx")
    pytest.importorskip("numpy")
    from ecmapper.brenda_escher import from_pwl, pwl_directory
    from ecmapper.pwl import parse

    doc = parse(PWL)
    log = {}
    reactions, data = from_pwl(doc, methods=["hints"], log=log)
    assert log == {"layout": "hints"}
    by_id = {r.id: r for r in reactions}
    r1 = by_id["r1"]
    assert [m.token for m in r1.main[-1]] == ["A"] and [m.token for m in r1.main[1]] == ["B"] and r1.oriented
    assert r1.main[-1][0].name == "compound A"
    assert [m.token for m in r1.cof[-1]] == ["NAD"] and [m.token for m in r1.cof[1]] == ["NADH"]
    assert sorted(m.token for side in by_id["r7"].main.values() for m in side) == ["G", "H"]  # sides of the layout
    assert [sorted(c["nodes"]) for c in data["circles"]] == [["P", "Q", "R", "S", "T"]]
    result = build(reactions, data, "demo")
    body = result[1]
    assert len(body["reactions"]) == 11
    assert [t["text"] for t in body["text_labels"].values()] == ["demo"]
    centre, radius = data["circles"][0]["centre"], data["circles"][0]["radius"]
    assert all(math.dist(by_id[q].mid, centre) == pytest.approx(radius) for q in ("q1", "q2", "q3", "q4", "q5"))
    (tmp_path / "demo.pwl").write_text(PWL, encoding="utf-8")
    report = pwl_directory(tmp_path, tmp_path / "out", methods=["hints"], log=lambda text: None)
    assert report["demo"]["escher"]["reactions"] == 11 and (tmp_path / "out" / "demo.escher.json").exists()
