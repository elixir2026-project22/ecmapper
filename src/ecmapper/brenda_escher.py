"""Escher maps of the subsystems of a BRENDA map, with curved reactions.

The BRENDA map draws a reaction as an enzyme box with straight links to its metabolite boxes; 4 of 10 enzyme boxes
sit at a corner or beside the line between their metabolites. Escher draws a reaction as one smooth line from the
substrates through a midmarker to the products, with the cofactors on branches. :func:`build` turns the first into
the second (see BRENDA_ESCHER.md):

1. the sides of each reaction: from its PWL description if there is one (:mod:`ecmapper.pwl`), else from the
   drawing (:func:`drawn_sides`);
2. rows and columns: nodes that are almost in a row or a column are put into it (:func:`align`), and the large
   cycles of the drawing on circles (:func:`circle_rings`);
3. the midmarker and the axis: a reaction with its metabolites opposite each other is a straight line, one on a
   circle an arc of it, and otherwise the corner at the enzyme box is rounded off (:func:`place_axes`); the
   tangents at the main metabolites: a pathway runs through a metabolite without a kink, but straight and
   upright lines are not bent for it (:func:`tangents`);
4. the main segments: cubic Bezier curves with these tangents (:func:`main_segment`);
5. the cofactors: beside the multimarker of their side, on the side of the axis where BRENDA draws them, where
   that is free (:func:`place_cofactors`);
6. the labels, where they cover the least (as in :mod:`ecmapper.draft_map.escher`).

Command line: ``python -m ecmapper.brenda_escher BRENDA_map_reactions.json OUT --pwl BRENDAmap/pwl``.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"


@dataclass
class Met:
    """A metabolite of a reaction: its drawn node (`node`: its ID in the export, `x`, `y`: its centre)."""

    node: str
    token: str  # the identifier on the Escher map
    name: str
    x: float
    y: float
    coef: float = 1.0


@dataclass
class Rxn:
    """A reaction of the map (an enzyme node). `main` and `cof`: {-1: substrates, 1: products}."""

    id: str
    subsystem: str
    label: str
    name: str
    x: float
    y: float
    main: dict[int, list[Met]]
    cof: dict[int, list[Met]]
    reversible: bool = True
    oriented: bool = False  # whether the direction is known (from PWL)
    # set by place_axes
    mid: tuple[float, float] = (0.0, 0.0)
    axis: tuple[float, float] = (1.0, 0.0)
    markers: dict[int, tuple[float, float]] = field(default_factory=dict)
    straight: frozenset = frozenset()  # the nodes of the two metabolites that it joins in a straight line
    ring: dict = field(default_factory=dict)  # on a cycle: {node: the direction in which it leaves the node}
    ring_handle: float = 0.0  # on a cycle: the length of the handles of its arcs, as a part of their chords


def _unit(x: float, y: float) -> tuple[float, float] | None:
    d = math.hypot(x, y)
    return (x / d, y / d) if d > 1e-9 else None


def _centre(mets: list[Met]) -> tuple[float, float]:
    return sum(m.x for m in mets) / len(mets), sum(m.y for m in mets) / len(mets)


def _token(name: str) -> str:
    return re.sub(r"\W+", "_", name).strip("_") or "x"


def drawn_sides(x: float, y: float, mets: list[Met]) -> tuple[list[Met], list[Met]]:
    """The main metabolites of an enzyme node at (x, y) as two sides, by where they are drawn.

    Seen from the enzyme node, each side is a range of directions: the metabolites are cut at the two widest gaps
    between their directions.
    """
    if len(mets) < 3:
        return mets[:1], mets[1:]
    order = sorted(mets, key=lambda m: math.atan2(m.y - y, m.x - x))
    angles = [math.atan2(m.y - y, m.x - x) for m in order]
    n = len(order)
    gaps = sorted(range(n), key=lambda i: -((angles[(i + 1) % n] - angles[i]) % (2 * math.pi)))
    i, j = sorted(gaps[:2])  # cuts after i and after j
    return order[i + 1:j + 1], order[j + 1:] + order[:i + 1]


def load(map_json, pwl_dir=None, scale: float = 2.0, circles: bool = True) -> tuple[list[Rxn], dict]:
    """The reactions of a BRENDA map export, their coordinates times `scale`, and the export.

    With `pwl_dir` (the PWL files of the map: ``subsystems/*.pwl``, ``node_index.json``), the sides, directions,
    coefficients and compound identifiers of the reactions come from their PWL descriptions. Without it, and for
    reactions without sides there, the sides are those of the drawing (:func:`drawn_sides`) and the identifiers
    are made from the displayed names. With `circles`, the cycles of the drawing (the ``ring`` statements of
    PWL) are put on circles (:func:`circle_rings`), which are added to the export as ``circles``. A reaction
    without a direction gets the one in which its consumed
    cofactors are drawn nearer to its substrates, or else runs from left to right (top to bottom).
    """
    data = json.loads(Path(map_json).read_text(encoding="utf-8"))
    nodes = {m["node_id"]: m for m in data["metabolite_nodes"]}
    docs, node_index, rings = {}, {}, {}
    if pwl_dir is not None:
        from .pwl import parse_file

        pwl_dir = Path(pwl_dir)
        node_index = json.loads((pwl_dir / "node_index.json").read_text(encoding="utf-8"))
        for path in sorted((pwl_dir / "subsystems").glob("*.pwl")):
            doc = parse_file(path, validate=False)
            ligands = {x: token for token, c in doc.compounds.items() for x in c.xref("brenda")}
            docs[path.stem] = (doc.reaction_index(), ligands)
            rings[path.stem] = _cycles(doc, node_index.get(path.stem, {}).get)

    def met(entry, token=None, coef=None):
        n = nodes[entry["node_id"]]
        try:
            coef = float(coef) if coef else 1.0
        except ValueError:  # a symbolic coefficient (n)
            coef = 1.0
        return Met(entry["node_id"], token or _token(entry["display_name"]),
                   (entry.get("brenda_name") or entry["display_name"]).strip(), n["x"] * scale, n["y"] * scale, coef)

    reactions = []
    for r in data["reactions"]:
        x, y = r["layout"]["x"] * scale, r["layout"]["y"] * scale
        described, ligands = docs.get(r["subsystem"], ({}, {}))
        p = described.get(r["node_id"])
        entries = {m["node_id"]: m for m in r["main_metabolites"]}
        sides, oriented, reversible = None, False, True
        if p is not None and p.sided:
            index = node_index.get(r["subsystem"], {})
            left, right = p.left, p.right
            if p.substrates is not None:
                left, right, oriented = p.substrates, p.products, True
                reversible = p.reversible is not False
            try:
                sides = [[met(entries[index[q.node]], q.token, q.coef) for q in side] for side in (left, right)]
            except KeyError:
                sides, oriented = None, False
            if sides and {m.node for side in sides for m in side} != set(entries):
                sides, oriented = None, False
        if sides is None:
            sides = list(drawn_sides(x, y, [met(e) for e in r["main_metabolites"]]))
        cof = {}
        for sign, key in ((-1, "consumed"), (1, "produced")):
            coefs = {}
            if p is not None:
                coefs = {q.token: q.coef for q in (p.consumed if sign < 0 else p.produced)}
            cof[sign] = []
            for e in r["cofactors"][key]:
                token = ligands.get(str(e.get("brenda_ligand_id")))
                cof[sign].append(met(e, token, coefs.get(token)))
        if not oriented and _backwards(sides, cof):
            sides.reverse()
        reactions.append(Rxn(r["node_id"], r["subsystem"], r["ec_label"], r.get("enzyme_name") or r["ec_label"],
                             x, y, {-1: sides[0], 1: sides[1]}, cof, reversible, oriented))
    align(reactions)
    data["circles"] = circle_rings(reactions, [c for cycles in rings.values() for c in cycles]) if circles else []
    if data["circles"]:  # rows and columns that the warp has bent, around the circles
        near = [r for r in reactions if any(math.dist((r.x, r.y), c["centre"]) < 3.5 * c["radius"]
                                            for c in data["circles"])]
        align(near, 60.0, 150.0, fixed={i for c in data["circles"] for i in c["nodes"]})
    return reactions, data


def circle_rings(reactions: list[Rxn], rings: list[list[str]], upright: float = 0.5, uneven: float = 0.35,
                 least: int = 7) -> list[dict]:
    """Put the metabolites of the cycles of the drawing on circles; returns the circles (nodes, centre, radius).

    The centre of the circle is the centre of the metabolites of the cycle, its radius their mean distance from
    it. The map is warped so that the polygon of the cycle becomes the circle: seen from the centre, a point
    inside the polygon keeps its direction and its share of the way to the polygon; a point outside moves as
    far as the polygon in its direction does, less with its distance, and not at all beyond two radii from the
    polygon. Everything in and around the cycle so keeps its place relative to it.

    Not for cycles with fewer than `least` metabolites unless they are drawn as circles already (the warp bends
    the straight lines around a cycle, which a small cycle is not worth), cycles drawn as rectangles (at least
    `upright` of their sides horizontal or vertical), cycles that do not run once around their centre, cycles
    with a metabolite more than `uneven` of the radius off the circle, and cycles next to a circle made before
    (the largest first).
    """
    position = {m.node: (m.x, m.y) for r in reactions for m in r.main[-1] + r.main[1]}
    items = [x for r in reactions for x in [r, *r.main[-1], *r.main[1], *r.cof[-1], *r.cof[1]]]
    circles = []
    for cycle in sorted(rings, key=lambda c: -len(c)):
        n = len(cycle)
        if n < 3 or any(i not in position for i in cycle):
            continue
        points = [position[i] for i in cycle]
        cx, cy = sum(p[0] for p in points) / n, sum(p[1] for p in points) / n
        spokes = [(p[0] - cx, p[1] - cy) for p in points]
        sides = [_unit(spokes[(k + 1) % n][0] - spokes[k][0], spokes[(k + 1) % n][1] - spokes[k][1]) for k in range(n)]
        if None in sides or sum(_upright(t) for t in sides) >= upright * n:
            continue
        turns = [math.atan2(spokes[k][0] * spokes[(k + 1) % n][1] - spokes[k][1] * spokes[(k + 1) % n][0],
                            spokes[k][0] * spokes[(k + 1) % n][0] + spokes[k][1] * spokes[(k + 1) % n][1])
                 for k in range(n)]
        if not (all(t > 0.02 for t in turns) or all(t < -0.02 for t in turns)):
            continue
        distances = [math.hypot(*v) for v in spokes]
        radius = sum(distances) / n
        off = max(abs(d - radius) for d in distances)
        if off > uneven * radius or (n < least and off > 0.05 * radius):
            continue
        reach = max(distances) + 2 * radius
        if any(math.dist((cx, cy), position[i]) < reach + c["radius"] for c in circles for i in c["nodes"]):
            continue

        def edge(x, y):
            """The distance of the polygon from the centre in the direction of (x, y)."""
            for k in range(n):
                (ax, ay), (bx, by) = spokes[k], spokes[(k + 1) % n]
                wx, wy = bx - ax, by - ay
                cross = x * wy - y * wx
                if abs(cross) < 1e-12:
                    continue
                rho, t = (ax * wy - ay * wx) / cross, (ax * y - ay * x) / cross
                if rho > 0 and -1e-9 <= t <= 1 + 1e-9:
                    return rho
            return None

        for item in items:
            dx, dy = item.x - cx, item.y - cy
            d = math.hypot(dx, dy)
            if d < 1e-9 or d > reach:
                continue
            rho = edge(dx / d, dy / d)
            if rho is None:
                continue
            new = d * radius / rho if d <= rho else d + (radius - rho) * max(0.0, 1 - (d - rho) / (2 * radius))
            item.x, item.y = cx + dx / d * new, cy + dy / d * new
        position = {m.node: (m.x, m.y) for r in reactions for m in r.main[-1] + r.main[1]}
        circles.append(dict(nodes=cycle, centre=(cx, cy), radius=radius))
    return circles


def _cycles(doc, node_of) -> list[list[str]]:
    """The cycles of the drawing of a PWL document (its ``ring`` statements) as lists of nodes; `node_of` gives
    the node of a drawn compound (``token#n``). As in :mod:`ecmapper.pwl.layout`: the nodes of the chains of a
    ring, without chords and groups."""
    from .pwl.syntax import Chain, Ring

    cycles = []
    for ring in (t for t in doc.statements if isinstance(t, Ring)):
        cycle = []
        for run in ring.statements:
            if isinstance(run, Chain) and not run.chord and not any(n.group or n.open for n in run.nodes):
                cycle += [node_of(n.participants[0].node) for n in run.nodes]
        cycle = list(dict.fromkeys(cycle))
        if len(cycle) >= 3 and None not in cycle:
            cycles.append(cycle)
    return cycles


def from_pwl(doc, node_ids: dict | None = None, scale: float = 2.0, methods=None, circles: bool = True,
             log: dict | None = None) -> tuple[list[Rxn], dict]:
    """The reactions of a PWL document, laid out from its description alone, for :func:`build`.

    The positions of the enzyme and metabolite nodes are those of the layout of :mod:`ecmapper.pwl.layout`
    (the hints of the description replayed and refined, the candidate with the fewest crossings; `methods`:
    its layouts to try, default all), times `scale`. Sides, directions, coefficients and identifiers are those
    of the description; reactions without sides (``react``) get those of the layout (:func:`drawn_sides`).
    Then as :func:`load`: rows and columns (:func:`align`, with a wider tolerance) and, with `circles`, the
    cycles on circles
    (:func:`circle_rings`; the layout draws rings as circles already, so cycles of any size count). `node_ids`
    names the drawn nodes (``node_index.json``; default: ``token#n``). `log` (a dict) gets the layout chosen.
    Needs the packages of the layout (extra ``draft-map``).
    """
    from .pwl.layout import METHODS, choose, layout_candidates

    node_ids = node_ids or {}
    candidates = layout_candidates(doc, node_ids, methods or METHODS, share=False)
    chosen = choose(candidates)
    dm = candidates[chosen]
    if log is not None:
        log["layout"] = chosen

    def name(token):
        c = doc.compounds.get(token)
        return c.label if c is not None and c.label else token

    def met(p, node):
        try:
            coef = float(p.coef) if p.coef else 1.0
        except ValueError:  # a symbolic coefficient (n)
            coef = 1.0
        return Met(node.id, p.token, name(p.token), node.x * scale, node.y * scale, coef)

    def main(p):
        return met(p, dm.nodes[node_ids.get(p.node, p.node)])

    reactions = []
    for eid, p in dm.pwl_reactions:
        enzyme = dm.nodes[eid]
        x, y = enzyme.x * scale, enzyme.y * scale
        oriented, reversible = False, True
        if p.sided:
            sides = [p.left, p.right]
            if p.substrates is not None:
                sides, oriented, reversible = [p.substrates, p.products], True, p.reversible is not False
            sides = [[main(q) for q in side] for side in sides]
        else:
            sides = list(drawn_sides(x, y, [main(q) for q in p.unsided]))
        cof, j = {}, 0
        for sign, participants in ((-1, p.consumed), (1, p.produced)):  # the order of the cofactor nodes
            cof[sign] = []
            for q in participants:
                cof[sign].append(met(q, dm.nodes[f"{eid}~{j}"]))
                j += 1
        if not oriented and _backwards(sides, cof):
            sides.reverse()
        label = p.enzyme_label if p.enzymes else ("spontaneous" if p.spontaneous else eid)
        reactions.append(Rxn(eid, doc.id, label, p.comment.split(";")[0].strip() or label, x, y,
                             {-1: sides[0], 1: sides[1]}, cof, reversible, oriented))
    align(reactions, 80.0, 240.0)  # the refined layouts move nodes out of the rows and columns of the hints
    data = dict(subsystems=[dict(id=doc.id, label_text=doc.title)], circles=[])
    if circles:
        data["circles"] = circle_rings(reactions, _cycles(doc, lambda node: node_ids.get(node, node)), least=3)
    return reactions, data


def pwl_directory(pwl_dir, outdir, subsystems=None, png: bool = False, scale: float = 2.0, methods=None,
                  circles: bool = True, labels: bool = True, log=print) -> dict:
    """Escher maps of the PWL files of `pwl_dir` (``subsystems/*.pwl`` or ``*.pwl``; ``node_index.json`` if there
    is one), laid out from the descriptions (:func:`from_pwl`): ``<id>.escher.json`` for each subsystem (those
    in `subsystems`), with `png` a preview, and ``quality.json``. Returns the counts of the maps."""
    from .pwl import parse_file

    pwl_dir, outdir = Path(pwl_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    index_path = pwl_dir / "node_index.json"
    node_index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else {}
    folder = pwl_dir / "subsystems" if (pwl_dir / "subsystems").is_dir() else pwl_dir
    report = {}
    for path in sorted(folder.glob("*.pwl")):
        if subsystems and path.stem not in subsystems:
            continue
        try:
            doc = parse_file(path)
            chosen = {}
            reactions, data = from_pwl(doc, node_index.get(doc.id, {}), scale, methods, circles, chosen)
            result = build(reactions, data, doc.id, labels=labels)
            (outdir / f"{path.stem}.escher.json").write_text(json.dumps(result), encoding="utf-8")
            report[path.stem] = dict(layout=chosen["layout"], circles=len(data["circles"]), escher=quality(result))
            if png:
                render(result, outdir / f"{path.stem}.png")
            log(f"{path.stem}: {chosen['layout']}")
        except Exception as e:  # noqa: BLE001 - report and go on with the other subsystems
            report[path.stem] = dict(error=str(e))
            log(f"{path.stem}: FAILED {e}")
    (outdir / "quality.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report


def align(reactions: list[Rxn], tolerance: float = 16.0, spread: float = 40.0, fixed=frozenset()) -> int:
    """Put enzyme nodes and main metabolites that are almost in a row or a column exactly into it.

    The links of the BRENDA map are mostly horizontal or vertical, but many of them only nearly (the centres of
    the boxes differ by a few units). The ends of a link that is less than `tolerance` off the horizontal get the
    same y, and so on along chains of such links (their mean y, unless they are more than `spread` apart); the
    same for vertical links and x. A row with one of the nodes `fixed` takes its place; one with several of them
    at different places is left. Returns the number of moved nodes.
    """
    moved = set()
    for here, across in (("x", "y"), ("y", "x")):
        parent, places = {}, {}

        def find(i):
            while parent.setdefault(i, i) != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for r in reactions:
            for m in r.main[-1] + r.main[1]:
                places.setdefault(m.node, []).append(m)
                off = abs(getattr(m, here) - getattr(r, here))
                if off < tolerance and abs(getattr(m, across) - getattr(r, across)) > 4 * off:
                    parent[find(r.id)] = find(m.node)
        rows = {}
        by_id = {r.id: r for r in reactions}
        for i in list(parent):
            rows.setdefault(find(i), []).append(i)
        for members in rows.values():
            items = [x for i in members for x in ([by_id[i]] if i in by_id else places[i])]
            values = [getattr(x, here) for x in items]
            if len(members) < 2 or max(values) - min(values) > spread or max(values) - min(values) < 1e-9:
                continue
            mean = sum(values) / len(values)
            anchors = {round(getattr(places[i][0], here), 3) for i in members if i in fixed}
            if len(anchors) > 1:
                continue
            if anchors:
                mean = anchors.pop()
            for x in items:
                setattr(x, here, mean)
            moved.update(members)
    return len(moved)


def _backwards(sides: list[list[Met]], cof: dict[int, list[Met]]) -> bool:
    """Whether a reaction without a direction reads better with its sides swapped."""
    a, b = sides
    if not a or not b:
        return not a and bool(b)
    (ax, ay), (bx, by) = _centre(a), _centre(b)
    score = 0.0
    for sign in (-1, 1):
        for c in cof[sign]:
            score += sign * (math.hypot(c.x - ax, c.y - ay) - math.hypot(c.x - bx, c.y - by))
    if abs(score) > 1.0:
        return score < 0  # consumed cofactors nearer to b, produced ones nearer to a
    return (bx - ax) + (by - ay) < 0


class _Grid:
    """Weighted points with owners: what is near a place."""

    def __init__(self, cell: float = 80.0):
        self.cell, self.cells = cell, {}

    def add(self, x: float, y: float, weight: float = 1.0, owner=None) -> None:
        self.cells.setdefault((int(x // self.cell), int(y // self.cell)), []).append((x, y, weight, owner))

    def near(self, x: float, y: float, r: float, skip=None) -> float:
        total = 0.0
        for i in range(int((x - r) // self.cell), int((x + r) // self.cell) + 1):
            for j in range(int((y - r) // self.cell), int((y + r) // self.cell) + 1):
                for px, py, weight, owner in self.cells.get((i, j), ()):
                    if (px - x) ** 2 + (py - y) ** 2 <= r * r and (skip is None or owner != skip):
                        total += weight
        return total


def _trunk(r: Rxn) -> tuple[Met | None, Met | None]:
    """The main substrate and the main product of a reaction that are drawn most nearly opposite each other."""
    best = None
    for m in r.main[-1] or [None]:
        for q in r.main[1] or [None]:
            if m is None or q is None:
                dot = 0.0
            else:
                a, b = _unit(m.x - r.x, m.y - r.y), _unit(q.x - r.x, q.y - r.y)
                dot = a[0] * b[0] + a[1] * b[1] if a and b else 0.0
            if best is None or dot < best[0]:
                best = (dot, m, q)
    return best[1], best[2]


def round_rings(reactions: list[Rxn], circles: list[dict], lane: float = 180.0) -> set:
    """Draw the reactions between neighbours on a circle (:func:`circle_rings`) as arcs; returns their IDs.

    A reaction whose trunk (:func:`_trunk`) joins two neighbours is an arc of the circle, with its midmarker in
    the middle of the arc. Several reactions between the same two metabolites are arcs through both, `lane`
    apart at their midmarkers, on both sides of the circle (moved by half a lane or a lane if a midmarker would
    lie on a metabolite), in the order of their enzyme nodes from the inside to the outside. The enzyme nodes
    are not used otherwise.
    """
    position = {m.node: (m.x, m.y) for r in reactions for m in r.main[-1] + r.main[1]}
    pieces = {}
    for c in circles:
        nodes = c["nodes"]
        if all(i in position for i in nodes):
            for k, i in enumerate(nodes):
                pieces.setdefault(frozenset((i, nodes[(k + 1) % len(nodes)])), c)
    groups = {}
    for r in reactions:
        s, p = _trunk(r)
        if s and p and frozenset((s.node, p.node)) in pieces:
            groups.setdefault((s.node, p.node) if s.node < p.node else (p.node, s.node), []).append((r, s))
    done = set()
    for (first, last), members in groups.items():
        c = pieces[frozenset((first, last))]
        a, b = position[first], position[last]
        half = math.dist(a, b) / 2
        e = _unit(b[0] - a[0], b[1] - a[1])
        centre = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        out = _unit(centre[0] - c["centre"][0], centre[1] - c["centre"][1]) or (-e[1], e[0])
        sagitta = c["radius"] - math.sqrt(max(c["radius"] ** 2 - half ** 2, 0.0))
        members.sort(key=lambda x: (x[0].x - centre[0]) * out[0] + (x[0].y - centre[1]) * out[1])
        rises = [sagitta + (i - (len(members) - 1) / 2) * lane for i in range(len(members))]
        for move in (0.0, 0.5, -0.5, 1.0, -1.0):  # no midmarker on a metabolite
            moved = [x + move * lane for x in rises]
            if not any(math.dist((centre[0] + out[0] * x, centre[1] + out[1] * x), q) < 0.5 * lane
                       for x in moved for q in position.values()):
                rises = moved
                break
        for (r, s), rise in zip(members, rises):
            r.mid = (centre[0] + out[0] * rise, centre[1] + out[1] * rise)
            r.axis = e if s.node == first else (-e[0], -e[1])
            beta = 2 * math.atan2(rise, half)  # the angle between the arc and its chord at the metabolites
            r.ring = {first: (e[0] * math.cos(beta) + out[0] * math.sin(beta),
                              e[1] * math.cos(beta) + out[1] * math.sin(beta)),
                      last: (-e[0] * math.cos(beta) + out[0] * math.sin(beta),
                             -e[1] * math.cos(beta) + out[1] * math.sin(beta))}
            r.ring_handle = 2 / 3 * math.tan(beta / 4) / math.sin(beta / 2) if abs(beta) > 1e-3 else 1 / 3
            done.add(r.id)
    return done


def place_axes(reactions: list[Rxn], radius: float = 450.0, marker: float = 30.0, clearance: float = 110.0,
               straight: float = 165.0, circles=()) -> int:
    """Set the midmarker, the axis and the multimarkers of each reaction; returns the number of rounded corners.

    The trunk of a reaction runs from the main substrate S to the main product P that are drawn most nearly
    opposite each other, seen from the enzyme node c (:func:`_trunk`); further main metabolites join it at the
    multimarkers. a and b are the directions from c to S and to P.

    - If the two links are at least `straight` degrees apart, the reaction is a straight line from S to P: the
      midmarker is the point of that line nearest to c, the axis its direction (`straight` of the reaction: the
      nodes of S and P). Not if another midmarker is that close to the line.
    - Otherwise the corner at c is rounded off like a fillet: the links are cut to the length R (the shorter one,
      at most `radius`), and the reaction follows the parabola through the two cut points with c as its control
      point. Its apex is the midmarker, c + R/4 (a + b), and its tangent there the axis, along b - a. A reaction
      with substrates and products on the same side of the node (a U) moves up to R/2 towards them. The midmarker
      moves only half as far, or stays at c, if it would come closer than `clearance` to a midmarker placed before
      or to a main metabolite (reactions that move least are placed first).

    The reactions on `circles` (:func:`circle_rings`) are placed before, as arcs (:func:`round_rings`). The
    multimarkers are on the axis, `marker` from the midmarker (less for short
    reactions).
    """
    for r in reactions:
        r.ring, r.ring_handle, r.straight = {}, 0.0, frozenset()
    on_ring = round_rings(reactions, circles) if circles else set()
    taken = _Grid()
    for r in reactions:
        if r.id in on_ring:
            taken.add(*r.mid)
    for r in reactions:
        for m in r.main[-1] + r.main[1]:
            taken.add(m.x, m.y, owner=r.id)
    plans = []
    limit = math.cos(math.radians(straight))
    for r in reactions:
        c = (r.x, r.y)
        s, p = _trunk(r)
        a = _unit(s.x - r.x, s.y - r.y) if s else None
        b = _unit(p.x - r.x, p.y - r.y) if p else None
        if a is None and b is None:
            a, b = (-1.0, 0.0), (1.0, 0.0)
        a = a or (-b[0], -b[1])
        b = b or (-a[0], -a[1])
        legs = [math.dist(c, (m.x, m.y)) for m in r.main[-1] + r.main[1]]
        size = min(min(legs, default=radius), radius)
        shift = (size / 4 * (a[0] + b[0]), size / 4 * (a[1] + b[1]))
        axis = _unit(b[0] - a[0], b[1] - a[1]) or (-a[1], a[0])
        line = None
        if s and p and a[0] * b[0] + a[1] * b[1] <= limit:
            chord = _unit(p.x - s.x, p.y - s.y)
            length = math.dist((s.x, s.y), (p.x, p.y))
            t = min(max(((r.x - s.x) * chord[0] + (r.y - s.y) * chord[1]) / length, 0.2), 0.8)
            line = ((s.x + chord[0] * t * length - r.x, s.y + chord[1] * t * length - r.y), chord,
                    frozenset((s.node, p.node)))
        plans.append((0.0 if line else math.hypot(*shift), r, shift, axis, line,
                      min(marker, 0.25 * min(legs, default=4 * marker))))
    rounded = 0
    for length, r, shift, axis, line, d in sorted(plans, key=lambda p: p[0]):
        if r.id in on_ring:
            pass
        elif line and not taken.near(r.x + line[0][0], r.y + line[0][1], 0.5 * clearance, skip=r.id):
            r.mid, r.axis, r.straight = (r.x + line[0][0], r.y + line[0][1]), line[1], line[2]
        else:
            r.axis = axis
            for f in (1.0, 0.5, 0.0):
                r.mid = (r.x + f * shift[0], r.y + f * shift[1])
                if f == 0.0 or not taken.near(*r.mid, clearance, skip=r.id):
                    break
            rounded += f > 0 and math.hypot(*shift) >= 1.0
        if r.id not in on_ring:
            taken.add(*r.mid)
        r.markers = {s: (r.mid[0] + s * r.axis[0] * d, r.mid[1] + s * r.axis[1] * d) for s in (-1, 1)}
    return rounded


def _upright(t: tuple[float, float]) -> bool:
    """Whether a direction is horizontal or vertical."""
    return min(abs(t[0]), abs(t[1])) < 0.035


def tangents(reactions: list[Rxn], through: float = 75.0, hub: float = 165.0, turn: float = 135.0) -> dict:
    """The directions in which the main segments leave their metabolites: {(reaction, node): unit vector}.

    Needs the axes (:func:`place_axes`). A segment of a straight reaction runs along its line, one of a reaction
    on a cycle along the cycle (both keep their directions below). Any other leaves
    its metabolite towards the enzyme node, as the link of the drawing does. Two segments of different reactions
    at a metabolite can get a common tangent, so that the pathway runs through the metabolite without a kink: at
    a metabolite with no other segment if they are at least `through` degrees apart (closer ones form a tip), at
    one of more reactions pairs that are at least `hub` degrees apart, the straightest pairs first. Straight and
    upright lines are not bent for it:

    - a segment of a straight reaction keeps its direction, and the other one takes it if the two are at least
      `turn` degrees apart;
    - of two segments that are at least `turn` degrees apart, a horizontal or vertical one gives the direction,
      else both turn by the same angle;
    - two segments at a corner (less than `turn` apart) get a common tangent only if none of them is horizontal
      or vertical; else the corner stays at the metabolite.
    """
    at = {}
    for r in reactions:
        for m in r.main[-1] + r.main[1]:
            fixed = m.node in r.straight or m.node in r.ring
            target = r.mid if fixed else (r.x, r.y)
            t = r.ring.get(m.node) or _unit(target[0] - m.x, target[1] - m.y) or (1.0, 0.0)
            at.setdefault(m.node, []).append((r.id, t, fixed))
    result = {}
    wide = math.cos(math.radians(turn))
    for node, ends in at.items():
        for rid, t, _ in ends:
            result[rid, node] = t
        limit = math.cos(math.radians(through if len(ends) == 2 else hub))
        pairs = sorted((a[1][0] * b[1][0] + a[1][1] * b[1][1], i, j) for i, a in enumerate(ends)
                       for j, b in enumerate(ends) if i < j and a[0] != b[0])
        used = set()
        for dot, i, j in pairs:
            if dot > limit or i in used or j in used:
                continue
            used.update((i, j))
            (ra, ta, fixed_a), (rb, tb, fixed_b) = ends[i], ends[j]
            if fixed_a and fixed_b:
                continue
            lead_a, lead_b = (fixed_a, fixed_b) if fixed_a or fixed_b else (_upright(ta), _upright(tb))
            if dot > wide:  # a corner
                if lead_a or lead_b:
                    continue
                t = _unit(ta[0] - tb[0], ta[1] - tb[1])
            elif lead_a != lead_b:
                t = ta if lead_a else (-tb[0], -tb[1])
            elif fixed_a or fixed_b:
                continue
            else:
                t = _unit(ta[0] - tb[0], ta[1] - tb[1])
            if t:
                result[ra, node], result[rb, node] = t, (-t[0], -t[1])
    return result


def main_segment(m: tuple[float, float], t: tuple[float, float], k: tuple[float, float], w: tuple[float, float],
                 handle: float = 0.4, max_lead: float = 150.0, arc: float = 0.0):
    """The control points (at the metabolite, at the multimarker) of the segment from a metabolite at `m`, which
    it leaves along `t`, to the multimarker at `k`, where it arrives against `w` (the axis, outwards).

    Each handle is `handle` of the distance long, shorter if its tangent points away from the other end, and the
    one at the multimarker at most `max_lead`: a long link keeps its course and bends only near the reaction.
    With `arc` (a part of the distance), both handles have that length: the segment is an arc of a circle. None
    if the segment is straight.
    """
    dx, dy = k[0] - m[0], k[1] - m[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return None
    ex, ey = dx / length, dy / length
    cos_m, cos_k = t[0] * ex + t[1] * ey, -(w[0] * ex + w[1] * ey)
    h_m = handle * length * max(cos_m, 0.3)
    h_k = min(handle * length * max(cos_k, 0.3), max_lead)
    if arc:
        h_m = h_k = arc * length
    c_m = (m[0] + t[0] * h_m, m[1] + t[1] * h_m)
    c_k = (k[0] + w[0] * h_k, k[1] + w[1] * h_k)
    off = max(abs((c[0] - m[0]) * ey - (c[1] - m[1]) * ex) for c in (c_m, c_k))
    return None if off < 2.0 and cos_m > 0 and cos_k > 0 else (c_m, c_k)


def bezier(p0, p1, p2, p3, n: int = 12) -> list[tuple[float, float]]:
    out = []
    for i in range(n + 1):
        t = i / n
        c = ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3)
        out.append((c[0] * p0[0] + c[1] * p1[0] + c[2] * p2[0] + c[3] * p3[0],
                    c[0] * p0[1] + c[1] * p1[1] + c[2] * p2[1] + c[3] * p3[1]))
    return out


def _polyline(p0, p3, controls, step: float = 25.0) -> list[tuple[float, float]]:
    p1, p2 = controls or (p0, p3)
    return bezier(p0, p1, p2, p3, max(4, int(math.dist(p0, p3) / step)))


def place_cofactors(reactions: list[Rxn], curves: dict, along: float = 50.0, step: float = 40.0,
                    lateral: float = 75.0, slots: int = 5) -> dict:
    """Positions and curves of the cofactors: {(reaction, side, index): ((x, y), control points)}.

    A cofactor goes beside the multimarker of its side: `along` (+ `step` for each further place) outwards along
    the axis and `lateral` to one side of it. The place is the one with the fewest nodes and lines of the map
    around it; the side of the axis where the cofactor is drawn in the BRENDA map and the places near the
    multimarker are preferred, places beyond the main metabolites of the side avoided. Its branch leaves the
    multimarker along the axis and bends towards it. `curves`: {(reaction, node): polyline} of the main segments.
    """
    taken = _Grid()
    for r in reactions:
        taken.add(*r.mid, 6.0, r.id)
        for m in r.main[-1] + r.main[1]:
            taken.add(m.x, m.y, 10.0)
            for d in (-45.0, 45.0):  # room for its label
                taken.add(m.x + d, m.y - 20, 2.0)
    for (rid, _), points in curves.items():
        for p in points[1:-1]:
            taken.add(*p, 1.5)
    placed = {}
    for r in sorted(reactions, key=lambda r: -len(r.cof[-1]) - len(r.cof[1])):
        ux, uy = r.axis
        nx, ny = -uy, ux
        for side in (-1, 1):
            kx, ky = r.markers[side]
            wx, wy = side * ux, side * uy
            reach = min((((m.x - kx) * wx + (m.y - ky) * wy) for m in r.main[side]), default=along + 2 * step)
            used = set()
            for index, c in enumerate(r.cof[side]):
                drawn = 1 if (c.x - r.x) * nx + (c.y - r.y) * ny >= 0 else -1
                best = None
                for s in (drawn, -drawn):
                    for j in range(slots):
                        if (s, j) in used:
                            continue
                        a = along + j * step
                        x, y = kx + wx * a + nx * s * lateral, ky + wy * a + ny * s * lateral
                        cost = taken.near(x, y, 48.0) + 0.6 * j + (1.0 if s != drawn else 0.0)
                        cost += 4.0 * max(0.0, (a + 60.0 - reach) / step)
                        # lines and nodes that the outer half of its branch would pass
                        branch = bezier((kx, ky), (kx + wx * 0.75 * a, ky + wy * 0.75 * a),
                                        (kx + wx * a + nx * s * 0.35 * lateral, ky + wy * a + ny * s * 0.35 * lateral),
                                        (x, y), 8)
                        cost += sum(taken.near(*p, 22.0) for p in branch[4:-1])
                        if best is None or cost < best[0]:
                            best = (cost, s, j, a, x, y)
                _, s, j, a, x, y = best
                used.add((s, j))
                b_marker = (kx + wx * 0.75 * a, ky + wy * 0.75 * a)
                b_met = (kx + wx * a + nx * s * 0.35 * lateral, ky + wy * a + ny * s * 0.35 * lateral)
                placed[r.id, side, index] = ((x, y), (b_marker, b_met))
                taken.add(x, y, 10.0)
                taken.add(x + 30 * s * nx, y + 30 * s * ny - 10, 2.0)
                for p in bezier((kx, ky), b_marker, b_met, (x, y), 6)[2:-1]:
                    taken.add(*p, 1.5)
    return placed


def build(reactions: list[Rxn], data: dict, map_name: str, subsystems=None, titles: bool = True,
          labels: bool = True, **options) -> list:
    """The Escher map of `reactions` (those of `subsystems`, if given), see the module.

    One Escher reaction per enzyme node (BiGG ID: its EC label) with a primary node for each drawn main
    metabolite and a secondary node of its own for each cofactor. `options`: the arguments of :func:`place_axes`
    (`radius`, `marker`, `clearance`, `straight`) and :func:`tangents` (`through`, `hub`, `turn`).
    """
    reactions = [r for r in reactions if subsystems is None or r.subsystem in subsystems]
    axes = {k: options[k] for k in ("radius", "marker", "clearance", "straight") if k in options}
    place_axes(reactions, circles=data.get("circles", ()), **axes)
    leaving = tangents(reactions, **{k: options[k] for k in ("through", "hub", "turn") if k in options})
    controls, curves = {}, {}
    for r in reactions:
        for side in (-1, 1):
            w = (side * r.axis[0], side * r.axis[1])
            for m in r.main[side]:
                c = main_segment((m.x, m.y), leaving[r.id, m.node], r.markers[side], w,
                                 arc=r.ring_handle if m.node in r.ring else 0.0)
                controls[r.id, m.node] = c
                curves[r.id, m.node] = _polyline((m.x, m.y), r.markers[side], c)
    cofactors = place_cofactors(reactions, curves)

    count = 0

    def new_id() -> str:
        nonlocal count
        count += 1
        return str(count)

    def point(p):
        return dict(x=p[0], y=p[1])

    nodes, out, drawn = {}, {}, {}
    for r in reactions:
        mid = new_id()
        nodes[mid] = dict(node_type="midmarker", x=r.mid[0], y=r.mid[1])
        segments, coefficients = {}, {}
        for side in (-1, 1):
            marker = new_id()
            nodes[marker] = dict(node_type="multimarker", x=r.markers[side][0], y=r.markers[side][1])
            pair = (marker, mid) if side < 0 else (mid, marker)
            segments[new_id()] = dict(from_node_id=pair[0], to_node_id=pair[1], b1=None, b2=None)
            ends = []
            for m in r.main[side]:
                if m.node not in drawn:
                    drawn[m.node] = new_id()
                    nodes[drawn[m.node]] = dict(node_type="metabolite", x=m.x, y=m.y, bigg_id=m.token, name=m.name,
                                                label_x=m.x + 20, label_y=m.y - 20, node_is_primary=True)
                ends.append((m, drawn[m.node], controls[r.id, m.node]))
            for index, m in enumerate(r.cof[side]):
                (x, y), (b_marker, b_met) = cofactors[r.id, side, index]
                node = new_id()
                nodes[node] = dict(node_type="metabolite", x=x, y=y, bigg_id=m.token, name=m.name,
                                   label_x=x + 10, label_y=y - 10, node_is_primary=False)
                ends.append((m, node, (b_met, b_marker)))
            for m, node, c in ends:
                b_met, b_marker = (point(c[0]), point(c[1])) if c else (None, None)
                if side < 0:
                    segments[new_id()] = dict(from_node_id=node, to_node_id=marker, b1=b_met, b2=b_marker)
                else:
                    segments[new_id()] = dict(from_node_id=marker, to_node_id=node, b1=b_marker, b2=b_met)
                coefficients[m.token] = coefficients.get(m.token, 0.0) + side * m.coef
        out[new_id()] = dict(
            name=r.name, bigg_id=r.label, reversibility=r.reversible, label_x=r.mid[0] + 20, label_y=r.mid[1] + 40,
            gene_reaction_rule="", genes=[],
            metabolites=[dict(bigg_id=t, coefficient=c) for t, c in coefficients.items() if c],
            segments=segments)
    text_labels = {}
    if titles:
        extent = {}  # subsystem -> left, right, top of its nodes
        for r in reactions:
            xs = [r.mid[0]] + [m.x for m in r.main[-1] + r.main[1]]
            top = min([r.mid[1]] + [m.y for m in r.main[-1] + r.main[1]])
            x0, x1, y0 = extent.get(r.subsystem, (min(xs), max(xs), top))
            extent[r.subsystem] = (min(x0, *xs), max(x1, *xs), min(y0, top))
        for s in data["subsystems"]:
            if s["id"] in extent:
                text = s.get("label_text") or s["id"]
                x0, x1, y0 = extent[s["id"]]
                text_labels[new_id()] = dict(text=text, x=(x0 + x1) / 2 - 14.0 * len(text), y=y0 - 170)
    body = dict(reactions=out, nodes=nodes, text_labels=text_labels)
    if labels:
        place_labels(body)
    xs = [n["x"] for n in nodes.values()] + [t["x"] for t in text_labels.values()]
    xs += [t["x"] + 28.0 * len(t["text"]) for t in text_labels.values()]
    ys = [n["y"] for n in nodes.values()] + [t["y"] for t in text_labels.values()]
    margin = 250
    body["canvas"] = dict(x=min(xs) - margin, y=min(ys) - margin, width=max(xs) - min(xs) + 2 * margin,
                          height=max(ys) - min(ys) + 2 * margin)
    head = dict(map_name=map_name, map_id=re.sub(r"\W+", "_", map_name), map_description="from the BRENDA map",
                homepage="https://escher.github.io", schema=SCHEMA)
    return [head, body]


def place_labels(body: dict) -> None:
    """Put the labels of the reactions and metabolites where they cover the fewest lines, nodes and labels."""
    from .draft_map.escher import ESCHER_FONTS, _choose, _map_points, place_node_labels

    nodes = body["nodes"]
    points = _map_points(body, node_labels=False)
    for rid, r in body["reactions"].items():
        mid = next(nodes[s[k]] for s in r["segments"].values() for k in ("from_node_id", "to_node_id")
                   if nodes[s[k]]["node_type"] == "midmarker")
        x, y, w = mid["x"], mid["y"], ESCHER_FONTS["reaction"][1] * len(r["bigg_id"])
        candidates = [(x + 22, y + 10), (x - 22 - w, y + 10), (x - w / 2, y + 48), (x - w / 2, y - 22),
                      (x + 18, y - 18), (x + 18, y + 40), (x - 18 - w, y - 18), (x - 18 - w, y + 40),
                      (x + 60, y + 10), (x - 60 - w, y + 10), (x - w / 2, y + 85), (x - w / 2, y - 60)]
        r["label_x"], r["label_y"] = _choose(points, ("label", "reaction", rid), "reaction", r["bigg_id"], candidates)
    place_node_labels(body)


def _crossings(lines: list[tuple[frozenset, list[tuple[float, float]]]]) -> int:
    """The number of crossings between polylines that share no end node (`lines`: (end nodes, points))."""
    cell, grid, pieces = 200.0, {}, []
    for ends, points in lines:
        for p, q in zip(points, points[1:]):
            index = len(pieces)
            pieces.append((ends, p, q))
            for i in range(int(min(p[0], q[0]) // cell), int(max(p[0], q[0]) // cell) + 1):
                for j in range(int(min(p[1], q[1]) // cell), int(max(p[1], q[1]) // cell) + 1):
                    grid.setdefault((i, j), []).append(index)

    def side(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    found = set()
    for members in grid.values():
        for x, i in enumerate(members):
            for j in members[x + 1:]:
                (e1, p1, q1), (e2, p2, q2) = pieces[i], pieces[j]
                if e1 & e2 or (i, j) in found:
                    continue
                if side(p1, q1, p2) * side(p1, q1, q2) < 0 and side(p2, q2, p1) * side(p2, q2, q1) < 0:
                    found.add((i, j))
    return len({(pieces[i][0], pieces[j][0]) for i, j in found})


def quality(data: list) -> dict:
    """Counts of an Escher map: reactions, curved segments, kinks, crossings and nodes on top of each other.

    `kinks`: primary metabolites of exactly two segments that meet at an angle of 75 to 170 degrees (straighter
    than 170: no kink; below 75: a tip). `main_crossings`: pairs of segments to primary metabolites of different
    reactions that cross; `crossings`: of all segments.
    """
    body = data[1]
    nodes = body["nodes"]
    lines, curved, ends = [], 0, {}
    for rid, r in body["reactions"].items():
        for s in r["segments"].values():
            a, b = nodes[s["from_node_id"]], nodes[s["to_node_id"]]
            p0, p3 = (a["x"], a["y"]), (b["x"], b["y"])
            c = ((s["b1"]["x"], s["b1"]["y"]), (s["b2"]["x"], s["b2"]["y"])) if s["b1"] else None
            curved += c is not None
            primary = any(n.get("node_is_primary") for n in (a, b))
            lines.append((primary, frozenset((rid, s["from_node_id"], s["to_node_id"])), _polyline(p0, p3, c, 40.0)))
            for key, n, q in ((s["from_node_id"], a, c[0] if c else p3), (s["to_node_id"], b, c[1] if c else p0)):
                if n.get("node_is_primary"):
                    ends.setdefault(key, []).append(_unit(q[0] - n["x"], q[1] - n["y"]))
    kinks = 0
    for directions in ends.values():
        if len(directions) == 2 and all(directions):
            (ax, ay), (bx, by) = directions
            angle = math.degrees(math.acos(max(-1.0, min(1.0, ax * bx + ay * by))))
            kinks += 75 <= angle < 170
    mets = [(n["x"], n["y"]) for n in nodes.values() if n["node_type"] == "metabolite"]
    grid = _Grid()
    for p in mets:
        grid.add(*p)
    return dict(reactions=len(body["reactions"]), nodes=len(nodes), curved_segments=curved, kinks=kinks,
                main_crossings=_crossings([(e, p) for primary, e, p in lines if primary]),
                crossings=_crossings([(e, p) for _, e, p in lines]),
                close_nodes=int(sum(grid.near(*p, 35.0) - 1 for p in mets) / 2))


def drawing_quality(reactions: list[Rxn]) -> dict:
    """The same counts for the drawing of the BRENDA map (straight links between the enzyme and metabolite nodes)."""
    lines, at = [], {}
    for r in reactions:
        for m in r.main[-1] + r.main[1]:
            lines.append((frozenset((r.id, m.node)), [(r.x, r.y), (m.x, m.y)]))
            at.setdefault(m.node, []).append(_unit(r.x - m.x, r.y - m.y))
    kinks = corners = 0
    for directions in at.values():
        if len(directions) == 2 and all(directions):
            (ax, ay), (bx, by) = directions
            angle = math.degrees(math.acos(max(-1.0, min(1.0, ax * bx + ay * by))))
            kinks += 75 <= angle < 170
    for r in reactions:
        if r.main[-1] and r.main[1]:
            a = _unit(*(p - q for p, q in zip(_centre(r.main[-1]), (r.x, r.y))))
            b = _unit(*(p - q for p, q in zip(_centre(r.main[1]), (r.x, r.y))))
            if a and b:
                corners += math.degrees(math.acos(max(-1.0, min(1.0, a[0] * b[0] + a[1] * b[1])))) < 170
    return dict(reactions=len(reactions), corners_at_enzymes=corners, kinks=kinks, main_crossings=_crossings(lines))


def render(data: list, path, width: float = 1800.0) -> None:
    """A preview of an Escher map as PNG (matplotlib; Escher's default colours and sizes, without arrowheads)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, PathPatch
    from matplotlib.path import Path as MPath

    body = data[1]
    nodes, canvas = body["nodes"], body["canvas"]
    px = min(width / canvas["width"], 0.6)  # pixels per map unit
    fig = plt.figure(figsize=(canvas["width"] * px / 100, canvas["height"] * px / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(canvas["x"], canvas["x"] + canvas["width"])
    ax.set_ylim(canvas["y"] + canvas["height"], canvas["y"])
    ax.axis("off")
    pt = px * 0.72  # points per map unit
    for r in body["reactions"].values():
        for s in r["segments"].values():
            a, b = nodes[s["from_node_id"]], nodes[s["to_node_id"]]
            if s["b1"]:
                verts = [(a["x"], a["y"]), (s["b1"]["x"], s["b1"]["y"]), (s["b2"]["x"], s["b2"]["y"]), (b["x"], b["y"])]
                codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
            else:
                verts, codes = [(a["x"], a["y"]), (b["x"], b["y"])], [MPath.MOVETO, MPath.LINETO]
            ax.add_patch(PathPatch(MPath(verts, codes), fc="none", ec="#334E75", lw=10 * pt, capstyle="round"))
        ax.text(r["label_x"], r["label_y"], r["bigg_id"], fontsize=30 * pt, color="#191978", style="italic",
                weight="bold", va="baseline")
    for n in nodes.values():
        if n["node_type"] == "metabolite":
            radius = 20 if n.get("node_is_primary") else 10
            ax.add_patch(Circle((n["x"], n["y"]), radius, fc="#E0865B", ec="#A24510", lw=2 * pt, zorder=3))
            ax.text(n["label_x"], n["label_y"], n["bigg_id"], fontsize=20 * pt, style="italic", weight="bold",
                    va="baseline", zorder=4)
        else:
            ax.add_patch(Circle((n["x"], n["y"]), 5, fc="white", ec="#334E75", lw=1 * pt, zorder=3))
    for t in body["text_labels"].values():
        ax.text(t["x"], t["y"], t["text"], fontsize=50 * pt, style="italic", weight="bold", va="baseline")
    fig.savefig(path)
    plt.close(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m ecmapper.brenda_escher", description=__doc__.split("\n")[0])
    parser.add_argument("map_json", type=Path, help="BRENDA_map_reactions.json")
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--pwl", type=Path, help="folder with the PWL files of the map (sides and directions)")
    parser.add_argument("--subsystem", action="append", help="only these subsystems (repeatable)")
    parser.add_argument("--whole", action="store_true", help="also write the whole map as one Escher map")
    parser.add_argument("--png", action="store_true", help="also write previews (needs matplotlib)")
    parser.add_argument("--scale", type=float, default=2.0, help="Escher units per map unit (default: 2)")
    parser.add_argument("--radius", type=float, default=450.0, help="largest rounding of a corner (default: 450)")
    parser.add_argument("--no-labels", action="store_true", help="do not place the labels (faster)")
    parser.add_argument("--no-circles", action="store_true", help="leave the cycles as they are drawn")
    args = parser.parse_args(argv)
    reactions, data = load(args.map_json, args.pwl, args.scale, circles=not args.no_circles)
    ids = [s["id"] for s in data["subsystems"]]
    selected = args.subsystem or ids
    unknown = sorted(set(selected) - set(ids))
    if unknown:
        parser.error("unknown subsystems: " + ", ".join(unknown))
    args.outdir.mkdir(parents=True, exist_ok=True)
    report = {}
    jobs = [(sid, [sid]) for sid in selected] + ([("BRENDA_map", selected)] if args.whole else [])
    for name, members in jobs:
        mine = [r for r in reactions if r.subsystem in members]
        if not mine:
            continue
        result = build(reactions, data, name, subsystems=set(members), radius=args.radius,
                       labels=not args.no_labels)
        (args.outdir / f"{name}.escher.json").write_text(json.dumps(result), encoding="utf-8")
        report[name] = dict(escher=quality(result), drawing=drawing_quality(mine))
        if args.png:
            render(result, args.outdir / f"{name}.png")
    (args.outdir / "quality.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    total = {k: sum(v["escher"][k] for n, v in report.items() if n != "BRENDA_map")
             for k in ("reactions", "curved_segments", "kinks", "main_crossings", "crossings", "close_nodes")}
    before = {k: sum(v["drawing"][k] for n, v in report.items() if n != "BRENDA_map")
              for k in ("corners_at_enzymes", "kinks", "main_crossings")}
    print(json.dumps(dict(maps=len(report), escher=total, drawing=before), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
