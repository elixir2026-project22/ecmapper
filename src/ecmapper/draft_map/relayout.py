"""Step 3b: lay out the subsystems with added reactions again (networkx and Graphviz layouts, the best one is kept).

The completion (step 3) places added reactions between their drawn
metabolites, which gives long, crossing edges where many reactions are
added. By default (:data:`CONSERVATIVE`), the BRENDA drawing stays as it
is and the added reactions are placed next to their metabolites only where
they cross nothing; the others are left out (:func:`place_conservative`),
and no other layout follows. With :data:`TEMPLATE`, only the added
reactions that go where a node of the BRENDA map with their EC number or
reaction ID was are drawn, at its place (:func:`place_template`). Otherwise,
for each subsystem with added
reactions, several layouts of its
enzyme and main metabolite nodes are computed, each subsystem is condensed
as in step 4, and the layout with the fewest crossings of main-metabolite
links is kept:

- ``place``: the added nodes placed again one by one next to their
  metabolites, at free points whose links cross least, the nodes of the
  BRENDA map fixed (:func:`place_added`); the other layouts start from it;
- ``spring``: networkx spring layout (Fruchterman-Reingold) of the added
  nodes, the nodes of the BRENDA map fixed;
- ``fdp``: the same with Graphviz ``fdp``;
- ``dot``: Graphviz ``dot`` (layered, Sugiyama), along the reaction directions;
- ``neato``: Graphviz ``neato`` (Kamada-Kawai), starting from the current
  drawing, the nodes of the BRENDA map pinned softly: each is tied to a
  pinned ghost node at its place by a short, heavy edge, so it may drift a
  little to make room (:func:`graphviz_spec`).

``place``, ``spring`` and ``fdp`` keep the drawing of the BRENDA map; ``dot``
and ``neato`` are only taken if they have at least `margin` fewer crossings,
and at least `min_gain` fewer, than the best of the current drawing,
``place``, ``spring`` and ``fdp`` (a relative margin alone lets one crossing
fewer discard the BRENDA drawing). Each layout is also tried untangled
locally (``<layout>+untangle``, :func:`untangle_layout`; in those that keep
the BRENDA drawing only added nodes move). Cofactors move with their enzyme
node. In subsystems whose BRENDA drawing has a cycle of main metabolites,
added reactions that touch at most two nodes of the BRENDA map are kept
outside the cycle (with ``place`` a rule of :func:`place_added`; without it,
they start outside (:func:`outside_cycle`), the moved drawing is the
candidate ``outside``, and ``spring``, ``fdp``, ``neato`` start from it). Subsystems
named in `circle` get the cycle of their main metabolites on a circle
(:func:`circle_layout`), whatever the crossings.
One more layout is only taken where it is chosen by hand (``--layout-of``,
the GUI), in any subsystem, also one without added reactions:
``autolayout``, the layout of the package sbml-autolayout for small pathways
(a straight main line, known cycles on circles; :func:`autolayout_layout`).
It redraws the whole subsystem.
Graphviz layouts need pygraphviz (2.0 or later ships Graphviz), networkx
layouts networkx, ``autolayout`` sbml_autolayout; missing ones are skipped.

Graphviz ``fdp`` keeps state between layouts in a process: after a first
``fdp`` layout, later ones vary between runs. Each Graphviz layout is
therefore computed in a process of its own (:mod:`.graphviz_worker`, several
in parallel), where it is reproducible.
"""

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .core import DraftMap, lane_room
from .layout import condense_subsystems, internal_metrics

#: node status of nodes of the BRENDA map (enzyme and metabolite nodes)
MAP_STATUS = {"map", "on-map", "vicinity"}
#: candidate layouts that keep the drawing of the BRENDA map, and those that do not
KEEPING = ("current", "outside", "place", "spring", "fdp")
REDRAWING = ("dot", "neato")
METHODS = tuple(m for m in KEEPING + REDRAWING if m != "outside")  # the layouts of --relayout that compete
#: the conservative placement (:func:`place_conservative`), the default of --relayout: alone, no other layout
CONSERVATIVE = "place_conservative"
#: only added reactions at the place of a node of the BRENDA map (:func:`place_template`): alone, no other layout
TEMPLATE = "place_template"
#: the layout of the package sbml-autolayout (:func:`autolayout_layout`): redraws the subsystem, only chosen by hand
AUTOLAYOUT = "autolayout"
#: suffix of a layout untangled locally (:func:`untangle_layout`)
UNTANGLED = "+untangle"


def subsystem_map(dm: DraftMap, k: str) -> DraftMap:
    """A copy of subsystem `k` (nodes and the edges inside it) as a draft map of its own."""
    s = DraftMap()
    s.subsystems = {k: copy.deepcopy(dm.subsystems[k])}
    s.nodes = {i: copy.deepcopy(n) for i, n in dm.nodes.items() if n.subsystem == k}
    s.edges = [copy.copy(e) for e in dm.edges if e.enzyme in s.nodes and e.met in s.nodes]
    return s


def _owners(s: DraftMap) -> dict[str, str]:
    """{cofactor node: its (first) enzyme node}."""
    owner = {}
    for e in s.edges:
        if e.kind == "cofactor":
            owner.setdefault(e.met, e.enzyme)
    return owner


def _groups(s: DraftMap) -> tuple[dict[str, list[float]], dict[str, str]]:
    """Boxes of the main nodes; an enzyme node with its cofactors is one box, with room for its reactions."""
    owner = _owners(s)
    room = lane_room(s)
    box = {}
    for i, n in s.nodes.items():
        if n.kind != "cofactor":
            w, h = n.w + room.get(i, (0.0, 0.0))[0], n.h + room.get(i, (0.0, 0.0))[1]
            box[i] = [n.x - w / 2, n.y - h / 2, n.x + w / 2, n.y + h / 2]
    for c, enzyme in owner.items():
        n, b = s.nodes[c], box[enzyme]
        b[0], b[1] = min(b[0], n.x - n.w / 2), min(b[1], n.y - n.h / 2)
        b[2], b[3] = max(b[2], n.x + n.w / 2), max(b[3], n.y + n.h / 2)
    return box, owner


def _move(s: DraftMap, centres: dict[str, tuple[float, float]], box, owner) -> None:
    """Move the main nodes so that their boxes are centred at `centres`; cofactors follow their enzyme node."""
    before = {i: (n.x, n.y) for i, n in s.nodes.items()}
    for i, (cx, cy) in centres.items():
        b = box[i]
        s.nodes[i].x = before[i][0] + cx - (b[0] + b[2]) / 2
        s.nodes[i].y = before[i][1] + cy - (b[1] + b[3]) / 2
    for c, enzyme in owner.items():
        s.nodes[c].x = s.nodes[enzyme].x + before[c][0] - before[enzyme][0]
        s.nodes[c].y = s.nodes[enzyme].y + before[c][1] - before[enzyme][1]


def _directed(s: DraftMap) -> list[tuple[str, str]]:
    """Main-metabolite edges along the reaction direction; reversible ones downwards as currently drawn."""
    out = []
    for e in s.edges:
        if e.kind != "main":
            continue
        enzyme, met = s.nodes[e.enzyme], s.nodes[e.met]
        if e.arrow_enzyme and not e.arrow_met:
            out.append((e.met, e.enzyme))  # substrate
        elif e.arrow_met and not e.arrow_enzyme:
            out.append((e.enzyme, e.met))  # product
        else:
            out.append((e.met, e.enzyme) if met.y <= enzyme.y else (e.enzyme, e.met))
    return out


def spring_layout(s: DraftMap, k: float = 150.0, iterations: int = 200, seed: int = 1) -> None:
    """networkx spring layout of the enzyme and main metabolite nodes, nodes of the BRENDA map fixed."""
    import networkx as nx

    owner = _owners(s)
    before = {i: (n.x, n.y) for i, n in s.nodes.items()}
    g = nx.Graph()
    g.add_nodes_from(i for i, n in s.nodes.items() if n.kind != "cofactor")
    g.add_edges_from(_directed(s))
    fixed = [i for i in g if s.nodes[i].status in MAP_STATUS]
    pos = nx.spring_layout(g, k=k, pos={i: before[i] for i in g}, fixed=fixed or None, iterations=iterations, seed=seed)
    for i, (x, y) in pos.items():
        s.nodes[i].x, s.nodes[i].y = float(x), float(y)
    for c, enzyme in owner.items():
        s.nodes[c].x = s.nodes[enzyme].x + before[c][0] - before[enzyme][0]
        s.nodes[c].y = s.nodes[enzyme].y + before[c][1] - before[enzyme][1]


def metabolite_graph(s: DraftMap, brenda: bool = False):
    """networkx graph of the main metabolites, joined if an enzyme node turns one into the other.

    Pairs of an enzyme node: its substrates with its products, and the
    metabolites of reversible edges (arrows at both ends) with all others.
    Parallel enzyme nodes give one edge; its attribute "enzymes" lists them.
    With `brenda`, only nodes of the BRENDA map (no added ones).
    """
    import itertools

    import networkx as nx

    sides = {}  # enzyme node -> ([substrates], [products], [reversible])
    for e in s.edges:
        if e.kind != "main":
            continue
        if brenda and not (s.nodes[e.enzyme].status in MAP_STATUS and s.nodes[e.met].status in MAP_STATUS):
            continue
        side = 2 if e.arrow_enzyme == e.arrow_met else 0 if e.arrow_enzyme else 1
        sides.setdefault(e.enzyme, ([], [], []))[side].append(e.met)
    g = nx.Graph()
    for enzyme, (sub, prod, rev) in sides.items():
        pairs = set(itertools.product(sub, prod)) | {(a, b) for a in rev for b in sub + prod + rev if a != b}
        for a, b in pairs:
            if g.has_edge(a, b):
                g[a][b]["enzymes"].add(enzyme)
            else:
                g.add_edge(a, b, enzymes={enzyme})
    return g


def main_cycle(s: DraftMap, min_length: int = 5, max_length: int = 14, limit: int = 100000):
    """The longest cycle of main metabolites of the BRENDA map (at least `min_length`), or None.

    Only nodes of the BRENDA map: cycles closed by added reactions (step 3)
    are mostly artefacts of reactions with several substrates. Among cycles
    of the same length, the one that encloses the largest area in the
    current drawing (the outer loop rather than one through a metabolite
    drawn inside it). Returns (metabolite nodes in cycle order, metabolite
    graph of all nodes, for the parallel reactions).
    """
    import networkx as nx

    g = metabolite_graph(s, brenda=True)
    best, key = None, None
    for i, cycle in enumerate(nx.simple_cycles(g, length_bound=max_length)):
        if i >= limit:
            break
        if len(cycle) < min_length:
            continue
        xy = [(s.nodes[m].x, s.nodes[m].y) for m in cycle]
        area = abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(xy, xy[1:] + xy[:1]))) / 2
        if key is None or (len(cycle), area) > key:
            best, key = cycle, (len(cycle), area)
    return (best, metabolite_graph(s)) if best else None


def _ray_exit(poly, cx, cy, dx, dy) -> float:
    """Distance from (cx, cy) along the unit direction (dx, dy) to the last crossing of the polygon's boundary."""
    best = 0.0
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        ex, ey = x2 - x1, y2 - y1
        den = dx * ey - dy * ex
        if abs(den) < 1e-9:
            continue
        t = ((x1 - cx) * ey - (y1 - cy) * ex) / den
        u = ((x1 - cx) * dy - (y1 - cy) * dx) / den
        if t > 0 and 0 <= u <= 1:
            best = max(best, t)
    return best


def outside_cycle(s: DraftMap, max_members: int = 2, clearance: float = 120.0) -> bool:
    """Added reactions that touch at most `max_members` nodes of the BRENDA map start outside its cycle.

    The completion (step 3) places an added reaction between its drawn
    metabolites, so a reaction between two neighbours of a cycle (or at one
    of its members) lands inside the ring and its edges cross it. With a
    cycle of main metabolites of the BRENDA map (:func:`main_cycle`), each
    added enzyme node whose main-metabolite partners of the BRENDA map number
    at most `max_members` is moved onto the ray from the centre of the cycle
    through the middle of these partners, `clearance` beyond the polygon of
    the cycle, at a free place (:func:`..core.place_free`). Its cofactors and
    its metabolite nodes that no other enzyme node uses move with it (as the
    PWL layouts put triangles outside a ring, LAYOUT_LEARNINGS.md 5). Added
    nodes joined to more BRENDA nodes (chords) stay. False if nothing moved.
    """
    import math

    from .core import place_free

    found = main_cycle(s)
    if not found:
        return False
    cycle, _ = found
    poly = [(s.nodes[m].x, s.nodes[m].y) for m in cycle]
    cx, cy = sum(p[0] for p in poly) / len(poly), sum(p[1] for p in poly) / len(poly)
    owner = _owners(s)
    partners: dict[str, set[str]] = {}  # added enzyme node -> its main metabolites of the BRENDA map
    users: dict[str, set[str]] = {}  # metabolite node -> enzyme nodes using it
    for e in s.edges:
        if e.kind != "main":
            continue
        users.setdefault(e.met, set()).add(e.enzyme)
        if s.nodes[e.enzyme].status == "added" and s.nodes[e.met].status in MAP_STATUS:
            partners.setdefault(e.enzyme, set()).add(e.met)
    moved = 0
    for enzyme, mets in partners.items():
        if len(mets) > max_members:
            continue
        en = s.nodes[enzyme]
        own = [c for c, o in owner.items() if o == enzyme]
        own += [m for m, es in users.items() if es == {enzyme} and s.nodes[m].status not in MAP_STATUS]
        px = sum(s.nodes[m].x for m in mets) / len(mets)
        py = sum(s.nodes[m].y for m in mets) / len(mets)
        d = math.hypot(px - cx, py - cy) or 1.0
        dx, dy = (px - cx) / d, (py - cy) / d
        r = _ray_exit(poly, cx, cy, dx, dy) + clearance
        boxes = [n.box for i, n in s.nodes.items() if i != enzyme and i not in own]
        x, y = place_free(cx + r * dx, cy + r * dy, en.w, en.h, boxes)
        ddx, ddy = x - en.x, y - en.y
        en.x, en.y = x, y
        boxes.append(en.box)
        for i in own:
            n = s.nodes[i]
            n.x, n.y = place_free(n.x + ddx, n.y + ddy, n.w, n.h, boxes, step=40)
            boxes.append(n.box)
        moved += 1
    return moved > 0


def _inside_polygon(poly, x, y):
    """Whether the points (x, y) (arrays) lie inside the polygon (even-odd rule)."""
    import numpy as np

    x, y = np.asarray(x, float), np.asarray(y, float)
    inside = np.zeros(np.broadcast_shapes(x.shape, y.shape), bool)
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        if y1 == y2:
            continue
        inside ^= ((y1 > y) != (y2 > y)) & (x < x1 + (y - y1) * (x2 - x1) / (y2 - y1))
    return inside


def _orient(o, s, t):
    """Sign of the cross product (s - o) x (t - o), broadcast over the leading axes (0 within 1e-9)."""
    import numpy as np

    v = (s[..., 0] - o[..., 0]) * (t[..., 1] - o[..., 1]) - (s[..., 1] - o[..., 1]) * (t[..., 0] - o[..., 0])
    return np.where(np.abs(v) < 1e-9, 0, np.sign(v))


def _crossings(p, q, a, b):
    """Proper crossings of the segments p[i]-q with the segments a[j]-b[j]: counts per i.

    As :func:`ecmapper.pwl.layout._cross`: segments with a common end do not
    cross. `p`: (k, 2), `q`: (2,) or (k, 2), `a`, `b`: (m, 2).
    """
    import numpy as np

    if not len(a):
        return np.zeros(len(p), int)
    q = np.asarray(q, float)
    p, q, a, b = p[:, None, :], q[:, None, :] if q.ndim == 2 else q[None, None, :], a[None], b[None]
    hit = (_orient(a, b, p) * _orient(a, b, q) < 0) & (_orient(p, q, a) * _orient(p, q, b) < 0)
    common = np.zeros(hit.shape, bool)
    for u, v in ((p, a), (p, b), (q, a), (q, b)):
        common |= np.linalg.norm(u - v, axis=-1) < 1e-6
    return (hit & ~common).sum(axis=1)


def _crosses_any(p, q, a, b):
    """Whether the segments p[i]-q properly cross any of the segments a[j]-b[j] (common ends not checked).

    `p`: (k, 2), `q`: (2,) or (k, 2), `a`, `b`: (m, 2); returns (k,) bool.
    """
    import numpy as np

    if not len(a):
        return np.zeros(len(p), bool)
    q = np.asarray(q, float)
    p, q, a, b = p[:, None, :], q[:, None, :] if q.ndim == 2 else q[None, None, :], a[None], b[None]
    return ((_orient(a, b, p) * _orient(a, b, q) < 0) & (_orient(p, q, a) * _orient(p, q, b) < 0)).any(1)


def _through(p, q, boxes, pad=4.0):
    """Whether the segments p-q pass through the boxes (x0, y0, x1, y1), broadcast (Liang-Barsky).

    As :func:`ecmapper.pwl.layout._segment_hits_box`, with corner boxes.
    """
    import numpy as np

    p, q, boxes = np.asarray(p, float), np.asarray(q, float), np.asarray(boxes, float)
    d = q - p
    shape = np.broadcast_shapes(p.shape[:-1], q.shape[:-1], boxes.shape[:-1])
    t0, t1 = np.zeros(shape), np.ones(shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        for k in (0, 1):
            lo, hi, s = boxes[..., k] - pad, boxes[..., k + 2] + pad, p[..., k]
            parallel = np.abs(d[..., k]) < 1e-12
            a, b = (lo - s) / d[..., k], (hi - s) / d[..., k]
            t0 = np.where(parallel, np.where((s < lo) | (s > hi), 2.0, t0), np.maximum(t0, np.minimum(a, b)))
            t1 = np.where(parallel, t1, np.minimum(t1, np.maximum(a, b)))
    return t1 - t0 > 1e-6


def _ideal_point(v, done, pos, partners, placed, centre):
    """Where node `v` goes (:func:`place_added`): the centroid of its placed partners `done`, half a step beyond one.

    With one placed partner, the flow continues away from the partner's
    other placed neighbours (else away from `centre`, else downwards);
    without, `v` stays where it is.
    """
    import numpy as np

    from ..pwl.layout import STEP

    if len(done) >= 2:
        return np.mean([pos[q] for q in done], axis=0)
    if not done:
        return np.asarray(pos[v], float)
    q = done[0]
    away = np.sum([np.subtract(pos[q], pos[n]) for n in partners[q] & placed if n != v] or [[0.0, 0.0]], axis=0)
    if np.hypot(*away) < 1e-6:
        away = np.subtract(pos[q], centre)
    if np.hypot(*away) < 1e-6:
        away = np.array([0.0, 1.0])
    return np.asarray(pos[q]) + STEP / 2 * away / np.hypot(*away)


def place_added(s: DraftMap, outside: bool = True, weight: float | None = None, step: float = 45.0,
                rings: int = 20, margin: float = 12.0, passes: int = 2) -> bool:
    """Place the added nodes again, one by one, next to their metabolites; the nodes of the BRENDA map never move.

    The BRENDA drawing is the prior; the added enzyme nodes (with their
    cofactors) and the metabolite nodes not on the map are placed into it
    greedily: next is the node with the most main-metabolite partners already
    placed (nodes of the map count as placed; enzyme nodes before metabolite
    nodes, then more partners, then the ID). Its ideal point is the centroid
    of its placed partners or, with one, half a reaction step (``STEP / 2``)
    beyond it, away from the partner's other placed neighbours (the flow
    continues); without placed partner its current position. Around it,
    rings of points as in :func:`..core.place_free` are tried (`rings` of
    16 directions, `step` apart); the box of the node (an enzyme node with
    its cofactors and room for its reactions, as in :func:`_groups`) must
    keep `margin` from the boxes placed so far. The point with the smallest
    ``(crossings + through + inside) * weight + distance to the ideal point``
    is taken (as :func:`ecmapper.pwl.layout.place_side_metabolites`):
    crossings of its main-metabolite links with those drawn, links through
    other nodes (its own through them and theirs through its box), and with
    `outside`, 1 for a point inside the longest cycle of the BRENDA drawing
    (:func:`main_cycle`) for enzyme nodes joined to at most two nodes of the
    map and for metabolite nodes whose placed partners are all outside it
    (the rule of :func:`outside_cycle`). `weight` defaults to 1.5 ``STEP``.
    The greedy order sees only the links to nodes placed before; in
    `passes` rounds of refinement, each added node (by ID) is placed again
    the same way with all its links, its current point among the candidates.
    Deterministic. False if no node moved.
    """
    import math

    import numpy as np

    from ..pwl.layout import STEP
    from .core import place_free

    weight = 1.5 * STEP if weight is None else weight
    box, owner = _groups(s)
    movable = {i for i, n in s.nodes.items() if n.kind == "enzyme" and n.status == "added"}
    movable |= {i for i, n in s.nodes.items() if n.kind == "main" and n.status not in MAP_STATUS}
    if not movable:
        return False
    cofactors: dict[str, list[str]] = {}
    for c, enzyme in owner.items():
        cofactors.setdefault(enzyme, []).append(c)
    partners: dict[str, set[str]] = {i: set() for i in box}
    for e in s.edges:
        if e.kind == "main" and e.enzyme in box and e.met in box:
            partners[e.enzyme].add(e.met)
            partners[e.met].add(e.enzyme)
    placed = {i for i in box if i not in movable}
    pos = {i: (s.nodes[i].x, s.nodes[i].y) for i in box}
    before = {i: (n.x, n.y) for i, n in s.nodes.items()}
    # unit boxes relative to the node centre: for overlaps (with cofactors and room), for links through nodes
    unit = {i: (b[0] - pos[i][0], b[1] - pos[i][1], b[2] - pos[i][0], b[3] - pos[i][1]) for i, b in box.items()}
    own = {i: (-s.nodes[i].w / 2, -s.nodes[i].h / 2, s.nodes[i].w / 2, s.nodes[i].h / 2) for i in box}

    def shifted(offsets, x, y):
        return np.stack([x + offsets[0], y + offsets[1], x + offsets[2], y + offsets[3]], axis=-1)

    obstacles = {i: box[i] for i in placed}
    links = {tuple(sorted((a, b))) for a in placed for b in partners[a] if b in placed}
    poly, map_partners = None, {}
    if outside:
        try:
            found = main_cycle(s)
        except ImportError:  # networkx
            found = None
        if found:
            poly = [pos[m] for m in found[0]]
            map_partners = {i: sum(s.nodes[q].status in MAP_STATUS for q in partners[i]) for i in movable}
    centre = np.mean([pos[i] for i in (placed or box)], axis=0)
    angles = 2 * math.pi * np.arange(16) / 16
    ring = np.concatenate([[[0.0, 0.0]]] + [np.stack([r * step * np.cos(angles) * 1.6, r * step * np.sin(angles)], 1)
                                           for r in range(1, rings + 1)])

    def place(v, keep=False):
        """Place `v` (not in `placed`) at the best free point near its ideal one (or, with `keep`, where it is)."""
        done = sorted(partners[v] & placed)
        ideal = _ideal_point(v, done, pos, partners, placed, centre)
        points = ideal + ring
        if keep:
            points = np.vstack([points, [pos[v]]])
        free = np.ones(len(points), bool)
        if obstacles:
            b, o = shifted(unit[v], points[:, :1], points[:, 1:2]), np.array(list(obstacles.values()))[None]
            free = ~np.any((b[..., 0] < o[..., 2] + margin) & (o[..., 0] - margin < b[..., 2])
                           & (b[..., 1] < o[..., 3] + margin) & (o[..., 1] - margin < b[..., 3]), axis=1)
        if free.any():
            points = points[free]
            penalty = np.zeros(len(points))
            if links:
                a = np.array([pos[x] for x, _ in links])
                b = np.array([pos[y] for _, y in links])
                # links drawn through the box of the node, (points, links)
                penalty += _through(a[None], b[None], shifted(own[v], points[:, :1], points[:, 1:2])).sum(1)
            for q in done:
                if links:
                    penalty += _crossings(points, pos[q], a, b)
                others = [i for i in sorted(placed) if i != q]
                if others:  # its links through other nodes, (points, nodes)
                    boxes = np.array([shifted(own[i], *pos[i]) for i in others])
                    penalty += _through(points[:, None], np.asarray(pos[q])[None, None], boxes[None]).sum(1)
            if poly is not None and (map_partners[v] <= 2 if s.nodes[v].kind == "enzyme"
                                     else done and not _inside_polygon(poly, *np.array([pos[q] for q in done]).T).any()):
                penalty += _inside_polygon(poly, points[:, 0], points[:, 1])
            x, y = points[np.argmin(penalty * weight + np.hypot(*(points - ideal).T))]
        else:  # no free point near: the nearest free place
            b = shifted(unit[v], *ideal)
            cx, cy = place_free((b[0] + b[2]) / 2, (b[1] + b[3]) / 2, b[2] - b[0], b[3] - b[1],
                                list(obstacles.values()), margin=margin)
            x, y = cx - (b[0] + b[2]) / 2 + ideal[0], cy - (b[1] + b[3]) / 2 + ideal[1]
        x, y = float(x), float(y)
        dx, dy = x - pos[v][0], y - pos[v][1]
        for i in [v, *cofactors.get(v, [])]:
            s.nodes[i].x += dx
            s.nodes[i].y += dy
        pos[v] = (x, y)
        obstacles[v] = tuple(shifted(unit[v], x, y))
        links.update(tuple(sorted((v, q))) for q in done)
        placed.add(v)

    while unplaced := movable - placed:
        place(min(unplaced, key=lambda i: (-len(partners[i] & placed), s.nodes[i].kind != "enzyme",
                                           -len(partners[i]), i)))
    # refinement: each node again, now with all its links (the greedy order sees only those to placed nodes)
    for _ in range(passes):
        for v in sorted(movable):
            placed.discard(v)
            del obstacles[v]
            links.difference_update(tuple(sorted((v, q))) for q in partners[v])
            place(v, keep=True)
    return any((n.x, n.y) != before[i] for i, n in s.nodes.items())


def place_conservative(s: DraftMap, outside: bool = True, weight: float | None = None, step: float = 45.0,
                       margin: float = 12.0, reach: float | None = None, tries: int = 30) -> set[str]:
    """Add the added reactions to the BRENDA drawing only where they cross nothing; leave out the others.

    As in :func:`place_added`, the nodes of the BRENDA map never move and
    the added enzyme nodes are placed one by one (the most placed partners
    first), each with its metabolite nodes that are not placed yet and its
    cofactors. Here a node only goes where none of its links crosses a link
    drawn so far (main metabolites and cofactors) or passes through a node,
    no drawn link passes through it, and its box (with room for its
    reactions) keeps `margin` from the others. The points tried lie on a
    grid (`step` apart) at most `reach` (default 2 ``STEP``) from the ideal
    point of the node, also outside the drawing, so it may grow (no long
    links across the drawing instead). An enzyme node goes to the one
    nearest to its ideal point (:func:`place_added`; with `outside`, a point
    inside the longest cycle of the BRENDA drawing costs `weight`, default
    1.5 ``STEP``), its new metabolite nodes to the nearest to theirs and its
    cofactors to the nearest to where they were relative to it (among the
    512 and 256 nearest grid points); if one of them does not fit, the next
    of `tries` points of the enzyme node is tried. Reactions that fit
    nowhere are left out: their enzyme nodes, their cofactors and the
    metabolite nodes only they use are removed. Reactions put into the
    subsystem by hand (enzyme nodes marked `assigned`,
    :func:`..assignment.apply_assignment`) and metabolite nodes of enzyme
    nodes of the map that fit nowhere stay where they are. Deterministic.
    Returns the removed nodes.
    """
    import numpy as np

    from ..pwl.layout import STEP

    weight = 1.5 * STEP if weight is None else weight
    reach = 2 * STEP if reach is None else reach
    enzymes = {i for i, n in s.nodes.items() if n.kind == "enzyme" and n.status == "added"}
    mets = {i for i, n in s.nodes.items() if n.kind == "main" and n.status not in MAP_STATUS}
    if not enzymes and not mets:
        return set()
    cofactors: dict[str, list[str]] = {}
    for c, enzyme in _owners(s).items():
        if enzyme in enzymes:
            cofactors.setdefault(enzyme, []).append(c)
    unit = enzymes | mets | {c for cs in cofactors.values() for c in cs}
    partners: dict[str, set[str]] = {i: set() for i, n in s.nodes.items() if n.kind != "cofactor"}
    adj: dict[str, set[str]] = {i: set() for i in s.nodes}  # all links
    for e in s.edges:
        adj[e.enzyme].add(e.met)
        adj[e.met].add(e.enzyme)
        if e.kind == "main":
            partners[e.enzyme].add(e.met)
            partners[e.met].add(e.enzyme)
    room = lane_room(s)
    # boxes relative to the node centre: the node, and with room for its reactions (for overlaps)
    own = {i: np.array([-n.w / 2, -n.h / 2, n.w / 2, n.h / 2]) for i, n in s.nodes.items()}
    over = {i: o + np.array([-1, -1, 1, 1]) * np.tile(room.get(i, (0.0, 0.0)), 2) / 2 for i, o in own.items()}
    fixed = [i for i in s.nodes if i not in unit]
    pos = {i: (n.x, n.y) for i, n in s.nodes.items()}
    links = sorted({tuple(sorted((i, j))) for i in fixed for j in adj[i] if j not in unit})
    state = (pos, fixed, set(fixed), links)  # positions, placed nodes (list and set), drawn links
    poly, map_partners = None, {}
    if outside:
        try:
            found = main_cycle(s)
        except ImportError:  # networkx
            found = None
        if found:
            poly = [pos[m] for m in found[0]]
            map_partners = {i: sum(s.nodes[q].status in MAP_STATUS for q in partners[i]) for i in enzymes | mets}
    centre = np.mean([pos[i] for i in fixed if i in partners] or [pos[i] for i in partners], axis=0)

    def arrays(st):
        p, placed, _, drawn = st
        boxes = np.array([own[i] + np.tile(p[i], 2) for i in placed]).reshape(-1, 4)
        obstacles = np.array([over[i] + np.tile(p[i], 2) for i in placed]).reshape(-1, 4)
        a = np.array([p[i] for i, _ in drawn], float).reshape(-1, 2)
        b = np.array([p[j] for _, j in drawn], float).reshape(-1, 2)
        incident: dict[str, list[int]] = {}
        for k, (i, j) in enumerate(drawn):
            incident.setdefault(i, []).append(k)
            incident.setdefault(j, []).append(k)
        return placed, boxes, obstacles, a, b, incident

    def valid(points, i, st, arr):
        """Which points (n, 2) take node `i` and its links to the placed nodes: (n,) bool.

        The tests run one after the other on the points that passed the earlier ones, cheap ones first.
        """
        placed, boxes, obstacles, a, b, incident = arr
        segs = sorted(j for j in adj[i] if j in st[2])

        def free(idx):  # its box keeps `margin` from the placed ones
            bx = over[i] + np.tile(points[idx], 2)
            return ~np.any((bx[:, None, 0] < obstacles[None, :, 2] + margin)
                           & (obstacles[None, :, 0] - margin < bx[:, None, 2])
                           & (bx[:, None, 1] < obstacles[None, :, 3] + margin)
                           & (obstacles[None, :, 1] - margin < bx[:, None, 3]), axis=1)

        def crossing(j):  # no crossing with the drawn links (but those at its placed end)
            rest = np.ones(len(a), bool)
            rest[incident.get(j, [])] = False
            return lambda idx: ~_crosses_any(points[idx], st[0][j], a[rest], b[rest])

        def through(j):  # not through other placed nodes
            keep = np.array([k != j for k in placed], bool)
            return lambda idx: ~_through(points[idx][:, None], np.asarray(st[0][j], float)[None, None],
                                         boxes[None, keep]).any(1)

        def drawn(idx):  # no drawn link through it
            return ~_through(a[None], b[None], (own[i] + np.tile(points[idx], 2))[:, None]).any(1)

        tests = ([free] if len(obstacles) else []) + [crossing(j) for j in segs] + [through(j) for j in segs]
        tests += [drawn] if len(a) else []
        idx = np.arange(len(points))
        for f in tests:
            if not len(idx):
                break
            idx = idx[f(idx)]
        ok = np.zeros(len(points), bool)
        ok[idx] = True
        return ok

    def search(ideal, i, st, penalise, limit=None):
        """The free points for node `i` without crossings, best first, among the `limit` best of the grid."""
        arr = arrays(st)
        offsets = step * np.arange(-np.floor(reach / step), np.floor(reach / step) + 1)
        grid = np.stack(np.meshgrid(ideal[0] + offsets, ideal[1] + offsets), -1).reshape(-1, 2)
        cost = np.hypot(*(grid - ideal).T)
        grid, cost = grid[cost <= reach], cost[cost <= reach]
        if penalise:
            cost = cost + weight * _inside_polygon(poly, grid[:, 0], grid[:, 1])
        order = np.argsort(cost, kind="stable")[:limit]
        start, size = 0, 64
        while start < len(order):  # in growing chunks, best first
            chunk = order[start:start + size]
            yield from grid[chunk[valid(grid[chunk], i, st, arr)]]
            start, size = start + size, min(4 * size, 4096)

    def commit(st, nodes):
        p, placed, pset, drawn = st
        for i, xy in nodes.items():
            p[i] = (float(xy[0]), float(xy[1]))
            placed.append(i)
            pset.add(i)
        drawn.extend(sorted({tuple(sorted((i, j))) for i in nodes for j in adj[i] if j in pset} - set(drawn)))

    def inside_rule(v, st):
        """Whether a point inside the cycle costs `weight` for node `v` (the rule of :func:`place_added`)."""
        if poly is None:
            return False
        if v in enzymes:
            return map_partners[v] <= 2
        done = [st[0][q] for q in sorted(partners[v] & st[2])]
        return bool(done) and not _inside_polygon(poly, *np.array(done).T).any()

    def place_next_to(i, st, ideal=None, limit=512):
        """Place node `i` at the free point without crossings nearest to `ideal` (default as in place_added)."""
        if ideal is None:
            done = sorted(partners[i] & st[2])
            ideal = _ideal_point(i, done, st[0], partners, st[2], centre)
        found = next(search(ideal, i, st, i in partners and inside_rule(i, st), limit=limit), None)
        if found is not None:
            commit(st, {i: found})
        return found is not None

    def place_enzyme(e, st):
        """The state with `e`, its metabolite nodes not placed yet and its cofactors placed, or None."""
        done = sorted(partners[e] & st[2])
        ideal = _ideal_point(e, done, st[0], partners, st[2], centre)
        offsets = {c: np.subtract(st[0][c], st[0][e]) for c in cofactors.get(e, [])}
        for n, point in enumerate(search(ideal, e, st, inside_rule(e, st))):
            if n >= tries:
                break
            trial = (dict(st[0]), list(st[1]), set(st[2]), list(st[3]))
            commit(trial, {e: point})
            # the metabolites next to it, the cofactors near where they were relative to it
            if (all(place_next_to(m, trial) for m in sorted(partners[e] & mets - trial[2]))
                    and all(place_next_to(c, trial, point + d, limit=256) for c, d in offsets.items())):
                return trial
        return None

    rejected = set()
    while True:
        pset = state[2]
        todo = [i for i in enzymes if i not in pset and i not in rejected]
        todo += [m for m in mets if m not in pset and partners[m] & pset]  # metabolites of enzyme nodes of the map
        if not todo:
            break
        v = min(todo, key=lambda i: (-len(partners[i] & pset), i not in enzymes, -len(partners[i]), i))
        if v in enzymes:
            placed = place_enzyme(v, state)
            if placed is not None:
                state = placed
            elif getattr(s.nodes[v], "assigned", False):  # put here by hand (assignment table): stays, crossing
                commit(state, {i: state[0][i] for i in [v, *sorted(partners[v] & mets - pset), *cofactors.get(v, [])]})
            else:
                rejected.add(v)
        elif not place_next_to(v, state):
            commit(state, {v: state[0][v]})  # part of a reaction of the map: stays
    pos, _, pset, _ = state
    for i in unit & pset:
        s.nodes[i].x, s.nodes[i].y = pos[i]
    removed = rejected | (mets - pset)
    removed |= {c for e in rejected for c in cofactors.get(e, []) if adj[c] <= removed}
    for i in removed:
        del s.nodes[i]
    s.edges = [e for e in s.edges if e.enzyme in s.nodes and e.met in s.nodes]
    s.links = [l for l in s.links if l.a in s.nodes and l.b in s.nodes]
    return removed


def place_template(s: DraftMap) -> set[str]:
    """Keep only the added reactions that go where a node of the BRENDA map was; leave out the others.

    Step 3 records for an added enzyme node the place of a node of the
    BRENDA map with its reaction ID or EC number that is not drawn any more
    (`template`, :meth:`.core.Completion.template_spot`; e.g. the node of the
    reaction with another cofactor pair, whose match was not accepted), and
    for its new metabolite nodes and cofactors the places of the template
    node's metabolites (:meth:`.core.Completion.template_met`). An added
    enzyme node is kept if it has such a place, each of its new metabolite
    nodes has a free place of the template, too, and its links from there
    (to the metabolites drawn already and to its cofactors, at their places
    of the template or with their offsets) cross no link drawn so far and
    pass through no node, nor a drawn link through its nodes (a reaction
    of the template drawn at a node far from its metabolites would cross
    the drawing). It then goes exactly there. Enzyme nodes are taken by ID;
    of those at the same place, the first. The others are left out as in
    :func:`place_conservative` (their enzyme nodes, cofactors and the
    metabolite nodes only they use); reactions put into the subsystem by
    hand (`assigned`) stay where they are. Overlaps (other label widths) are
    left to step 4. Deterministic. Returns the removed nodes.
    """
    import numpy as np

    enzymes = sorted(i for i, n in s.nodes.items() if n.kind == "enzyme" and n.status == "added")
    if not enzymes:
        return set()
    cofactors: dict[str, list[str]] = {}
    for c, enzyme in _owners(s).items():
        cofactors.setdefault(enzyme, []).append(c)
    new = {i for i, n in s.nodes.items() if n.kind == "main" and n.status not in MAP_STATUS}
    unit = set(enzymes) | new | {c for e in enzymes for c in cofactors.get(e, [])}
    adj: dict[str, set[str]] = {i: set() for i in s.nodes}
    for e in s.edges:
        adj[e.enzyme].add(e.met)
        adj[e.met].add(e.enzyme)
    pos = {i: (n.x, n.y) for i, n in s.nodes.items()}
    placed = {i for i in s.nodes if i not in unit}
    drawn = sorted({tuple(sorted((i, j))) for i in placed for j in adj[i] if j in placed})

    def free_spot(i):  # the template place of node i, if that node of the map is not drawn here
        spot = getattr(s.nodes[i], "template", None)
        return spot if spot is not None and spot[0] not in s.nodes else None

    def box(i):
        (x, y), n = pos[i], s.nodes[i]
        return x - n.w / 2, y - n.h / 2, x + n.w / 2, y + n.h / 2

    def clear(nodes):
        """Whether `nodes` at their new places cross nothing with their links to the placed nodes and each other."""
        a = np.array([pos[i] for i, _ in drawn], float).reshape(-1, 2)
        b = np.array([pos[j] for _, j in drawn], float).reshape(-1, 2)
        links = sorted({tuple(sorted((i, j))) for i in nodes for j in adj[i] if j in placed or j in nodes})
        others = [k for k in placed | set(nodes)]
        boxes = np.array([box(k) for k in others], float).reshape(-1, 4)
        for i, j in links:
            p, q = np.array([pos[i]], float), np.array(pos[j], float)
            rest = np.array([i not in d and j not in d for d in drawn], bool)
            if _crosses_any(p, q, a[rest], b[rest])[0]:
                return False
            ends = np.array([k not in (i, j) for k in others], bool)
            if ends.any() and _through(p[0], q, boxes[ends]).any():
                return False
        return not (len(a) and any(_through(a, b, np.array(box(i))).any() for i in nodes))

    kept, taken = [], set()
    for e in enzymes:
        spot = free_spot(e)
        if spot is None or spot[0] in taken or any(m in new and m not in placed and free_spot(m) is None
                                                   for m in adj[e]):
            continue
        before = {i: pos[i] for i in [e, *cofactors.get(e, []), *(m for m in adj[e] if m in new)]}
        dx, dy = spot[1] - pos[e][0], spot[2] - pos[e][1]
        nodes = {e: (spot[1], spot[2])}
        for m in sorted(adj[e]):
            if m in new and m not in placed:  # a metabolite node of several of them: the first place
                nodes[m] = free_spot(m)[1:3]
        for c in cofactors.get(e, []):
            cof = getattr(s.nodes[c], "template", None)
            nodes[c] = (cof[1], cof[2]) if cof is not None else (pos[c][0] + dx, pos[c][1] + dy)
        pos.update(nodes)
        if not clear(nodes):
            pos.update(before)
            continue
        kept.append(e)
        taken.add(spot[0])
        drawn.extend(sorted({tuple(sorted((i, j))) for i in nodes for j in adj[i] if j in placed or j in nodes}
                            - set(drawn)))
        placed.update(nodes)
    for i in placed & unit:
        s.nodes[i].x, s.nodes[i].y = pos[i]
    stay = set(kept) | {e for e in enzymes if getattr(s.nodes[e], "assigned", False)}
    removed = set(enzymes) - stay
    removed |= {c for e in removed for c in cofactors.get(e, [])}
    removed |= {m for m in new if adj[m] <= removed}
    for i in removed:
        del s.nodes[i]
    s.edges = [e for e in s.edges if e.enzyme in s.nodes and e.met in s.nodes]
    s.links = [l for l in s.links if l.a in s.nodes and l.b in s.nodes]
    return removed


def circle_layout(s: DraftMap, gap: float = 30.0, k: float = 150.0, seed: int = 1) -> bool:
    """The longest cycle of main metabolites on a circle (see :func:`main_cycle`), the rest around it.

    The metabolites of the cycle and one enzyme node between each two of them
    (the one nearest to the middle between them) go on a circle, in the
    order and orientation closest to the current drawing; their boxes (with
    cofactors and room for parallel reactions) `gap` apart. Further enzyme
    nodes between the same two metabolites (parallel reactions) go on the
    same radius, outside. The other nodes follow in a spring layout
    (networkx); those joined to two or more nodes of the circle (chords)
    start inside it. The cofactors of the enzyme nodes on the circle turn
    with their reaction; condensing keeps the circle in shape (`rigid`).
    False (nothing changed) if there is no such cycle.
    """
    import math

    import networkx as nx

    found = main_cycle(s)
    if not found:
        return False
    cycle, g = found
    n = len(cycle)
    box, owner = _groups(s)
    pos = {i: (sum(b[0::2]) / 2, sum(b[1::2]) / 2) for i, b in box.items()}  # box centres
    edges = []  # per step of the cycle: [enzyme on the circle, parallel ones]
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        mid = ((pos[a][0] + pos[b][0]) / 2, (pos[a][1] + pos[b][1]) / 2)
        edges.append(sorted(g[a][b]["enzymes"], key=lambda e: (math.dist(pos[e], mid), e)))
    ring = [x for m, enzymes in zip(cycle, edges) for x in (m, enzymes[0])]
    cx = sum(pos[m][0] for m in cycle) / n
    cy = sum(pos[m][1] for m in cycle) / n
    current = [math.atan2(pos[x][1] - cy, pos[x][0] - cx) for x in ring]

    def place(box):
        """Angles of the ring (order and orientation closest to the current drawing) and the radius."""
        def extents(angles):  # along the circle
            return [abs((box[x][2] - box[x][0]) * math.sin(a)) + abs((box[x][3] - box[x][1]) * math.cos(a))
                    for x, a in zip(ring, angles)]

        ext = extents(current)
        for _ in range(2):
            total = sum(ext) + gap * len(ext)
            frac, run = [], 0.0
            for e in ext:
                frac.append((run + e / 2) / total)
                run += e + gap
            best = None
            for direction in (1, -1):
                offsets = [phi - direction * 2 * math.pi * f for phi, f in zip(current, frac)]
                start = math.atan2(sum(math.sin(o) for o in offsets), sum(math.cos(o) for o in offsets))
                residual = sum(1 - math.cos(o - start) for o in offsets)
                if best is None or residual < best[0]:
                    best = (residual, [start + direction * 2 * math.pi * f for f in frac])
            angles = best[1]
            ext = extents(angles)
        return angles, max(total / (2 * math.pi), 100.0)

    # the cofactors of the enzyme nodes on the circle turn with their reaction (from the one metabolite to the next)
    angles, radius = place(box)
    for i, enzymes in enumerate(edges):
        a, b = cycle[i], cycle[(i + 1) % n]
        old = math.atan2(pos[b][1] - pos[a][1], pos[b][0] - pos[a][0])
        now = math.atan2(math.sin(angles[2 * i + 2 - 2 * n * (i == n - 1)]) - math.sin(angles[2 * i]),
                         math.cos(angles[2 * i + 2 - 2 * n * (i == n - 1)]) - math.cos(angles[2 * i]))
        cos, sin = math.cos(now - old), math.sin(now - old)
        for e in enzymes:
            for c in (c for c, enzyme in owner.items() if enzyme == e):
                dx, dy = s.nodes[c].x - s.nodes[e].x, s.nodes[c].y - s.nodes[e].y
                s.nodes[c].x, s.nodes[c].y = s.nodes[e].x + dx * cos - dy * sin, s.nodes[e].y + dx * sin + dy * cos
    box, owner = _groups(s)
    angles, radius = place(box)
    new = {x: (cx + radius * math.cos(a), cy + radius * math.sin(a)) for x, a in zip(ring, angles)}

    def radial(x, a):  # extent of a box across the circle at angle a
        return abs((box[x][2] - box[x][0]) * math.cos(a)) + abs((box[x][3] - box[x][1]) * math.sin(a))

    for i, enzymes in enumerate(edges):
        # parallel enzyme nodes on the same radius, outside
        a = angles[2 * i + 1]
        r = radius + radial(enzymes[0], a) / 2
        for e in enzymes[1:]:
            new[e] = (cx + (r + gap + radial(e, a) / 2) * math.cos(a), cy + (r + gap + radial(e, a) / 2) * math.sin(a))
            r += gap + radial(e, a)
    # the other nodes: spring layout with the circle fixed; nodes joined to two or more nodes of the circle
    # (chords, e.g. the glyoxylate shunt) start inside it, the others outside
    graph = nx.Graph()
    graph.add_nodes_from(box)
    graph.add_edges_from(_directed(s))
    start_pos = dict(new)
    for x in graph:
        if x in new:
            continue
        on_ring = [y for y in graph[x] if y in new]
        if len(on_ring) >= 2:
            px, py = (sum(new[y][0] for y in on_ring) / len(on_ring), sum(new[y][1] for y in on_ring) / len(on_ring))
            start_pos[x] = (cx + 0.6 * (px - cx), cy + 0.6 * (py - cy))
        else:
            dx, dy = pos[x][0] - cx, pos[x][1] - cy
            d = math.hypot(dx, dy) or 1.0
            start_pos[x] = (cx + dx / d * max(d, radius + 150), cy + dy / d * max(d, radius + 150))
    layout = nx.spring_layout(graph, k=k, pos=start_pos, fixed=list(new), iterations=200, seed=seed)
    _move(s, {x: tuple(map(float, layout[x])) for x in box}, box, owner)
    # condensing keeps the circle (with its cofactors) in shape
    for sub in s.subsystems.values():
        sub.rigid = [list(new) + [c for c, enzyme in owner.items() if enzyme in new]]
    return True


WORKER = Path(__file__).with_name("graphviz_worker.py")

#: Graphviz layouts: program and graph attributes
GRAPHVIZ = {
    "fdp": dict(prog="fdp", K="2.0"),
    "dot": dict(prog="dot", nodesep="0.35", ranksep="0.55"),
    # soft pinning (ghost nodes); stress majorisation (neato's default mode) moves pinned nodes far away here
    "neato": dict(prog="neato", pin=False, soft=5.0, drift=72.0, mode="KK", len="2.0", sep="+6"),
}
#: prefix of the pinned ghost nodes of soft pinning (:func:`graphviz_spec`)
GHOST = "ghost:"


def graphviz_spec(s: DraftMap, prog: str, pin: bool = True, soft: float | None = None, drift: float = 72.0,
                  **graph_attr) -> tuple[dict, dict, dict]:
    """The graph of the enzyme and main metabolite nodes for :mod:`.graphviz_worker`, and the node boxes.

    ``neato`` and ``fdp`` start from the current drawing; with `pin`, the
    nodes of the BRENDA map keep their positions. With `soft` (and without
    `pin`), they are pinned softly: each gets a pinned ghost node (a point,
    :data:`GHOST` + its ID) at its position, joined by an edge of length
    `drift` (points) and weight `soft`, so it may drift a little to make
    room but not wander across the subsystem. Graphviz reads input positions in inches unless
    ``inputscale`` is set, and ``overlap=true`` keeps overlaps (removing
    them would move pinned nodes, and separate the nodes from their ghosts;
    step 4 removes them anyway).
    """
    box, owner = _groups(s)
    nodes, edges = [], list(_directed(s))
    for i, b in box.items():
        attrs = dict(width=f"{(b[2] - b[0]) / 72:.3f}", height=f"{(b[3] - b[1]) / 72:.3f}")
        if prog != "dot":
            x, y = (b[0] + b[2]) / 2, -(b[1] + b[3]) / 2  # Graphviz: y upwards
            on_map = s.nodes[i].status in MAP_STATUS
            attrs["pos"] = f"{x:.1f},{y:.1f}" + ("!" if pin and on_map else "")
            if soft and not pin and on_map:
                # the node starts 1 point beside its ghost (neato fails on nodes at the same place)
                attrs["pos"] = f"{x + 1:.1f},{y:.1f}"
                nodes.append([GHOST + i, dict(width="0.01", height="0.01", pos=f"{x:.1f},{y:.1f}!")])
                # (below 0.1 inch, neato fails)
                edges.append([i, GHOST + i, dict(len=f"{max(drift, 7.2) / 72:.3f}", weight=f"{soft:g}")])
        nodes.append([i, attrs])
    spec = dict(prog=prog, directed=prog == "dot",
                graph={"start": "random1", "inputscale": "72", "overlap": "true", **graph_attr},
                node=dict(shape="box", fixedsize="true", label=""), nodes=nodes, edges=edges)
    return spec, box, owner


def graphviz_positions(spec: dict) -> dict[str, list[float]]:
    """Run a Graphviz layout in a process of its own (see the module doc)."""
    result = subprocess.run([sys.executable, str(WORKER)], input=json.dumps(spec), capture_output=True, text=True,
                            check=False)
    if result.returncode:
        raise RuntimeError(f"Graphviz {spec['prog']} failed: {result.stderr.strip()}")
    return json.loads(result.stdout)


def apply_graphviz(s: DraftMap, positions: dict[str, list[float]], box, owner) -> None:
    """Move the nodes to the positions of a Graphviz layout; with ghost nodes (soft pinning), back to where
    those were (Graphviz may move the whole drawing)."""
    ghosts = [i.removeprefix(GHOST) for i in positions if i.startswith(GHOST)]
    dx = sum(positions[GHOST + i][0] - (box[i][0] + box[i][2]) / 2 for i in ghosts) / max(len(ghosts), 1)
    dy = sum(-positions[GHOST + i][1] - (box[i][1] + box[i][3]) / 2 for i in ghosts) / max(len(ghosts), 1)
    _move(s, {i: (x - dx, -y - dy) for i, (x, y) in positions.items() if i in box}, box, owner)


def graphviz_layout(s: DraftMap, prog: str, pin: bool = True, soft: float | None = None, drift: float = 72.0,
                    **graph_attr) -> None:
    """A Graphviz layout of the enzyme and main metabolite nodes (pygraphviz, in a process of its own)."""
    spec, box, owner = graphviz_spec(s, prog, pin, soft, drift, **graph_attr)
    apply_graphviz(s, graphviz_positions(spec), box, owner)


def untangle_layout(s: DraftMap, pin: bool = True) -> bool:
    """Resolve crossings of main-metabolite links locally (:func:`ecmapper.pwl.layout.untangle`).

    Graphviz ``neato`` moves only the nodes at the ends of crossing links (the
    boxes of :func:`_groups`), the rest keeps its place; with `pin`, the nodes
    of the BRENDA map never move. Cofactors move with their enzyme node.
    False if nothing changed.
    """
    from ..pwl.layout import untangle

    box, owner = _groups(s)
    pos = {i: ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2) for i, b in box.items()}
    sizes = {i: (b[2] - b[0], b[3] - b[1]) for i, b in box.items()}
    pinned = {i for i in box if s.nodes[i].status in MAP_STATUS} if pin else set()
    new = untangle(s, pos, sizes=sizes, pinned=pinned)
    # pinned nodes stay exactly where they are (Graphviz rounds the input positions)
    new = {i: p for i, p in new.items() if i not in pinned and p != pos[i]}
    if not new:
        return False
    _move(s, new, box, owner)
    return True


def autolayout_layout(s: DraftMap, scale: float = 3.0, max_reactions: int = 150):
    """sbml-autolayout (package ``sbml_autolayout``) of the enzyme and main metabolite nodes; redraws all of them.

    Each enzyme node is one reaction of the main metabolite nodes at its
    main edges (directions as in :func:`_directed`), each main metabolite
    node a metabolite of its own, with its name on the BRENDA map and one of
    its model metabolites as ID (known pathways and cycles are found by
    these). Only the node centres are taken, `scale`
    times the coordinates of the package (110 between metabolites there);
    cofactors move with their enzyme node, and condensing keeps the circles
    of known cycles in shape (`rigid`). Nodes without main edges go in a
    row below the drawing. False (no layout) without main edges or with
    more than `max_reactions` enzyme nodes (the run time grows fast).
    """
    from sbml_autolayout.layout import Params, compute_layout
    from sbml_autolayout.model import pathway_from_dict

    box, owner = _groups(s)
    reactions = {}
    for a, b in _directed(s):
        enzyme, met, side = (a, b, "products") if s.nodes[a].kind == "enzyme" else (b, a, "substrates")
        r = reactions.setdefault(enzyme, dict(id=enzyme, name=s.nodes[enzyme].label, substrates={}, products={},
                                              reversible=False, main=[], main_explicit=True))
        if met not in r["main"]:
            r[side][met] = 1.0
            r["main"].append(met)
    if not reactions or len(reactions) > max_reactions:
        return False
    for e in s.edges:
        if e.kind == "main" and e.arrow_enzyme == e.arrow_met:
            reactions[e.enzyme]["reversible"] = True
    # the name of the BRENDA map (the title; the label may be an abbreviation) and one model metabolite, the one
    # with the shortest ID (a BiGG ID rather than one of ModelSEED)
    metabolites = {i: dict(name=(n.title or getattr(n, "full_label", n.label)).split(";")[0],
                           base=min(n.model_mets, key=lambda m: (len(m), m), default=i))
                   for i, n in s.nodes.items() if n.kind == "main"}
    result = compute_layout(pathway_from_dict(dict(id="subsystem", metabolites=metabolites,
                                                   reactions=list(reactions.values()))),
                            Params(pseudo="draw", reaction_labels=False))
    centres = {i: (x * scale, y * scale) for i, (x, y) in {**result.met_pos, **result.rxn_pos}.items() if i in box}
    # nodes that are not in the layout: in a row below it
    x, y = min(c[0] for c in centres.values()), max(c[1] for c in centres.values())
    for i in sorted(set(box) - set(centres)):
        w, h = box[i][2] - box[i][0], box[i][3] - box[i][1]
        centres[i] = (x + w / 2, y + 110 * scale / 2 + h / 2)
        x += w + 20
    _move(s, centres, box, owner)
    # condensing keeps the circles of the known cycles (with their cofactors) in shape
    rings = [[i for _, i in c.slots if i in box] for c in result.cycles]
    for sub in s.subsystems.values():
        sub.rigid = [ring + [c for c, enzyme in owner.items() if enzyme in ring] for ring in rings]


LAYOUTS = {"spring": spring_layout, "circle": circle_layout, AUTOLAYOUT: autolayout_layout,
           **{m: (lambda s, a=a: graphviz_layout(s, **a)) for m, a in GRAPHVIZ.items()}}


def available_methods(methods=METHODS) -> list[str]:
    """The methods whose packages are installed."""
    result = []
    for m in methods:
        package = {"spring": "networkx", "circle": "networkx", "fdp": "pygraphviz", "dot": "pygraphviz",
                   "neato": "pygraphviz", AUTOLAYOUT: "sbml_autolayout"}.get(m)
        if package:
            try:
                __import__(package)
            except ImportError:
                continue
        result.append(m)
    return result


def _base(m: str) -> str:
    """The layout a candidate starts from (``spring+untangle`` -> ``spring``)."""
    return m.removesuffix(UNTANGLED)


def layout_key(m: str) -> tuple[int, bool]:
    """Sort key of layouts: the current drawing, the placements, the other layouts; each before its untangled form."""
    order = ["current", CONSERVATIVE, TEMPLATE, "place", "outside", *METHODS, "circle", AUTOLAYOUT]
    return order.index(_base(m)) if _base(m) in order else len(order), m != _base(m)


def layout_methods(m: str) -> tuple[list[str], bool, bool]:
    """The methods, whether `circle` and whether `untangle`, with which :func:`relayout_subsystems` computes
    layout `m` as :func:`compare_layouts` does (``spring+untangle``: ``place`` and ``spring``, untangled)."""
    base, untangled = _base(m), m.endswith(UNTANGLED)
    if base in (CONSERVATIVE, TEMPLATE):
        return [base], False, False
    if base == "circle":
        return ["place"], True, False
    if base == "outside":
        return [], False, untangled
    if base == AUTOLAYOUT:  # redraws the whole subsystem: no placement before it
        return [base], False, untangled
    return ["place", *([base] if base not in ("current", "place") else [])], False, untangled


def relayout_subsystems(dm: DraftMap, methods=METHODS, margin: float = 0.25, circle=(), untangle: bool = True,
                        min_gain: int = 3, outside: bool = True, found: dict | None = None,
                        by_hand: dict[str, str] | None = None, only: set[str] | None = None,
                        uniform: bool = False, **condense_args) -> list[dict]:
    """Step 3b: lay out each subsystem with added reactions with `methods`, keep the best layout (see module doc).

    The layouts are compared after condensing the subsystem as in step 4
    (`condense_args`). A redrawing layout (``dot``, ``neato``) is taken only
    with at least `margin` fewer and at least `min_gain` fewer crossings of
    main-metabolite links than the best keeping layout. With `untangle`, each
    layout (but the circle) is also tried untangled (:func:`untangle_layout`,
    as ``<layout>+untangle``; in layouts that keep the BRENDA drawing only
    the added nodes move). With ``place`` among the `methods`, the added
    nodes are placed again one by one next to their metabolites
    (:func:`place_added`, `outside`: not inside a cycle of the BRENDA
    drawing): the placed drawing is the candidate ``place``, and the other
    layouts (but ``current``) start from it. Without ``place``, with
    `outside`, added reactions at a cycle of the BRENDA drawing start outside
    it (:func:`outside_cycle`), as the candidate ``outside`` from which the
    other layouts start. The subsystems in `circle` get the circular
    layout of their cycle (:func:`circle_layout`) if they have one. With
    :data:`CONSERVATIVE` (alone), the added reactions are placed only where
    they cross nothing (:func:`place_conservative`), no other layout is
    tried and the reactions that do not fit are removed from `dm` (with
    their inter-pathway links); the rows list them (``left_out``). The same
    with :data:`TEMPLATE` (alone): only the added reactions at the place of a
    node of the BRENDA map are drawn (:func:`place_template`). With
    `uniform`, the subsystems of these two placements are only shrunk
    uniformly instead of condensed (:func:`.layout.shrink_uniform`). Changes
    the node positions in `dm` (before step 4) and the `layout` and `rigid`
    of the subsystems. Returns one row per subsystem: the crossings of each
    layout and the one kept. `found`: gets {subsystem: {layout: dict(
    crossings_main, crossings_all, area, map: the subsystem condensed,
    left_out: model reactions not drawn)}} for a comparison
    (:func:`compare_layouts`). `by_hand`: {subsystem: layout} chosen by
    hand (``--layout-of``, the GUI): these subsystems get that layout,
    computed as :func:`compare_layouts` does (:func:`layout_methods`), the
    rows say ``by_hand``. `only`: lay out only these subsystems.
    """
    by_hand = {k: m for k, m in (by_hand or {}).items() if k in dm.subsystems and (only is None or k in only)}
    if by_hand and set(by_hand) != only:  # each subsystem laid out by hand on its own, with its methods
        args = dict(margin=margin, min_gain=min_gain, outside=outside, found=found, uniform=uniform, **condense_args)
        rest = set(dm.subsystems if only is None else only) - set(by_hand)
        rows = relayout_subsystems(dm, methods, circle=circle, untangle=untangle, only=rest, **args)
        for k, m in sorted(by_hand.items()):
            ms, on_circle, untangled = layout_methods(m)
            rows += relayout_subsystems(dm, ms, circle=[k] if on_circle else (), untangle=untangled,
                                        by_hand={k: m}, only={k}, **args)
        return rows
    asked = [m for m in available_methods(methods) if m != "current"]
    conservative, templated = CONSERVATIVE in asked, TEMPLATE in asked
    alone = CONSERVATIVE if conservative else TEMPLATE if templated else None  # a placement used alone
    if alone and len(asked) > 1:
        raise ValueError(f"{alone} is not combined with other layouts: {asked}")
    untangle = untangle and not alone and available_methods(["neato"]) == ["neato"]  # needs pygraphviz
    placing = "place" in asked
    methods = [m for m in asked if m not in ("place", CONSERVATIVE, TEMPLATE)]  # the layouts started from the base
    circle = set(circle) & set(dm.subsystems) if available_methods(["circle"]) else set()  # needs networkx
    # needs networkx; no use before a layout that redraws everything
    outside = outside and available_methods(["circle"]) == ["circle"] and asked != [AUTOLAYOUT]
    keys = [k for k in sorted(dm.subsystems) if (only is None or k in only)
            and (k in circle or by_hand.get(k) is not None and _base(by_hand[k]) == AUTOLAYOUT
                 or any(n.kind == "enzyme" and n.status == "added" for n in dm.members(k)))]
    currents = {k: subsystem_map(dm, k) for k in keys}
    bases = {k: copy.deepcopy(s) for k, s in currents.items()}
    # the candidate of the base, if it differs from the current drawing
    extra = alone or ("place" if placing else "outside")
    removed = {k: place_conservative(bases[k], outside=outside) if conservative
               else place_template(bases[k]) if templated else set() for k in keys}
    moved = {k: bool(alone) or (place_added(bases[k], outside=outside) if placing
                                else outside and outside_cycle(bases[k])) for k in keys}
    # the Graphviz layouts of all subsystems, several processes at a time
    jobs = {(k, m): graphviz_spec(bases[k], **GRAPHVIZ[m]) for k in keys for m in methods if m in GRAPHVIZ}
    with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        computed = dict(zip(jobs, pool.map(graphviz_positions, [spec for spec, _, _ in jobs.values()])))
    rows = []
    for k in keys:
        base = bases[k]
        results = {}

        def evaluate(m, s):
            positions = {i: (n.x, n.y) for i, n in s.nodes.items()}
            rigid = copy.deepcopy(getattr(s.subsystems[k], "rigid", []))
            condense_subsystems(s, keep_crossings=bool(alone), uniform=uniform and bool(alone), **condense_args)
            metrics = internal_metrics(s)
            results[m] = (metrics["crossings_main"], metrics["crossings_all"], s.subsystems[k].w * s.subsystems[k].h,
                          positions, rigid, s)

        for m in ["current", *([extra] if moved[k] else []), *methods, *(["circle"] if k in circle else [])]:
            s = copy.deepcopy(currents[k] if m == "current" else base)
            if m in GRAPHVIZ:
                apply_graphviz(s, computed[k, m], *jobs[k, m][1:])
            elif m in LAYOUTS and LAYOUTS[m](s) is False:  # not "current" and the base
                continue  # no cycle
            if untangle and m != "circle":
                u = copy.deepcopy(s)
                if untangle_layout(u, pin=m in KEEPING):
                    evaluate(m + UNTANGLED, u)
            evaluate(m, s)
        rank = {m: r[:3] for m, r in results.items()}  # crossings (main, all), area
        best = alone or min((m for m in results if _base(m) in KEEPING), key=rank.get)
        redraw = [m for m in results if _base(m) in REDRAWING]
        if redraw:
            other = min(redraw, key=rank.get)
            if (results[other][0] < (1 - margin) * results[best][0]
                    and results[best][0] - results[other][0] >= min_gain):
                best = other
        if "circle" in results:  # asked for
            best = "circle"
        if k in by_hand:
            if by_hand[k] in results:
                best = by_hand[k]
            else:
                dm.log.append(f"step 3b: layout {by_hand[k]!r} of {k} (chosen by hand) not available, {best} kept")
        left_out = sorted(r for i in removed[k] for r in dm.nodes[i].reactions) if best != "current" else []
        if best != "current" and removed[k]:
            dm.remove_nodes(removed[k])
            dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]
        for i, (x, y) in results[best][3].items():
            dm.nodes[i].x, dm.nodes[i].y = x, y
        dm.subsystems[k].layout, dm.subsystems[k].rigid = best, results[best][4]
        rows.append(dict(subsystem=k, layout=best, **{m: dict(crossings_main=r[0], crossings_all=r[1], area=r[2])
                                                       for m, r in results.items()}))
        if alone:
            rows[-1]["left_out"] = left_out
        if k in by_hand:
            rows[-1]["by_hand"] = True
        if found is not None:
            gone = sorted(r for i in removed[k] for r in currents[k].nodes[i].reactions)
            found[k] = {m: dict(crossings_main=r[0], crossings_all=r[1], area=r[2], map=r[5],
                                left_out=gone if m != "current" else []) for m, r in results.items()}
    names = ["current", *([extra] if any(moved.values()) else []), *methods, *(["circle"] if circle else [])]
    names += [m + UNTANGLED for m in names if untangle and m != "circle"]
    chosen = {m: n for m in names if (n := sum(r["layout"] == m for r in rows))}
    before = sum(r["current"]["crossings_main"] for r in rows)
    after = sum(r[r["layout"]]["crossings_main"] for r in rows)
    circles = f" ({sum(r['layout'] == 'circle' for r in rows)} of them on a circle, --circle)" if circle else ""
    if conservative:
        rings = (f"; added reactions placed where they cross nothing, {sum(len(r['left_out']) for r in rows)} "
                 f"model reactions left out")
    elif templated:
        rings = (f"; added reactions placed where a node of the BRENDA map was, "
                 f"{sum(len(r['left_out']) for r in rows)} model reactions without such a place left out")
    elif placing:
        rings = f"; added nodes placed again next to their metabolites in {sum(moved.values())}"
    else:
        rings = f"; added reactions start outside the cycle in {sum(moved.values())}" if any(moved.values()) else ""
    dm.log.append(
        f"step 3b: {len(rows)} subsystems with added reactions laid out again ({', '.join(asked) or 'no layouts'}): "
        f"{chosen}; crossings of main-metabolite links in them {before} -> {after}{circles}{rings}"
    )
    return rows


def compare_layouts(s: DraftMap, untangle: bool = True, outside: bool = True, margin: float = 0.25,
                    min_gain: int = 3, uniform: bool = False, **condense_args) -> tuple[dict[str, dict], str | None]:
    """All layouts of step 3b for one subsystem, to compare them by eye.

    `s`: the subsystem before step 3b, as a draft map of its own
    (:func:`subsystem_map`). The layouts are those of
    :func:`relayout_subsystems` with all available :data:`METHODS` (with
    `untangle`, also untangled), the circular layout if the subsystem has a
    cycle, and :data:`CONSERVATIVE`; each condensed as in step 4
    (`condense_args`; with `uniform`, the placements used alone only
    shrunk uniformly). Returns {layout: dict(crossings_main, crossings_all,
    area, map: the condensed subsystem, left_out: model reactions not
    drawn)} and the layout step 3b takes from :data:`METHODS` (with `margin`
    and `min_gain`; None if the subsystem has no added reactions).
    """
    k = next(iter(s.subsystems))
    found: dict = {}
    rows = relayout_subsystems(copy.deepcopy(s), METHODS, margin=margin, untangle=untangle, min_gain=min_gain,
                               outside=outside, found=found, **condense_args)
    result = dict(found.get(k, {}))
    for methods, circle in (([CONSERVATIVE], ()), ([TEMPLATE], ()), (["place"], (k,))):
        found = {}
        relayout_subsystems(copy.deepcopy(s), methods, circle=circle, untangle=False, outside=outside, found=found,
                            uniform=uniform, **condense_args)
        result.update({m: r for m, r in found.get(k, {}).items() if m not in result})
    found = {}  # chosen by hand only, also without added reactions
    relayout_subsystems(copy.deepcopy(s), [AUTOLAYOUT], untangle=untangle, outside=False, found=found,
                        by_hand={k: AUTOLAYOUT}, only={k}, **condense_args)
    result.update({m: r for m, r in found.get(k, {}).items() if m not in result})
    if "current" not in result:  # no added reactions, no cycle
        c = copy.deepcopy(s)
        condense_subsystems(c, **condense_args)
        metrics = internal_metrics(c)
        result["current"] = dict(crossings_main=metrics["crossings_main"], crossings_all=metrics["crossings_all"],
                                 area=c.subsystems[k].w * c.subsystems[k].h, map=c, left_out=[])
    return result, rows[0]["layout"] if rows else None
