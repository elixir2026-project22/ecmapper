"""Layouts of PWL subsystems in the style of the BRENDA map (networkx and Graphviz).

A PWL document is turned into a one-subsystem :class:`~ecmapper.draft_map.core.DraftMap`
(:func:`to_draftmap`): an enzyme node per reaction, a main-metabolite node per
drawn compound (``token#n``), a cofactor node per cofactor of a reaction. Then
several layouts of the enzyme and main-metabolite nodes are computed
(:func:`layout_candidates`):

``hints``
    The drawing hints of the PWL file, replayed as a turtle walk on an
    orthogonal grid (:func:`hint_positions`): every chain, ring run and spoke
    starts with the heading of ``flow``; a link turns the heading by its turn
    letter (``R`` 90 degrees clockwise on screen, ``L`` counter-clockwise,
    ``B`` back) and takes one step of :data:`STEP` map units; star spokes are
    spread over the ``arc``; groups ``(A + B)`` are spread across the step;
    the reactions of a fan lie side by side.
``hints+neato``
    ``hints``, untangled by Graphviz ``neato`` with all nodes pinned except
    those at the ends of crossing links (:func:`untangle`).
``hints+spring``
    networkx ``spring_layout`` started from ``hints`` (few iterations), then
    links that are nearly horizontal or vertical are straightened.
``hints+spring+neato``
    ``hints+spring``, untangled like ``hints+neato``.
baselines
    ``dot`` (ranks along ``flow``), ``neato``, ``fdp``, ``sfdp`` (Graphviz),
    ``kamada_kawai`` and ``spring`` (networkx) and ``circle``
    (:func:`ecmapper.draft_map.relayout.circle_layout`, subsystems with a cycle)
    use only the graph and ``flow``.

All layouts are scaled to the BRENDA spacing; cofactors are then placed next
to their enzyme (:func:`place_cofactors`) on the side given by ``^L``/``^R``,
consumed ones upstream, produced ones downstream, like on the BRENDA map. The
drawn layout is the best ``hints`` layout by crossings, unless a baseline has
clearly fewer (:func:`choose`). :mod:`ecmapper.pwl.compare` compares the
layouts with the original drawing.

Needs numpy, networkx (and scipy for ``kamada_kawai``); the Graphviz layouts
need pygraphviz, the SVG needs the ``draft-map`` dependencies.
"""

from __future__ import annotations

import copy
import math
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from ..draft_map.core import DraftMap, Edge, Node, Subsystem, text_width
from .syntax import Chain, Document, React, Ring, Star

#: distance of two main metabolites of a reaction (map units; median of the BRENDA map: 327)
STEP = 330.0
#: distance of a cofactor node from its enzyme node (centre to centre)
COFACTOR_DISTANCE = 80.0
#: free space between nodes of a group, and between parallel reactions of a fan
GAP = 30.0
#: subsystem box: padding and title gap as in :mod:`ecmapper.draft_map.layout`
PADDING, TITLE_GAP = 35.0, 22.0
HEAD = {"down": (0, 1), "up": (0, -1), "right": (1, 0), "left": (-1, 0)}
HINT_METHODS = ("hints", "hints+neato", "hints+spring", "hints+spring+neato")
BASELINES = ("dot", "neato", "fdp", "sfdp", "kamada_kawai", "spring", "circle")
METHODS = HINT_METHODS + BASELINES


def _rot(h, turn):
    """The heading after a turn letter; screen coordinates (y down), R = 90 degrees clockwise."""
    x, y = h
    return {"": (x, y), "F": (x, y), "R": (-y, x), "L": (y, -x), "B": (-x, -y)}[turn]


def _unit(v):
    n = math.hypot(*v)
    return (v[0] / n, v[1] / n) if n > 1e-9 else (0.0, 1.0)


# ---------------------------------------------------------------------------------- PWL -> DraftMap


def to_draftmap(doc: Document, node_ids: dict[str, str] | None = None, fill: str = "#dedede") -> DraftMap:
    """A draft map of the subsystem of `doc`; positions are set by a layout later.

    `node_ids` maps drawn nodes (``token#n``) to node IDs, e.g. the BRENDA node
    IDs of ``node_index.json``; other nodes are named by their token. Enzyme
    nodes are named by the reaction ID (``r1``, ``r2``, ... if there is none).
    The reactions are kept in ``dm.pwl_reactions`` as (enzyme node, reaction).
    """
    node_ids = node_ids or {}
    dm = DraftMap()
    sid = doc.id
    dm.subsystems[sid] = Subsystem(sid, doc.title, "", fill, 0, 0, 0, 0, text_width(doc.title, 18))

    def label(token):
        c = doc.compounds.get(token)
        return c.label if c is not None and c.label else token

    def ligand(token):
        c = doc.compounds.get(token)
        ids = c.xref("brenda") if c is not None else []
        return int(ids[0]) if ids and ids[0].isdigit() else None

    def main_node(p):
        nid = node_ids.get(p.node, p.node)
        if nid not in dm.nodes:
            text = label(p.token)
            dm.nodes[nid] = Node(nid, "main", sid, 0.0, 0.0, text_width(text), 30.0, text, svg_class="nodeHM",
                                 ligand_id=ligand(p.token))
        return nid

    dm.pwl_reactions = []
    for k, r in enumerate(doc.reactions(), 1):
        eid = r.rid or f"r{k}"
        text = r.enzyme_label if r.enzymes else (r.rid or "")
        title = r.comment
        dm.nodes[eid] = Node(eid, "enzyme", sid, 0.0, 0.0, max(text_width(text), 40.0), 30.0, text, title,
                             "nodeEN", ec_label=r.enzyme_label if r.enzymes and r.enzymes != ["?"] else "")
        dm.pwl_reactions.append((eid, r))
        subs = {id(p) for p in (r.substrates or [])}
        reversible = r.reversible is True
        for p in r.participants:
            nid = main_node(p)
            if r.substrates is None or reversible:
                arrows = (True, True)  # unknown or both directions: double-headed, as on the BRENDA map
            elif id(p) in subs:
                arrows = (True, False)
            else:
                arrows = (False, True)
            dm.edges.append(Edge(eid, nid, "main", *arrows))
        for j, (role, p) in enumerate([("consumed", p) for p in r.consumed] + [("produced", p) for p in r.produced]):
            cid = f"{eid}~{j}"
            text = label(p.token)
            dm.nodes[cid] = Node(cid, "cofactor", sid, 0.0, 0.0, text_width(text), 30.0, text, svg_class="nodeCM",
                                 ligand_id=ligand(p.token))
            dm.edges.append(Edge(eid, cid, "cofactor", role == "consumed", role == "produced"))
    return dm


# ---------------------------------------------------------------------------------- secondary metabolites


def _main_components(dm: DraftMap, skip=frozenset()) -> int:
    import networkx as nx

    g = nx.Graph()
    g.add_nodes_from(i for i, n in dm.nodes.items() if n.kind != "cofactor" and i not in skip)
    g.add_edges_from((e.enzyme, e.met) for e in dm.edges if e.kind == "main" and e.met not in skip)
    return nx.number_connected_components(g)


def side_metabolites(doc: Document, dm: DraftMap, node_ids: dict[str, str] | None = None) -> set[str]:
    """Secondary main metabolites: drawn only as a partner in groups ``(A + B)`` of several reactions, never as a
    node of their own (e.g. acetyl-CoA in ``(Ac_CoA + oxaloacetate)``).

    They are laid out after the others (:func:`place_side_metabolites`), next to
    their enzymes: as group members they would tie distant reactions together
    and pull the drawing into knots. A metabolite stays in the layout if the
    drawing would fall apart without it, or if a group has no other member.
    """
    node_ids = node_ids or {}

    def nid(p):
        return node_ids.get(p.node, p.node)

    alone, groups = set(), []

    def visit(statements):
        for t in statements:
            if isinstance(t, Ring):
                visit(t.statements)
            elif isinstance(t, Star):
                alone.update(nid(p) for p in t.hub.participants)
                visit(t.spokes)
            elif isinstance(t, React):
                alone.update(nid(p) for p in t.reaction.unsided)
            elif isinstance(t, Chain):
                for node in t.nodes:
                    ids = [nid(p) for p in node.participants]
                    if len(set(ids)) == 1:
                        alone.update(ids)
                    elif ids:
                        groups.append(set(ids))

    visit(doc.statements)
    degree = defaultdict(set)
    for e in dm.edges:
        if e.kind == "main":
            degree[e.met].add(e.enzyme)
    side = {i for g in groups for i in g if i not in alone and len(degree[i]) >= 2}
    for g in sorted(groups, key=sorted):  # a group keeps a member: the one with the fewest reactions
        if g <= side:
            side.discard(min(g, key=lambda i: (len(degree[i]), i)))
    components = _main_components(dm)
    # the drawing must not fall apart without them
    side = {i for i in side if _main_components(dm, frozenset({i})) <= components}
    while side and _main_components(dm, frozenset(side)) > components:
        side.discard(min(side, key=lambda i: (len(degree[i]), i)))
    return side


def _segment_hits_box(p, q, box, pad=4.0) -> bool:
    """Whether segment p-q passes through the box (x, y, w, h centred), by Liang-Barsky clipping."""
    x, y, w, h = box
    x0, x1, y0, y1 = x - w / 2 - pad, x + w / 2 + pad, y - h / 2 - pad, y + h / 2 + pad
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for d, lo, hi, s in ((dx, x0, x1, p[0]), (dy, y0, y1, p[1])):
        if abs(d) < 1e-12:
            if s < lo or s > hi:
                return False
            continue
        a, b = (lo - s) / d, (hi - s) / d
        a, b = min(a, b), max(a, b)
        t0, t1 = max(t0, a), min(t1, b)
        if t0 > t1:
            return False
    return t1 - t0 > 1e-6


def place_side_metabolites(dm: DraftMap, pos: dict, side) -> dict:
    """Place each secondary metabolite near the enzymes it links to: at the free point of a grid around their
    centroid whose links cross the fewest main links and pass through the fewest nodes, then the closest one."""
    pos = dict(pos)
    enzymes = defaultdict(list)
    for e in dm.edges:
        if e.kind == "main" and e.met in side:
            enzymes[e.met].append(e.enzyme)
    for sid_ in sorted(side, key=lambda i: (-len(enzymes[i]), i)):
        es = [e for e in enzymes[sid_] if e in pos]
        if not es:
            continue
        centre = _centroid([pos[e] for e in es], (0.0, 0.0))
        segments = [(pos[a], pos[b]) for a, b in _segments(dm, pos)]
        n = dm.nodes[sid_]
        boxes = {i: (x, y, dm.nodes[i].w, dm.nodes[i].h) for i, (x, y) in pos.items()
                 if dm.nodes[i].kind != "cofactor"}
        best, key = None, None
        q = STEP / 4
        for i in range(-6, 7):
            for j in range(-6, 7):
                c = (centre[0] + i * q, centre[1] + j * q)
                if any(abs(c[0] - x) < (n.w + w) / 2 + GAP / 2 and abs(c[1] - y) < (n.h + h) / 2 + GAP / 2
                       for x, y, w, h in boxes.values()):
                    continue
                if min(math.dist(c, pos[e]) for e in es) < STEP / 4:
                    continue
                crossings = sum(_cross(c, pos[e], a, b) for e in es for a, b in segments)
                through = sum(_segment_hits_box(c, pos[e], box) for e in es for o, box in boxes.items() if o != e)
                k = (crossings + through) * 1.5 * STEP + sum(math.dist(c, pos[e]) for e in es) / len(es)
                if key is None or k < key:
                    best, key = c, k
        pos[sid_] = best if best is not None else (centre[0] + STEP / 3, centre[1])
    return pos


# ---------------------------------------------------------------------------------- hint layout


class _Turtle:
    """Places the main metabolites and enzymes of a document by its drawing hints."""

    def __init__(self, doc: Document, dm: DraftMap, node_ids: dict[str, str], circle_rings: bool = True):
        self.doc, self.dm, self.node_ids = doc, dm, node_ids
        self.side = getattr(dm, "side", set())  # placed later, see place_side_metabolites
        self.circle_rings = circle_rings
        self.head0 = HEAD.get(doc.flow or "down", (0, 1))
        self.pos: dict[str, tuple[float, float]] = {}
        self.links = []  # (reactions, from-side node IDs, to-side node IDs, direction, from point, to point)
        self.segments = []  # links drawn so far (point, point), to avoid crossings
        self.mids = []  # their middles, where the enzymes will go
        self.rid = {id(r): eid for eid, r in dm.pwl_reactions}

    def nid(self, p):
        return self.node_ids.get(p.node, p.node)

    # --- free places
    def _free(self, nid, point, src=None) -> bool:
        """Whether a node fits at `point` (and the enzyme of its link from `src` at the middle)."""
        n = self.dm.nodes[nid]
        boxes = [(x, y, self.dm.nodes[o].w, self.dm.nodes[o].h) for o, (x, y) in self.pos.items() if o != nid]
        boxes += [(x, y, 70.0, 30.0) for x, y in self.mids]
        tests = [(point, n.w, n.h)]
        if src is not None:
            tests.append((((src[0] + point[0]) / 2, (src[1] + point[1]) / 2), 70.0, 30.0))
        for (px, py), w, h in tests:
            for x, y, ow, oh in boxes:
                if abs(px - x) < (w + ow) / 2 + GAP / 2 and abs(py - y) < (h + oh) / 2 + GAP / 2:
                    return False
        return True

    def _crossings(self, p, q) -> int:
        return sum(_cross(p, q, a, b) for a, b in self.segments)

    def _place(self, nid, point, direction, src=None):
        """Place a node at `point` or, if that place is taken, at the free point of a half-step grid around it
        whose link from `src` crosses the fewest links drawn so far, then the nearest one."""
        if nid in self.pos:
            return
        if self._free(nid, point, src) and (src is None or not self._crossings(src, point)):
            self.pos[nid] = point
            return
        best, key = None, None
        h = STEP / 2
        for i in range(-4, 5):
            for j in range(-4, 5):
                cand = (point[0] + i * h, point[1] + j * h)
                if src is not None and math.dist(cand, src) < h:
                    continue
                if not self._free(nid, cand, src):
                    continue
                k = ((self._crossings(src, cand) if src is not None else 0) * 1.5 * STEP + math.dist(cand, point))
                if key is None or k < key:
                    best, key = cand, k
        self.pos[nid] = best if best is not None else point

    def _new_component(self):
        """A start point for a part not connected to the drawn ones: next to them, across the flow."""
        if not self.pos:
            return (0.0, 0.0)
        xs = [x for x, _ in self.pos.values()]
        ys = [y for _, y in self.pos.values()]
        if self.head0[1] != 0:  # vertical flow: to the right
            return (max(xs) + STEP, min(ys))
        return (min(xs), max(ys) + STEP)

    def _node_point(self, node, direction, at=None, src=None):
        """Point of a chain node; places its unplaced members around `at` (or the point of the placed ones)."""
        ids = [self.nid(p) for p in node.participants if self.nid(p) not in self.side]
        placed = [self.pos[i] for i in ids if i in self.pos]
        if placed:
            centre = (sum(x for x, _ in placed) / len(placed), sum(y for _, y in placed) / len(placed))
        elif at is not None:
            centre = at
        else:
            centre = self._new_component()
        free = [i for i in ids if i not in self.pos]
        if free:
            across = (-direction[1], direction[0])

            def size(i):  # extent of a node across the step, with a gap
                return abs(across[0]) * self.dm.nodes[i].w + abs(across[1]) * self.dm.nodes[i].h + GAP

            if not placed:
                # the whole group (or the single node) centred on the point, members side by side
                offset = -sum(size(i) for i in ids) / 2
                for i in ids:
                    c = offset + size(i) / 2
                    offset += size(i)
                    self._place(i, (centre[0] + across[0] * c, centre[1] + across[1] * c), direction, src)
            else:
                # members still free go next to the placed ones, alternately on both sides
                for k, i in enumerate(free):
                    c = (k // 2 + 1) * size(i) * (1 if k % 2 == 0 else -1)
                    self._place(i, (centre[0] + across[0] * c, centre[1] + across[1] * c), direction)
        if len(ids) == 1 and ids[0] in self.pos:
            centre = self.pos[ids[0]]  # where it went: the walk goes on from there
        return centre, ids

    def _touches(self, nodes) -> bool:
        return any(self.nid(p) in self.pos and self.nid(p) not in self.side
                   for node in nodes for p in node.participants)

    def chain(self, chain: Chain, head=None, first_dir=None, first_step=None):
        """Walk a chain; it starts at its first node that is drawn already (forwards and backwards from there)."""
        h = head or self.head0
        dirs = []
        for k, link in enumerate(chain.links):
            d = first_dir if (k == 0 and first_dir is not None) else _rot(h, link.turn)
            dirs.append(d)
            # the generator measured turns against the drawn step, rounded to four directions
            h = _dir4(d)
        nodes = chain.nodes
        steps = [STEP / 2 if (nodes[k].open or nodes[k + 1].open) else STEP for k in range(len(dirs))]
        if chain.wrap and chain.wrap[0] == "serpentine":
            # a turn back in a wrapped chain goes to the next row (the drawn step was mostly sideways)
            row = _rot(self.head0, "R")
            for k, link in enumerate(chain.links):
                if link.turn == "B":
                    dirs[k] = row
                    steps[k] = STEP * 0.45
        if first_step is not None and steps:
            steps[0] = first_step
        anchor = next((k for k, node in enumerate(nodes) if self._touches([node])), 0)
        points: list = [None] * len(nodes)
        ids: list = [None] * len(nodes)
        d0 = dirs[min(anchor, len(dirs) - 1)]
        points[anchor], ids[anchor] = self._node_point(nodes[anchor], d0)
        for k in range(anchor + 1, len(nodes)):
            at = (points[k - 1][0] + dirs[k - 1][0] * steps[k - 1], points[k - 1][1] + dirs[k - 1][1] * steps[k - 1])
            points[k], ids[k] = self._node_point(nodes[k], dirs[k - 1], at, points[k - 1])
            self._segment(points[k - 1], points[k])
        for k in range(anchor - 1, -1, -1):
            at = (points[k + 1][0] - dirs[k][0] * steps[k], points[k + 1][1] - dirs[k][1] * steps[k])
            points[k], ids[k] = self._node_point(nodes[k], dirs[k], at, points[k + 1])
            self._segment(points[k], points[k + 1])
        for k, link in enumerate(chain.links):
            self.links.append((link.reactions, ids[k], ids[k + 1], dirs[k], points[k], points[k + 1]))

    def _segment(self, p, q):
        if math.dist(p, q) > 1e-6:
            self.segments.append((p, q))
            self.mids.append(((p[0] + q[0]) / 2, (p[1] + q[1]) / 2))

    def star(self, star: Star):
        self._node_point(star.hub, self.head0)
        n = len(star.spokes)
        if not n:
            return
        a0 = math.atan2(*reversed(_rot(self.head0, star.spokes[0].links[0].turn)))
        arc = math.radians(min(star.arc, 360))
        step = arc / max(n - 1, 1) if star.arc < 360 else 2 * math.pi / n
        sign = 1 if star.order == "cw" else -1  # cw on screen: increasing angle with y down
        # spokes far enough from the hub that their first nodes are about 0.8 STEP apart
        radius = max(STEP, 0.8 * STEP / max(step, 1e-3))
        for i, spoke in enumerate(star.spokes):
            a = a0 + sign * i * step
            self.chain(spoke, first_dir=(math.cos(a), math.sin(a)), first_step=radius)

    def react(self, react: React):
        r = react.reaction
        ids = [self.nid(p) for p in r.unsided]
        placed = [self.pos[i] for i in ids if i in self.pos]
        centre = ((sum(x for x, _ in placed) / len(placed), sum(y for _, y in placed) / len(placed)) if placed
                  else self._new_component())
        free = [i for i in ids if i not in self.pos]
        for j, i in enumerate(free):
            a = 2 * math.pi * j / max(len(free), 1) + math.pi / 2
            self._place(i, (centre[0] + math.cos(a) * STEP / 2, centre[1] + math.sin(a) * STEP / 2), self.head0)
        self.links.append(([r], ids, [], self.head0, centre, centre))

    def ring(self, ring: Ring) -> None:
        """A ring, the cycle of the drawing: its metabolites go on a circle in the order and sense of the ring,
        starting in the direction of the first step; the other statements of the ring are then walked from them.

        Without `circle_rings`, the ring is walked like a chain unless the walk of a run hits itself (the four
        directions of the turn letters cannot express diagonal steps, e.g. a ring drawn as a diamond)."""
        runs = [t for t in ring.statements if isinstance(t, Chain) and not t.chord
                and not any(n.group or n.open for n in t.nodes)]
        cycle, seen = [], set()
        for run in runs:
            for node in run.nodes:
                i = self.nid(node.participants[0])
                if i not in seen:
                    seen.add(i)
                    cycle.append(i)
        collides = False
        for run in runs:
            h, point, points = self.head0, (0, 0), {(0, 0)}
            for link in run.links:
                d = _rot(h, link.turn)
                point = (point[0] + d[0], point[1] + d[1])
                collides |= point in points
                points.add(point)
                h = d
        if (collides or self.circle_rings) and len(cycle) >= 3 and runs:
            n = len(cycle)
            radius = STEP * (1.0 if self.circle_rings else 0.8) / (2 * math.sin(math.pi / n))
            sign = 1 if ring.sense == "cw" else -1  # cw on screen: increasing angle (y down)
            d = _rot(self.head0, runs[0].links[0].turn)
            theta0 = math.atan2(d[1], d[0]) - sign * (math.pi / 2 + math.pi / n)
            angles = [theta0 + sign * 2 * math.pi * k / n for k in range(n)]
            placed = [k for k, i in enumerate(cycle) if i in self.pos]
            if placed:
                k = placed[0]
                x, y = self.pos[cycle[k]]
                centre = (x - radius * math.cos(angles[k]), y - radius * math.sin(angles[k]))
            else:
                start = self._new_component()
                centre = (start[0] + radius - radius * math.cos(angles[0]), start[1] - radius * math.sin(angles[0]))
            for i, a in zip(cycle, angles):
                if i not in self.pos:
                    self.pos[i] = (centre[0] + radius * math.cos(a), centre[1] + radius * math.sin(a))
            self._ring_inner(ring, cycle, angles, centre, radius)
        for t in ring.statements:
            self.react(t) if isinstance(t, React) else self.chain(t)

    def _ring_inner(self, ring: Ring, cycle, angles, centre, radius) -> None:
        """Other metabolites of a ring that react with two or more metabolites of the circle: a metabolite
        bridging distant ones goes inside, at their centre (a chord); one next to two neighbouring ones goes
        outside, beyond their edge (a triangle on the circle)."""
        index = {i: k for k, i in enumerate(cycle)}
        n = len(cycle)
        partners = defaultdict(set)
        reactions = []
        for t in ring.statements:
            if isinstance(t, React):
                reactions.append(t.reaction)
            else:
                reactions += [r for link in t.links for r in link.reactions]
        for r in reactions:
            ids = {self.nid(p) for p in r.participants} - self.side
            for i in ids - set(index):
                partners[i] |= ids & set(index)
        for i, ps in sorted(partners.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            if i in self.pos or len(ps) < 2:
                continue
            ks = sorted(index[p] for p in ps)
            gaps = [(ks[(j + 1) % len(ks)] - ks[j]) % n for j in range(len(ks))]
            if max(gaps) >= n - 1 - (len(ks) - 2):  # all on one arc of neighbouring nodes: outside
                start = ks[(gaps.index(max(gaps)) + 1) % len(ks)]
                mid = start + (n - max(gaps)) / 2
                a = angles[0] + (angles[1] - angles[0]) * mid if n > 1 else 0.0
                point = (centre[0] + (radius + 0.6 * STEP) * math.cos(a), centre[1] + (radius + 0.6 * STEP) * math.sin(a))
            else:
                point = _centroid([self.pos[p] for p in ps], centre)
            self._place(i, point, self.head0)

    def run(self):
        """Process the statements in their order, but each one only once it touches the drawn part (if any does):
        the generator writes reactions with several metabolites last, though they often join two chains."""
        units = list(self.doc.statements)

        def nodes_of(u):
            if isinstance(u, Star):
                return [u.hub]
            if isinstance(u, Ring):
                return [n for t in u.statements for n in nodes_of(t)]
            if isinstance(u, React):
                from .syntax import Node as PNode
                return [PNode([p]) for p in u.reaction.unsided]
            return u.nodes

        todo = list(units)
        while todo:
            k = next((k for k, u in enumerate(todo) if self._touches(nodes_of(u))), 0)
            u = todo.pop(k)
            if isinstance(u, Star):
                self.star(u)
            elif isinstance(u, Ring):
                self.ring(u)
            elif isinstance(u, React):
                self.react(u)
            else:
                self.chain(u)
        return self.pos


def _cross(p1, p2, q1, q2) -> bool:
    """Whether two segments cross properly (not at a common end)."""
    def orient(a, b, c):
        v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return 0 if abs(v) < 1e-9 else (1 if v > 0 else -1)
    if min(math.dist(p1, q1), math.dist(p1, q2), math.dist(p2, q1), math.dist(p2, q2)) < 1e-6:
        return False
    return (orient(p1, p2, q1) * orient(p1, p2, q2) < 0) and (orient(q1, q2, p1) * orient(q1, q2, p2) < 0)


def _dir4(v):
    if abs(v[0]) >= abs(v[1]):
        return (1, 0) if v[0] >= 0 else (-1, 0)
    return (0, 1) if v[1] >= 0 else (0, -1)


def _sides(r, from_ids, to_ids, node_ids):
    """Main nodes of the substrate side and of the product side (or the two written sides)."""
    def ids(ps):
        return [node_ids.get(p.node, p.node) for p in ps]
    if r.substrates is not None and r.reversible is not True:
        return ids(r.substrates), ids(r.products)
    if r.left or r.right:
        return ids(r.left), ids(r.right)
    return from_ids, to_ids


def _centroid(points, default):
    points = [p for p in points if p is not None]
    if not points:
        return default
    return (sum(x for x, _ in points) / len(points), sum(y for _, y in points) / len(points))


def place_enzymes(dm: DraftMap, pos: dict, links, node_ids) -> None:
    """Enzyme nodes between the two sides of their reaction; the reactions of a fan side by side."""
    for reactions, from_ids, to_ids, d, p0, p1 in links:
        k = len(reactions)
        for j, r in enumerate(reactions):
            eid = next(e for e, rr in dm.pwl_reactions if rr is r)
            a, b = _sides(r, from_ids, to_ids, node_ids)
            ca = _centroid([pos.get(i) for i in a], p0)
            cb = _centroid([pos.get(i) for i in b], p1)
            mid = ((ca[0] + cb[0]) / 2, (ca[1] + cb[1]) / 2)
            if not b and not r.unsided:  # open end: closer to the drawn side
                mid = ((ca[0] * 2 + p1[0]) / 3, (ca[1] * 2 + p1[1]) / 3)
            axis = _unit((cb[0] - ca[0], cb[1] - ca[1])) if (a and b) else _unit(d)
            across = (-axis[1], axis[0])
            n = dm.nodes[eid]
            size = abs(across[0]) * n.w + abs(across[1]) * n.h + 2 * GAP
            off = (j - (k - 1) / 2) * size
            pos[eid] = (mid[0] + across[0] * off, mid[1] + across[1] * off)


def hint_positions(doc: Document, dm: DraftMap, node_ids: dict[str, str] | None = None,
                   circle_rings: bool = True) -> dict:
    """Positions of the enzyme and main-metabolite nodes from the drawing hints of `doc`
    (without the secondary metabolites of ``dm.side``)."""
    node_ids = node_ids or {}
    turtle = _Turtle(doc, dm, node_ids, circle_rings)
    pos = turtle.run()
    side = getattr(dm, "side", set())
    for i, n in dm.nodes.items():  # nodes no statement reaches (cannot happen in a valid document)
        if n.kind == "main" and i not in pos and i not in side:
            pos[i] = turtle._new_component()
    place_enzymes(dm, pos, turtle.links, node_ids)
    return pos


# ---------------------------------------------------------------------------------- cofactors and helpers


def _axis(dm: DraftMap, eid: str, r, pos, node_ids):
    """Unit vector from the substrate side to the product side of a reaction, and its position."""
    a, b = _sides(r, [], [], node_ids)
    e = pos[eid]
    ca = _centroid([pos.get(i) for i in a], None)
    cb = _centroid([pos.get(i) for i in b], None)
    if ca is not None and cb is not None:
        return _unit((cb[0] - ca[0], cb[1] - ca[1]))
    if ca is not None:
        return _unit((e[0] - ca[0], e[1] - ca[1]))
    if cb is not None:
        return _unit((cb[0] - e[0], cb[1] - e[1]))
    return (0.0, 1.0)


def place_cofactors(dm: DraftMap, pos: dict, node_ids=None) -> None:
    """Cofactor nodes around their enzyme: consumed upstream, produced downstream, on the side ``^L``/``^R``.

    Without a side, consumed cofactors go to the right of the reaction axis and
    produced ones to the left, as most often on the BRENDA map. Several
    cofactors of a kind are fanned out.
    """
    node_ids = node_ids or {}
    by_enzyme = defaultdict(list)
    for e in dm.edges:
        if e.kind == "cofactor":
            by_enzyme[e.enzyme].append(e)
    for eid, r in dm.pwl_reactions:
        edges = by_enzyme.get(eid)
        if not edges or eid not in pos:
            continue
        u = _axis(dm, eid, r, pos, node_ids)
        right = (-u[1], u[0])  # screen: to the right of the axis
        ex, ey = pos[eid]
        for role in ("consumed", "produced"):
            items = [e for e in edges if (e.arrow_enzyme if role == "consumed" else e.arrow_met)]
            if r.cofactor_side:
                sign = 1 if r.cofactor_side == "R" else -1
            else:
                sign = 1 if role == "consumed" else -1
            for j, e in enumerate(items):
                # angle from the axis: upstream (consumed) around 135 degrees, downstream around 45 degrees
                base = 125 if role == "consumed" else 55
                spread = (j - (len(items) - 1) / 2) * 32 * (1 if role == "consumed" else -1)
                a = math.radians(base + spread)
                dist = COFACTOR_DISTANCE + dm.nodes[e.met].w / 4
                v = (math.cos(a) * u[0] + math.sin(a) * right[0] * sign,
                     math.cos(a) * u[1] + math.sin(a) * right[1] * sign)
                # cos(a) < 0 for upstream: negative along the axis
                pos[e.met] = (ex + v[0] * dist, ey + v[1] * dist)


def cofactor_overlaps(dm: DraftMap, pos: dict) -> int:
    """Number of node pairs with a cofactor that overlap."""
    ids = [i for i in pos if i in dm.nodes]
    count = 0
    for a_i, a in enumerate(ids):
        na = dm.nodes[a]
        for b in ids[a_i + 1:]:
            nb = dm.nodes[b]
            if na.kind != "cofactor" and nb.kind != "cofactor":
                continue
            if (abs(pos[a][0] - pos[b][0]) < (na.w + nb.w) / 2 + 4
                    and abs(pos[a][1] - pos[b][1]) < (na.h + nb.h) / 2 + 4):
                count += 1
    return count


def spread_for_cofactors(dm: DraftMap, pos: dict, node_ids=None, factors=(1.0, 1.1, 1.2, 1.3)) -> dict:
    """Place the cofactors (:func:`place_cofactors`); if they overlap other nodes, spread the drawing (enzymes
    and main metabolites, about its centre) by the next factor to make room, and keep the first factor without
    overlaps (or the one with the fewest)."""
    main = {i: p for i, p in pos.items() if dm.nodes[i].kind != "cofactor"}
    cx, cy = _centroid(list(main.values()), (0.0, 0.0))
    best, key = None, None
    for f in factors:
        p = {i: (cx + (x - cx) * f, cy + (y - cy) * f) for i, (x, y) in main.items()}
        place_cofactors(dm, p, node_ids)
        k = cofactor_overlaps(dm, p)
        if key is None or k < key:
            best, key = p, k
        if not k:
            break
    return best


def _merged_point(dm: DraftMap, pos: dict, nodes, cof):
    """Where a merged cofactor node goes: the centre of the nodes or one of their places, the first one that
    overlaps no other node and whose links cross no more links than those of the separate nodes; else None."""
    others = [(pos[e.enzyme], pos[e.met]) for e in dm.edges
              if e.enzyme in pos and e.met in pos and e.met not in nodes]
    old = sum(_cross(pos[i], pos[cof[i]], a, b) for i in nodes for a, b in others)
    n = dm.nodes[min(nodes)]
    for point in [_centroid([pos[i] for i in nodes], None)] + sorted(pos[i] for i in nodes):
        if any(abs(point[0] - x) < (n.w + dm.nodes[o].w) / 2 + 4 and abs(point[1] - y) < (n.h + dm.nodes[o].h) / 2 + 4
               for o, (x, y) in pos.items() if o not in nodes and o in dm.nodes):
            continue
        new = sum(_cross(point, pos[cof[i]], a, b) for i in nodes for a, b in others)
        if new <= old:
            return point
    return None


def share_cofactors(dm: DraftMap, pos: dict, reach: float = 0.6 * STEP) -> dict:
    """Merge cofactor nodes of the same compound of nearby reactions (nodes at most `reach` apart) into one node
    linked to all of these enzymes, at their centre; a node is linked to an enzyme only once."""
    pos = dict(pos)
    cof = {e.met: e.enzyme for e in dm.edges if e.kind == "cofactor" and e.met in pos}
    by_label = defaultdict(list)
    for i in cof:
        by_label[dm.nodes[i].label].append(i)
    for ids in by_label.values():
        if len(ids) < 2:
            continue
        cluster = {i: {i} for i in ids}
        pairs = sorted((math.dist(pos[a], pos[b]), a, b) for k, a in enumerate(ids) for b in ids[k + 1:])
        for d, a, b in pairs:
            if d > reach or cluster[a] is cluster[b]:
                continue
            ca, cb = cluster[a], cluster[b]
            if {cof[i] for i in ca} & {cof[i] for i in cb}:
                continue
            # all members close to each other, not a long chain of merges
            if any(math.dist(pos[x], pos[y]) > reach * 1.5 for x in ca for y in cb):
                continue
            ca |= cb
            for i in cb:
                cluster[i] = ca
        for c in {id(c): c for c in cluster.values()}.values():
            if len(c) < 2:
                continue
            keep = min(c)
            point = _merged_point(dm, pos, c, cof)
            if point is None:
                continue
            pos[keep] = point
            for e in dm.edges:
                if e.kind == "cofactor" and e.met in c:
                    e.met = keep
            for i in c - {keep}:
                del dm.nodes[i]
                del pos[i]
    return pos


def _apply(dm: DraftMap, pos: dict) -> None:
    for i, (x, y) in pos.items():
        dm.nodes[i].x, dm.nodes[i].y = float(x), float(y)


def _main_graph(dm: DraftMap, directed=False):
    import networkx as nx

    g = nx.DiGraph() if directed else nx.Graph()
    side = getattr(dm, "side", set())
    g.add_nodes_from(i for i, n in dm.nodes.items() if n.kind != "cofactor" and i not in side)
    for e in dm.edges:
        if e.kind != "main" or e.met in side:
            continue
        if directed and e.arrow_enzyme and not e.arrow_met:
            g.add_edge(e.met, e.enzyme)
        else:
            g.add_edge(e.enzyme, e.met)
    return g


def scale_to_brenda(dm: DraftMap, pos: dict) -> dict:
    """Scale positions so that the median enzyme-metabolite distance is half of :data:`STEP`."""
    lengths = sorted(math.dist(pos[e.enzyme], pos[e.met]) for e in dm.edges
                     if e.kind == "main" and e.enzyme in pos and e.met in pos)
    if not lengths or lengths[len(lengths) // 2] < 1e-6:
        return pos
    f = (STEP / 2) / lengths[len(lengths) // 2]
    return {i: (x * f, y * f) for i, (x, y) in pos.items()}


def remove_overlaps(dm: DraftMap, pos: dict, rounds: int = 60, cofactors: bool = False) -> dict:
    """Push overlapping boxes apart along the smaller overlap; small moves keep the shape.

    Without `cofactors`, only enzyme and main nodes; with it, all nodes, and a
    cofactor that overlaps a main or enzyme node moves alone.
    """
    ids = [i for i in pos if cofactors or dm.nodes[i].kind != "cofactor"]
    pos = dict(pos)
    for _ in range(rounds):
        moved = False
        for a_i in range(len(ids)):
            for b_i in range(a_i + 1, len(ids)):
                a, b = ids[a_i], ids[b_i]
                na, nb = dm.nodes[a], dm.nodes[b]
                dx, dy = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
                gap = GAP / 2 if na.kind != "cofactor" and nb.kind != "cofactor" else 4.0
                ox = (na.w + nb.w) / 2 + gap - abs(dx)
                oy = (na.h + nb.h) / 2 + gap - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    wa = 0.0 if (na.kind != "cofactor" and nb.kind == "cofactor") else (
                        1.0 if (na.kind == "cofactor" and nb.kind != "cofactor") else 0.5)
                    wb = 1.0 - wa
                    if ox < oy:
                        sx = (1 if dx >= 0 else -1) * ox
                        pos[a] = (pos[a][0] - sx * wa, pos[a][1])
                        pos[b] = (pos[b][0] + sx * wb, pos[b][1])
                    else:
                        sy = (1 if dy >= 0 else -1) * oy
                        pos[a] = (pos[a][0], pos[a][1] - sy * wa)
                        pos[b] = (pos[b][0], pos[b][1] + sy * wb)
        if not moved:
            break
    return pos


def straighten(dm: DraftMap, pos: dict, tolerance: float = 20.0, rounds: int = 3) -> dict:
    """Make main links that are within `tolerance` degrees of horizontal or vertical exactly so."""
    pos = dict(pos)
    edges = [(e.enzyme, e.met) for e in dm.edges if e.kind == "main" and e.enzyme in pos and e.met in pos]
    t = math.tan(math.radians(tolerance))
    for _ in range(rounds):
        for a, b in edges:
            dx, dy = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
            if abs(dx) <= t * abs(dy):
                x = (pos[a][0] + pos[b][0]) / 2
                pos[a], pos[b] = (x, pos[a][1]), (x, pos[b][1])
            elif abs(dy) <= t * abs(dx):
                y = (pos[a][1] + pos[b][1]) / 2
                pos[a], pos[b] = (pos[a][0], y), (pos[b][0], y)
    return pos


# ---------------------------------------------------------------------------------- candidates


def _graphviz(dm: DraftMap, pos: dict, prog: str, flow: str, start=None, lengths=False) -> dict:
    """A Graphviz layout of the enzyme and main nodes (nodes as boxes of their size)."""
    from ..draft_map.relayout import graphviz_positions

    ids = list(_main_graph(dm).nodes)
    nodes = []
    for i in ids:
        n = dm.nodes[i]
        attrs = dict(width=f"{n.w / 72:.3f}", height=f"{n.h / 72:.3f}")
        if start is not None and prog != "dot":
            attrs["pos"] = f"{start[i][0]:.1f},{-start[i][1]:.1f}"
        nodes.append([i, attrs])
    g = _main_graph(dm, directed=prog == "dot")
    edges = [list(e) for e in g.edges()]
    graph = {"inputscale": "72", "start": "1"}
    if prog == "dot":
        graph.update(rankdir={"down": "TB", "up": "BT", "right": "LR", "left": "RL"}.get(flow, "TB"),
                     nodesep="0.35", ranksep="0.9")
    elif prog == "neato":
        graph.update(overlap="vpsc", sep="+6", len="3.0")
    elif prog == "fdp":
        graph.update(K="2.5", overlap="vpsc")
    elif prog == "sfdp":
        graph.update(K="2.5", overlap="vpsc")
    spec = dict(prog=prog, directed=prog == "dot", graph=graph, node=dict(shape="box", fixedsize="true", label=""),
                nodes=nodes, edges=edges)
    if lengths and start is not None:
        # edge lengths of the start layout (inches): stress majorization keeps its shape
        spec["edge_len"] = [math.dist(start[a], start[b]) / 72 for a, b in edges]
        spec["graph"]["mode"] = "major"
        return _graphviz_with_lengths(spec)
    out = graphviz_positions(spec)
    return {i: (x, -y) for i, (x, y) in out.items()}


def _segments(dm: DraftMap, pos: dict):
    return [(e.enzyme, e.met) for e in dm.edges if e.kind == "main" and e.enzyme in pos and e.met in pos]


def crossing_nodes(dm: DraftMap, pos: dict) -> tuple[int, set[str]]:
    """Number of crossings of the main links and the nodes at their ends."""
    seg = _segments(dm, pos)
    count, nodes = 0, set()
    for i in range(len(seg)):
        a, b = seg[i]
        for j in range(i + 1, len(seg)):
            c, d = seg[j]
            if _cross(pos[a], pos[b], pos[c], pos[d]):
                count += 1
                nodes |= {a, b, c, d}
    return count, nodes


def untangle(dm: DraftMap, pos: dict, rounds: int = 3) -> dict:
    """Graphviz ``neato`` with all nodes pinned except those at the ends of crossing links.

    The free nodes move to positions that fit the lengths of their links
    (taken from the current drawing), so crossings are resolved locally and
    the rest of the drawing keeps its shape. A round is kept only if it
    reduces the crossings.
    """
    import pygraphviz as pgv

    pos = dict(pos)
    ids = [i for i in _main_graph(dm).nodes if i in pos]
    edges = [list(e) for e in _main_graph(dm).edges()]
    lengths = sorted(math.dist(pos[a], pos[b]) for a, b in edges) or [STEP / 2]
    typical = lengths[len(lengths) // 2]
    for _ in range(rounds):
        count, free = crossing_nodes(dm, pos)
        if not count:
            break
        g = pgv.AGraph(directed=False, strict=False)
        g.graph_attr.update(inputscale="72", overlap="true", mode="major", start="1")
        g.node_attr.update(shape="box", fixedsize="true", label="")
        for i in ids:
            n = dm.nodes[i]
            g.add_node(i, width=f"{n.w / 72:.3f}", height=f"{n.h / 72:.3f}",
                       pos=f"{pos[i][0]:.1f},{-pos[i][1]:.1f}" + ("" if i in free else "!"))
        for a, b in edges:
            g.add_edge(a, b, len=f"{typical / 72:.3f}")
        g.layout(prog="neato")
        new = {i: (float(g.get_node(i).attr["pos"].split(",")[0].rstrip("!")),
                   -float(g.get_node(i).attr["pos"].split(",")[1].rstrip("!"))) for i in ids}
        # Graphviz translates the drawing: move it back by a pinned node
        anchor = next((i for i in ids if i not in free), None)
        if anchor is not None:
            dx, dy = pos[anchor][0] - new[anchor][0], pos[anchor][1] - new[anchor][1]
            new = {i: (x + dx, y + dy) for i, (x, y) in new.items()}
        new = remove_overlaps(dm, {**pos, **new})
        if crossing_nodes(dm, new)[0] >= count:
            break
        pos = new
    return pos


def _graphviz_with_lengths(spec: dict) -> dict:
    import pygraphviz as pgv

    g = pgv.AGraph(directed=False, strict=False)
    g.graph_attr.update(spec["graph"])
    g.node_attr.update(spec["node"])
    for name, attrs in spec["nodes"]:
        g.add_node(name, **attrs)
    for (a, b), length in zip(spec["edges"], spec["edge_len"]):
        g.add_edge(a, b, len=f"{max(length, 0.3):.3f}")
    g.layout(prog=spec["prog"])
    return {name: (float(g.get_node(name).attr["pos"].split(",")[0]),
                   -float(g.get_node(name).attr["pos"].split(",")[1])) for name, _ in spec["nodes"]}


def _networkx(dm: DraftMap, method: str, start=None, seed: int = 1) -> dict:
    import networkx as nx

    g = _main_graph(dm)
    if method == "kamada_kawai":
        pos = nx.kamada_kawai_layout(g, scale=STEP * max(1.0, math.sqrt(len(g)) / 2))
    elif start is not None:
        pos = nx.spring_layout(g, pos={i: start[i] for i in g}, k=STEP / 2, iterations=25, seed=seed)
    else:
        pos = nx.spring_layout(g, k=1.5 / math.sqrt(max(len(g), 1)), iterations=200, seed=seed)
    return {i: (float(x), float(y)) for i, (x, y) in pos.items()}


def _circle(dm: DraftMap, pos: dict) -> dict | None:
    from ..draft_map.relayout import circle_layout

    s = copy.deepcopy(dm)
    for i in getattr(dm, "side", set()):
        del s.nodes[i]
    s.edges = [e for e in s.edges if e.met in s.nodes]
    _apply(s, pos)
    for n in s.nodes.values():
        n.status = "map"  # circle_layout looks for cycles among nodes of the BRENDA map
    if not circle_layout(s):
        return None
    return {i: (n.x, n.y) for i, n in s.nodes.items() if n.kind != "cofactor"}


def available(method: str) -> bool:
    """Whether the packages of a layout method are installed."""
    from importlib.util import find_spec

    needs = ["networkx"]
    if method in ("dot", "neato", "fdp", "sfdp", "hints+neato", "hints+spring+neato"):
        needs.append("pygraphviz")
    if method in ("kamada_kawai", "circle"):
        needs.append("scipy")
    return all(find_spec(n) is not None for n in needs)


def layout_candidates(doc: Document, node_ids: dict[str, str] | None = None, methods=METHODS,
                      fill: str = "#dedede", side: bool = True, circle_rings: bool = True,
                      spread: bool = True, share: bool = True) -> dict[str, DraftMap]:
    """The draft maps of `doc` laid out by each method (only those whose packages are installed).

    Steps after the layout of enzymes and main metabolites (each can be switched off):
    `side` lays out secondary metabolites afterwards (:func:`side_metabolites`),
    `circle_rings` puts rings on circles, `spread` makes room for cofactors
    (:func:`spread_for_cofactors`) and `share` merges the same cofactor of
    nearby reactions (:func:`share_cofactors`).
    """
    node_ids = node_ids or {}
    base = to_draftmap(doc, node_ids, fill)
    base.side = side_metabolites(doc, base, node_ids) if side else set()
    hints = hint_positions(doc, base, node_ids, circle_rings)
    flow = doc.flow or "down"
    jobs = {}
    for m in methods:
        if not available(m):
            continue
        if m == "hints":
            jobs[m] = lambda: hints
        elif m == "hints+neato":
            jobs[m] = lambda: untangle(base, remove_overlaps(base, hints))
        elif m == "hints+spring":
            jobs[m] = lambda: straighten(base, scale_to_brenda(base, _networkx(base, "spring", start=hints)))
        elif m == "hints+spring+neato":
            jobs[m] = lambda: untangle(base, remove_overlaps(base, straighten(
                base, scale_to_brenda(base, _networkx(base, "spring", start=hints)))))
        elif m in ("dot", "neato", "fdp", "sfdp"):
            jobs[m] = lambda m=m: _graphviz(base, hints, m, flow)
        elif m in ("kamada_kawai", "spring"):
            jobs[m] = lambda m=m: _networkx(base, m)
        elif m == "circle":
            jobs[m] = lambda: _circle(base, hints)
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {m: pool.submit(f) for m, f in jobs.items()}
        results = {}
        for m, f in futures.items():
            try:
                results[m] = f.result()
            except Exception as e:  # a layout program can fail on odd graphs; the others still count
                base.log.append(f"{m} failed: {e}")
    out = {}
    for m, pos in results.items():
        if pos is None:
            continue
        pos = dict(pos)
        if m != "hints":
            pos = scale_to_brenda(base, pos)
        pos = remove_overlaps(base, pos)
        if base.side:
            pos = remove_overlaps(base, place_side_metabolites(base, pos, base.side))
        dm = copy.deepcopy(base)
        dm.pwl_reactions = base.pwl_reactions  # the same reaction objects (identity is used)
        if spread:
            pos = spread_for_cofactors(dm, pos, node_ids)
        else:
            place_cofactors(dm, pos, node_ids)
        if share:
            pos = share_cofactors(dm, pos)
        pos = remove_overlaps(dm, pos, rounds=200, cofactors=True)
        _apply(dm, pos)
        fit_subsystem(dm)
        dm.subsystems[doc.id].layout = m
        out[m] = dm
    return out


def fit_subsystem(dm: DraftMap) -> None:
    """Translate the drawing to the origin and fit the subsystem box around it (title on top)."""
    nodes = list(dm.nodes.values())
    if not nodes:
        return
    x0 = min(n.x - n.w / 2 for n in nodes)
    y0 = min(n.y - n.h / 2 for n in nodes)
    x1 = max(n.x + n.w / 2 for n in nodes)
    y1 = max(n.y + n.h / 2 for n in nodes)
    for s in dm.subsystems.values():
        top = s.title_h / 2 + TITLE_GAP
        width = max(x1 - x0 + 2 * PADDING, s.title_w + 2 * PADDING)
        s.x, s.y, s.w, s.h = 0.0, 0.0, width, (y1 - y0) + top + PADDING
        dx = (width - (x1 - x0)) / 2 - x0
        dy = top - y0
        for n in nodes:
            n.x += dx
            n.y += dy
        for e in dm.edges:
            if getattr(e, "route", None):
                e.route = [(x + dx, y + dy) for x, y in e.route]


def metrics(dm: DraftMap) -> dict:
    """Crossings and overlaps (:func:`ecmapper.draft_map.layout.internal_metrics`)."""
    from ..draft_map.layout import internal_metrics

    return internal_metrics(dm)


def choose(candidates: dict[str, DraftMap], margin: float = 0.5) -> str:
    """The layout to draw: the hint layout with the fewest crossings and overlaps, unless another one has at
    least `margin` fewer main crossings and no overlapping nodes."""
    scores = {m: metrics(dm) for m, dm in candidates.items()}

    def key(m):
        s = scores[m]
        return (s["overlapping_nodes"] > 0, s["crossings_main"], s["crossings_all"], s["overlapping_nodes"])

    hints = [m for m in candidates if m in HINT_METHODS]
    best = min(hints, key=key) if hints else min(candidates, key=key)
    others = [m for m in candidates if m not in HINT_METHODS and scores[m]["overlapping_nodes"] == 0]
    if others:
        other = min(others, key=key)
        if scores[other]["crossings_main"] < (1 - margin) * scores[best]["crossings_main"]:
            return other
    return best


def route(dm: DraftMap) -> int:
    """Route main edges that pass through other nodes around them (Graphviz), see :mod:`ecmapper.draft_map.routes`."""
    try:
        from ..draft_map.routes import route_edges
    except ImportError:
        return 0
    return route_edges(dm)


__all__ = ["BASELINES", "HINT_METHODS", "METHODS", "STEP", "choose", "fit_subsystem", "hint_positions",
           "layout_candidates", "metrics", "place_cofactors", "place_enzymes", "remove_overlaps", "route",
           "scale_to_brenda", "straighten", "to_draftmap"]
