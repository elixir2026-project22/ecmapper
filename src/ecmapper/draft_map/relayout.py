"""Step 3b: lay out the subsystems with added reactions again (networkx and Graphviz layouts, the best one is kept).

The completion (step 3) places added reactions between their drawn
metabolites, which gives long, crossing edges where many reactions are
added. For each subsystem with added reactions, several layouts of its
enzyme and main metabolite nodes are computed, each subsystem is condensed
as in step 4, and the layout with the fewest crossings of main-metabolite
links is kept:

- ``spring``: networkx spring layout (Fruchterman-Reingold) of the added
  nodes, the nodes of the BRENDA map fixed;
- ``fdp``: the same with Graphviz ``fdp``;
- ``dot``: Graphviz ``dot`` (layered, Sugiyama), along the reaction directions;
- ``neato``: Graphviz ``neato`` (stress majorisation), all nodes free,
  starting from the current drawing, overlaps removed with VPSC.

``spring`` and ``fdp`` keep the drawing of the BRENDA map; ``dot`` and
``neato`` are only taken if they have at least `margin` fewer crossings than
the best of the current drawing, ``spring`` and ``fdp``. Cofactors move with
their enzyme node. Subsystems named in `circle` get the cycle of their main
metabolites on a circle (:func:`circle_layout`), whatever the crossings.
Graphviz layouts need pygraphviz (2.0 or later ships Graphviz), networkx
layouts networkx; missing ones are skipped.

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
KEEPING = ("current", "spring", "fdp")
REDRAWING = ("dot", "neato")
METHODS = KEEPING + REDRAWING


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
    "neato": dict(prog="neato", pin=False, len="2.0", overlap="vpsc", sep="+6"),
}


def graphviz_spec(s: DraftMap, prog: str, pin: bool = True, **graph_attr) -> tuple[dict, dict, dict]:
    """The graph of the enzyme and main metabolite nodes for :mod:`.graphviz_worker`, and the node boxes.

    ``neato`` and ``fdp`` start from the current drawing; with `pin`, the
    nodes of the BRENDA map keep their positions. Graphviz reads input
    positions in inches unless ``inputscale`` is set, and ``overlap=true``
    keeps overlaps (removing them would move pinned nodes; step 4 removes
    them anyway).
    """
    box, owner = _groups(s)
    nodes = []
    for i, b in box.items():
        attrs = dict(width=f"{(b[2] - b[0]) / 72:.3f}", height=f"{(b[3] - b[1]) / 72:.3f}")
        if prog != "dot":
            pinned = "!" if pin and s.nodes[i].status in MAP_STATUS else ""
            attrs["pos"] = f"{(b[0] + b[2]) / 2:.1f},{-(b[1] + b[3]) / 2:.1f}{pinned}"  # Graphviz: y upwards
        nodes.append([i, attrs])
    spec = dict(prog=prog, directed=prog == "dot",
                graph={"start": "random1", "inputscale": "72", "overlap": "true", **graph_attr},
                node=dict(shape="box", fixedsize="true", label=""), nodes=nodes, edges=_directed(s))
    return spec, box, owner


def graphviz_positions(spec: dict) -> dict[str, list[float]]:
    """Run a Graphviz layout in a process of its own (see the module doc)."""
    result = subprocess.run([sys.executable, str(WORKER)], input=json.dumps(spec), capture_output=True, text=True,
                            check=False)
    if result.returncode:
        raise RuntimeError(f"Graphviz {spec['prog']} failed: {result.stderr.strip()}")
    return json.loads(result.stdout)


def apply_graphviz(s: DraftMap, positions: dict[str, list[float]], box, owner) -> None:
    _move(s, {i: (x, -y) for i, (x, y) in positions.items()}, box, owner)


def graphviz_layout(s: DraftMap, prog: str, pin: bool = True, **graph_attr) -> None:
    """A Graphviz layout of the enzyme and main metabolite nodes (pygraphviz, in a process of its own)."""
    spec, box, owner = graphviz_spec(s, prog, pin, **graph_attr)
    apply_graphviz(s, graphviz_positions(spec), box, owner)


LAYOUTS = {"spring": spring_layout, "circle": circle_layout,
           **{m: (lambda s, a=a: graphviz_layout(s, **a)) for m, a in GRAPHVIZ.items()}}


def available_methods(methods=METHODS) -> list[str]:
    """The methods whose packages are installed."""
    result = []
    for m in methods:
        package = {"spring": "networkx", "circle": "networkx", "fdp": "pygraphviz", "dot": "pygraphviz",
                   "neato": "pygraphviz"}.get(m)
        if package:
            try:
                __import__(package)
            except ImportError:
                continue
        result.append(m)
    return result


def relayout_subsystems(dm: DraftMap, methods=METHODS, margin: float = 0.25, circle=(), **condense_args) -> list[dict]:
    """Step 3b: lay out each subsystem with added reactions with `methods`, keep the best layout (see module doc).

    The layouts are compared after condensing the subsystem as in step 4
    (`condense_args`). The subsystems in `circle` get the circular layout of
    their cycle (:func:`circle_layout`) if they have one. Changes the node
    positions in `dm` (before step 4) and the `layout` and `rigid` of the
    subsystems. Returns one row per subsystem: the crossings of each layout
    and the one kept.
    """
    methods = [m for m in available_methods(methods) if m != "current"]
    circle = set(circle) & set(dm.subsystems) if available_methods(["circle"]) else set()  # needs networkx
    keys = [k for k in sorted(dm.subsystems)
            if k in circle or any(n.kind == "enzyme" and n.status == "added" for n in dm.members(k))]
    bases = {k: subsystem_map(dm, k) for k in keys}
    # the Graphviz layouts of all subsystems, several processes at a time
    jobs = {(k, m): graphviz_spec(bases[k], **GRAPHVIZ[m]) for k in keys for m in methods if m in GRAPHVIZ}
    with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        computed = dict(zip(jobs, pool.map(graphviz_positions, [spec for spec, _, _ in jobs.values()])))
    rows = []
    for k in keys:
        base = bases[k]
        results = {}
        for m in ["current", *methods, *(["circle"] if k in circle else [])]:
            s = copy.deepcopy(base)
            if m in GRAPHVIZ:
                apply_graphviz(s, computed[k, m], *jobs[k, m][1:])
            elif m != "current" and LAYOUTS[m](s) is False:
                continue  # no cycle
            positions = {i: (n.x, n.y) for i, n in s.nodes.items()}
            rigid = copy.deepcopy(getattr(s.subsystems[k], "rigid", []))
            condense_subsystems(s, **condense_args)
            metrics = internal_metrics(s)
            results[m] = (metrics["crossings_main"], metrics["crossings_all"], s.subsystems[k].w * s.subsystems[k].h,
                          positions, rigid)
        rank = {m: r[:3] for m, r in results.items()}  # crossings (main, all), area
        best = min((m for m in results if m in KEEPING), key=rank.get)
        redraw = [m for m in results if m in REDRAWING]
        if redraw:
            other = min(redraw, key=rank.get)
            if results[other][0] < (1 - margin) * results[best][0]:
                best = other
        if "circle" in results:  # asked for
            best = "circle"
        for i, (x, y) in results[best][3].items():
            dm.nodes[i].x, dm.nodes[i].y = x, y
        dm.subsystems[k].layout, dm.subsystems[k].rigid = best, results[best][4]
        rows.append(dict(subsystem=k, layout=best, **{m: dict(crossings_main=r[0], crossings_all=r[1], area=r[2])
                                                       for m, r in results.items()}))
    chosen = {m: sum(r["layout"] == m for r in rows) for m in ["current", *methods, *(["circle"] if circle else [])]}
    before = sum(r["current"]["crossings_main"] for r in rows)
    after = sum(r[r["layout"]]["crossings_main"] for r in rows)
    circles = f" ({sum(r['layout'] == 'circle' for r in rows)} of them on a circle, --circle)" if circle else ""
    dm.log.append(
        f"step 3b: {len(rows)} subsystems with added reactions laid out again ({', '.join(methods) or 'no layouts'}): "
        f"{chosen}; crossings of main-metabolite links in them {before} -> {after}{circles}"
    )
    return rows
