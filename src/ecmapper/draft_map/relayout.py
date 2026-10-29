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
their enzyme node. Graphviz layouts need pygraphviz (2.0 or later ships
Graphviz), networkx layouts networkx; missing ones are skipped.

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

from .core import DraftMap
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
    """Boxes of the main nodes; an enzyme node with its cofactors is one box."""
    owner = _owners(s)
    box = {i: [n.x - n.w / 2, n.y - n.h / 2, n.x + n.w / 2, n.y + n.h / 2]
           for i, n in s.nodes.items() if n.kind != "cofactor"}
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


LAYOUTS = {"spring": spring_layout, **{m: (lambda s, a=a: graphviz_layout(s, **a)) for m, a in GRAPHVIZ.items()}}


def available_methods(methods=METHODS) -> list[str]:
    """The methods whose packages are installed."""
    result = []
    for m in methods:
        package = {"spring": "networkx", "fdp": "pygraphviz", "dot": "pygraphviz", "neato": "pygraphviz"}.get(m)
        if package:
            try:
                __import__(package)
            except ImportError:
                continue
        result.append(m)
    return result


def relayout_subsystems(dm: DraftMap, methods=METHODS, margin: float = 0.25, **condense_args) -> list[dict]:
    """Step 3b: lay out each subsystem with added reactions with `methods`, keep the best layout (see module doc).

    The layouts are compared after condensing the subsystem as in step 4
    (`condense_args`). Changes the node positions in `dm` (before step 4).
    Returns one row per subsystem: the crossings of each layout and the one
    kept.
    """
    methods = [m for m in available_methods(methods) if m != "current"]
    keys = [k for k in sorted(dm.subsystems)
            if any(n.kind == "enzyme" and n.status == "added" for n in dm.members(k))]
    bases = {k: subsystem_map(dm, k) for k in keys}
    # the Graphviz layouts of all subsystems, several processes at a time
    jobs = {(k, m): graphviz_spec(bases[k], **GRAPHVIZ[m]) for k in keys for m in methods if m in GRAPHVIZ}
    with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        computed = dict(zip(jobs, pool.map(graphviz_positions, [spec for spec, _, _ in jobs.values()])))
    rows = []
    for k in keys:
        base = bases[k]
        results = {}
        for m in ["current", *methods]:
            s = copy.deepcopy(base)
            if m in GRAPHVIZ:
                apply_graphviz(s, computed[k, m], *jobs[k, m][1:])
            elif m != "current":
                LAYOUTS[m](s)
            positions = {i: (n.x, n.y) for i, n in s.nodes.items()}
            condense_subsystems(s, **condense_args)
            metrics = internal_metrics(s)
            results[m] = (metrics["crossings_main"], metrics["crossings_all"], s.subsystems[k].w * s.subsystems[k].h,
                          positions)
        rank = {m: r[:3] for m, r in results.items()}  # crossings (main, all), area
        best = min((m for m in results if m in KEEPING), key=rank.get)
        redraw = [m for m in results if m in REDRAWING]
        if redraw:
            other = min(redraw, key=rank.get)
            if results[other][0] < (1 - margin) * results[best][0]:
                best = other
        for i, (x, y) in results[best][3].items():
            dm.nodes[i].x, dm.nodes[i].y = x, y
        rows.append(dict(subsystem=k, layout=best, **{m: dict(crossings_main=r[0], crossings_all=r[1], area=r[2])
                                                       for m, r in results.items()}))
    chosen = {m: sum(r["layout"] == m for r in rows) for m in ["current", *methods]}
    before = sum(r["current"]["crossings_main"] for r in rows)
    after = sum(r[r["layout"]]["crossings_main"] for r in rows)
    dm.log.append(
        f"step 3b: {len(rows)} subsystems with added reactions laid out again ({', '.join(methods) or 'no layouts'}): "
        f"{chosen}; crossings of main-metabolite links in them {before} -> {after}"
    )
    return rows
