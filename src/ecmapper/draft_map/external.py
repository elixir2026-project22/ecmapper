"""Arrangements of the subsystem boxes with general graph layout tools (Graphviz, ELK) for comparison."""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np

from .layout import Ports, compact_boxes

ELK_SCRIPT = Path(__file__).with_name("elk_layout.js")


def _compact(x, y, ports, gap):
    """Top left positions from centres after an order preserving compaction (removes overlaps and gaps)."""
    cx, cy = compact_boxes(x, y, ports.w + gap, ports.h + gap, 0.0, 0.0, 0.5, rounds=3, max_scale=10.0)
    return cx - ports.w / 2, cy - ports.h / 2


def graphviz_arrangement(dm, ports: Ports, prog="neato", gap=80.0):
    """Graphviz (neato: stress majorisation, fdp: force directed) on the box graph, overlaps removed with VPSC."""
    import pygraphviz as pgv

    g = pgv.AGraph(strict=False, directed=False)
    g.graph_attr.update(overlap="vpsc", sep=f"+{gap / 2:.0f}", splines="false", outputorder="edgesfirst")
    for i, name in enumerate(ports.names):
        cx, cy = dm.orig_centres[name]
        g.add_node(i, shape="box", fixedsize="true", width=ports.w[i] / 72, height=ports.h[i] / 72,
                   pos=f"{cx:.0f},{-cy:.0f}", label="")
    weights = Counter(tuple(sorted((int(a), int(b)))) for a, b in zip(ports.a, ports.b))
    diag = np.hypot(ports.w, ports.h) / 2 / 72
    for (a, b), k in weights.items():
        g.add_edge(a, b, weight=k, len=float(diag[a] + diag[b]) + gap / 72)
    g.layout(prog=prog)
    x = np.array([float(g.get_node(i).attr["pos"].split(",")[0]) for i in range(len(ports.names))])
    y = np.array([-float(g.get_node(i).attr["pos"].split(",")[1]) for i in range(len(ports.names))])
    return _compact(x, y, ports, gap)


def elk_arrangement(dm, ports: Ports, algorithm="layered", gap=80.0, aspect=1.6):
    """ELK (layered: Sugiyama with crossing minimisation and fixed port positions; stress).

    Needs ``node`` and the npm package ``elkjs`` (``$NODE_PATH`` set to the ``node_modules`` directory with it).
    """
    children = []
    for i, name in enumerate(ports.names):
        children.append(dict(id=f"b{i}", width=float(ports.w[i]), height=float(ports.h[i]), ports=[],
                             layoutOptions={"elk.portConstraints": "FIXED_POS"}))
    edges = []
    for k, (a, b) in enumerate(zip(ports.a, ports.b)):
        pa, pb = f"b{a}_p{k}a", f"b{b}_p{k}b"
        if algorithm != "layered":
            edges.append(dict(id=f"e{k}", sources=[f"b{a}"], targets=[f"b{b}"]))
            continue
        children[a]["ports"].append(dict(id=pa, x=float(ports.pa[k, 0]), y=float(ports.pa[k, 1]), width=1, height=1))
        children[b]["ports"].append(dict(id=pb, x=float(ports.pb[k, 0]), y=float(ports.pb[k, 1]), width=1, height=1))
        edges.append(dict(id=f"e{k}", sources=[pa], targets=[pb]))
    options = {
        "elk.algorithm": algorithm,
        "elk.spacing.nodeNode": gap,
        "elk.aspectRatio": aspect,
        "elk.separateConnectedComponents": "true",
        "elk.spacing.componentComponent": gap,
    }
    if algorithm == "layered":
        options.update({
            "elk.direction": "RIGHT",
            "elk.layered.spacing.nodeNodeBetweenLayers": gap,
            "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
            "elk.layered.wrapping.strategy": "MULTI_EDGE",
            "elk.edgeRouting": "POLYLINE",
        })
    graph = dict(id="root", layoutOptions=options, children=children, edges=edges)
    result = subprocess.run(["node", str(ELK_SCRIPT)], input=json.dumps(graph), capture_output=True, text=True,
                            check=False)
    if result.returncode:
        raise RuntimeError(f"ELK failed: {result.stderr[-2000:]}")
    out = json.loads(result.stdout)
    pos = {c["id"]: (c["x"], c["y"]) for c in out["children"]}
    x = np.array([pos[f"b{i}"][0] for i in range(len(ports.names))]) + ports.w / 2
    y = np.array([pos[f"b{i}"][1] for i in range(len(ports.names))]) + ports.h / 2
    # layered drawings leave much space between the layers: compact them like the others
    return _compact(x, y, ports, gap)
