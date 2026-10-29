"""Edge routes around nodes: main edges that pass through other nodes as Graphviz splines.

A straight main edge can pass through other nodes of its subsystem (and look
as if it joined them), e.g. after a Graphviz ``dot`` layout (step 3b), whose
own routes go around the nodes of the layers between, and condensing (step
4). Graphviz ``neato -n2`` keeps the final node positions and routes the
edges as splines around the node boxes (``splines=true``); only edges that
pass through another node take their route, the others stay straight. The
routes are cubic Bezier control points from the enzyme node to the
metabolite node, in map coordinates; the SVG draws them as paths and the
Escher export follows them (:mod:`.escher`).
"""

from __future__ import annotations

import itertools
import math
from concurrent.futures import ThreadPoolExecutor

from .core import DraftMap
from .relayout import graphviz_positions


def route_spec(dm: DraftMap, k: str) -> tuple[dict, list]:
    """The Graphviz graph of subsystem `k`: all its nodes fixed (boxes), its main edges to route."""
    members = dm.members(k)
    ids = {n.id for n in members}
    edges = [e for e in dm.edges if e.kind == "main" and e.enzyme in ids and e.met in ids]
    nodes = [[n.id, dict(pos=f"{n.x:.2f},{-n.y:.2f}!", width=f"{n.w / 72:.4f}", height=f"{n.h / 72:.4f}")]
             for n in members]
    spec = dict(prog="neato", args="-n2", directed=False, routes=True,
                graph=dict(splines="true", esep="+4", inputscale="72"),
                node=dict(shape="box", fixedsize="true", label=""), nodes=nodes,
                edges=[[e.enzyme, e.met] for e in edges])
    return spec, edges


def _through(dm: DraftMap, e, boxes: dict, margin: float = 2.0) -> bool:
    """Whether the straight edge `e` passes through another node box of its subsystem."""
    a, b = dm.nodes[e.enzyme], dm.nodes[e.met]
    points = [(a.x + (b.x - a.x) * i / 60, a.y + (b.y - a.y) * i / 60) for i in range(61)]
    return any(x0 + margin <= x <= x1 - margin and y0 + margin <= y <= y1 - margin
               for i, (x0, y0, x1, y1) in boxes.items() if i not in (e.enzyme, e.met) for x, y in points)


def route_edges(dm: DraftMap, subsystems=None) -> int:
    """Route the main edges that pass through other nodes around them (Graphviz, several processes at a time).

    For the subsystems in `subsystems` (default all), sets `route` of each
    main edge: if the straight edge passes through another node of the
    subsystem, the control points of its route (map coordinates, from the
    enzyme node to the metabolite node), else None. Returns the number of
    routed edges.
    """
    subsystems = list(dm.subsystems) if subsystems is None else subsystems
    for e in dm.edges:
        e.route = None
    specs = {k: route_spec(dm, k) for k in subsystems if dm.members(k)}
    specs = {k: (spec, edges) for k, (spec, edges) in specs.items() if edges}
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = dict(zip(specs, pool.map(graphviz_positions, [spec for spec, _ in specs.values()])))
    bent = 0
    for k, (spec, edges) in specs.items():
        result = results[k]
        # Graphviz moves the drawing (not the scale): back to map coordinates, y downwards
        name, attrs = spec["nodes"][0]
        x_in, y_in = (float(v) for v in attrs["pos"].rstrip("!").split(","))
        x_out, y_out = result["positions"][name]
        dx, dy = x_in - x_out, y_in - y_out
        boxes = {n.id: n.box for n in dm.members(k)}
        for e, route in zip(edges, result["routes"]):
            if not _through(dm, e, boxes):
                continue
            points = [(x + dx, -(y + dy)) for x, y in route]
            if len(points) >= 4:
                enzyme = dm.nodes[e.enzyme]
                if math.dist(points[-1], (enzyme.x, enzyme.y)) < math.dist(points[0], (enzyme.x, enzyme.y)):
                    points.reverse()  # from the enzyme node to the metabolite node
            e.route = points if len(points) >= 4 else None
            bent += e.route is not None
    return bent


def bezier_path(points: list[tuple[float, float]]) -> str:
    """SVG path of cubic Bezier control points (3n + 1 points)."""
    out = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    for i in range(1, len(points) - 2, 3):
        (x1, y1), (x2, y2), (x3, y3) = points[i:i + 3]
        out.append(f"C {x1:.1f} {y1:.1f} {x2:.1f} {y2:.1f} {x3:.1f} {y3:.1f}")
    return " ".join(out)


def sample(points: list[tuple[float, float]], n: int = 12) -> list[tuple[float, float]]:
    """Points along a chain of cubic Bezier curves, `n` per curve."""
    out = [points[0]]
    for i in range(0, len(points) - 3, 3):
        p0, p1, p2, p3 = points[i:i + 4]
        for j in range(1, n + 1):
            t = j / n
            c = ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3)
            out.append((c[0] * p0[0] + c[1] * p1[0] + c[2] * p2[0] + c[3] * p3[0],
                        c[0] * p0[1] + c[1] * p1[1] + c[2] * p2[1] + c[3] * p3[1]))
    return out


def fit_cubic(p0, p3, points) -> tuple[tuple[float, float], tuple[float, float]]:
    """Control points b1, b2 of the cubic Bezier curve from p0 to p3 closest to `points` (least squares)."""
    pts = [p0, *points, p3]
    lengths = [0.0]
    for a, b in itertools.pairwise(pts):
        lengths.append(lengths[-1] + math.dist(a, b))
    total = lengths[-1] or 1.0
    s11 = s12 = s22 = 0.0
    r1, r2 = [0.0, 0.0], [0.0, 0.0]
    for (x, y), length in zip(pts, lengths):
        t = length / total
        a1, a2 = 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2
        rx = x - (1 - t) ** 3 * p0[0] - t ** 3 * p3[0]
        ry = y - (1 - t) ** 3 * p0[1] - t ** 3 * p3[1]
        s11, s12, s22 = s11 + a1 * a1, s12 + a1 * a2, s22 + a2 * a2
        r1 = [r1[0] + a1 * rx, r1[1] + a1 * ry]
        r2 = [r2[0] + a2 * rx, r2[1] + a2 * ry]
    det = s11 * s22 - s12 * s12
    if abs(det) < 1e-9:
        return ((2 * p0[0] + p3[0]) / 3, (2 * p0[1] + p3[1]) / 3), ((p0[0] + 2 * p3[0]) / 3, (p0[1] + 2 * p3[1]) / 3)
    b1 = tuple((s22 * r1[i] - s12 * r2[i]) / det for i in range(2))
    b2 = tuple((s11 * r2[i] - s12 * r1[i]) / det for i in range(2))
    return b1, b2
