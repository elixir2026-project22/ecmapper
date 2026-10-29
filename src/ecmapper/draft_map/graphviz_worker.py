"""Graphviz layout of one graph, run as a script in a process of its own (see :mod:`.relayout`).

Reads the graph as JSON from stdin ({"prog", "args", "directed", "graph", "node": attributes, "nodes": [[name,
attributes]], "edges": [[a, b]], "routes"}) and writes {name: [x, y]} (points, y upwards) to stdout; with "routes",
{"positions": {name: [x, y]}, "routes": [[[x, y], ...] for each edge]} (the control points of the edge splines).
Only needs pygraphviz.
"""

import json
import sys


def layout(spec: dict) -> dict:
    import pygraphviz as pgv

    g = pgv.AGraph(directed=spec["directed"], strict=False)
    g.graph_attr.update(spec["graph"])
    g.node_attr.update(spec["node"])
    for name, attrs in spec["nodes"]:
        g.add_node(name, **attrs)
    for i, (a, b) in enumerate(spec["edges"]):
        g.add_edge(a, b, key=str(i))
    g.layout(prog=spec["prog"], args=spec.get("args", ""))
    positions = {name: [float(v) for v in g.get_node(name).attr["pos"].rstrip("!").split(",")[:2]]
                 for name, _ in spec["nodes"]}
    if not spec.get("routes"):
        return positions
    routes = []
    for i, (a, b) in enumerate(spec["edges"]):
        pos = g.get_edge(a, b, key=str(i)).attr["pos"] or ""
        points = [p for p in pos.split() if not p.startswith(("e,", "s,"))]
        routes.append([[float(v) for v in p.split(",")[:2]] for p in points])
    return {"positions": positions, "routes": routes}


if __name__ == "__main__":
    json.dump(layout(json.load(sys.stdin)), sys.stdout)
