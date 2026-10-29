"""Graphviz layout of one graph, run as a script in a process of its own (see :mod:`.relayout`).

Reads the graph as JSON from stdin ({"prog", "directed", "graph", "node": attributes, "nodes": [[name, attributes]],
"edges": [[a, b]]}) and writes {name: [x, y]} (points, y upwards) to stdout. Only needs pygraphviz.
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
    for a, b in spec["edges"]:
        g.add_edge(a, b)
    g.layout(prog=spec["prog"])
    return {name: [float(v) for v in g.get_node(name).attr["pos"].rstrip("!").split(",")[:2]]
            for name, _ in spec["nodes"]}


if __name__ == "__main__":
    json.dump(layout(json.load(sys.stdin)), sys.stdout)
