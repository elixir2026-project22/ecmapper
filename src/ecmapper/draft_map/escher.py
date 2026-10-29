"""Export of a draft map as an Escher map (JSON, schema 1-0-0)."""

from __future__ import annotations

import json
import math
import re

from .core import PROTONS, DraftMap, ModelData, text_width

SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"


class _Ids:
    def __init__(self):
        self.n = 0

    def __call__(self) -> str:
        self.n += 1
        return str(self.n)


def to_escher(dm: DraftMap, md: ModelData, map_name: str, scale: float = 2.5, subsystem_labels: bool = True,
              protons: bool = False, description: str = "") -> list:
    """An Escher map of the model reactions of a draft map.

    Every model reaction of an enzyme node becomes an Escher reaction with its
    midmarker at the node (further reactions of the node below it) and a
    multimarker towards its substrates and one towards its products. Drawn
    metabolites are the metabolite nodes of the draft map: main metabolites
    (primary nodes) and cofactors (secondary nodes) whose names match the
    model metabolites; one Escher node per drawn node and model metabolite
    (compartment). Metabolites of a reaction that are not drawn get secondary
    nodes next to its multimarkers; protons are not drawn unless `protons`.
    Inter-pathway links and subsystem boxes are not exported; the names of
    the subsystems are text labels if `subsystem_labels`.
    """
    new_id = _Ids()
    nodes, reactions, text_labels = {}, {}, {}
    escher_node = {}  # (drawn node ID, BiGG ID) -> Escher node ID
    offsets = {}  # drawn node ID -> list of BiGG IDs (later ones are drawn below the first)
    edges = {}
    for e in dm.edges:
        edges.setdefault(e.enzyme, []).append(e)

    def metabolite_node(drawn, bigg_id, primary):
        key = (drawn.id, bigg_id)
        if key not in escher_node:
            ids = offsets.setdefault(drawn.id, [])
            ids.append(bigg_id)
            k = len(ids) - 1
            x, y = drawn.x * scale, (drawn.y + 34 * k) * scale
            escher_node[key] = node_id = new_id()
            nodes[node_id] = dict(
                node_type="metabolite", x=x, y=y, bigg_id=bigg_id, name=md.metabolites[bigg_id]["name"] or bigg_id,
                label_x=x + (15 if primary else 8) * scale / 2.5, label_y=y - (22 if primary else 12) * scale / 2.5,
                node_is_primary=primary,
            )
        return escher_node[key]

    def free_node(bigg_id, x, y):
        node_id = new_id()
        nodes[node_id] = dict(node_type="metabolite", x=x, y=y, bigg_id=bigg_id,
                              name=md.metabolites[bigg_id]["name"] or bigg_id, label_x=x + 8, label_y=y - 12,
                              node_is_primary=False)
        return node_id

    def assign(enzyme, met_ids):
        """Drawn metabolite nodes of the model metabolites of a reaction at an enzyme node (one to one).

        A metabolite matches a drawn node that was added for it, a main
        metabolite node anchored to it or a node with common keys (KEGG
        compounds, names); the best matches (most common keys, same role) are
        taken first, so that e.g. the two forms of a racemase get both nodes.
        """
        candidates = []
        for met_id in met_ids:
            base, keys, main = md.base[met_id], md.metabolites[met_id]["keys"], md.is_main(met_id)
            for e in edges.get(enzyme.id, []):
                node = dm.nodes[e.met]
                common = len(keys & md.node_met_keys.get((enzyme.id, node.id), frozenset()))
                if met_id in getattr(node, "model_met_ids", ()):
                    score = 1000
                elif main and e.kind == "main" and base in node.model_mets:
                    score = 100 + common
                elif common:
                    score = common + (5 if (e.kind == "main") == main else 0)
                else:
                    continue
                candidates.append((score, met_id, node))
        result, used = {}, set()
        for score, met_id, node in sorted(candidates, key=lambda c: -c[0]):
            if met_id not in result and node.id not in used:
                result[met_id] = node
                used.add(node.id)
        return result

    for enzyme in dm.nodes.values():
        if enzyme.kind != "enzyme" or not enzyme.reactions:
            continue
        for idx, rid in enumerate(enzyme.reactions):
            r = md.reactions[rid]
            mx, my = enzyme.x * scale, (enzyme.y + 34 * idx) * scale
            placed = {}  # metabolite -> Escher node
            missing = []
            mets = [m for m in r["metabolites"] if protons or md.base[m] not in PROTONS]
            drawn = assign(enzyme, mets)
            for met_id in mets:
                if met_id in drawn:
                    placed[met_id] = metabolite_node(drawn[met_id], met_id, drawn[met_id].kind == "main")
                else:
                    missing.append(met_id)
            # multimarkers towards the drawn substrates and products
            markers = {}
            for side, sign in (("sub", -1), ("prod", 1)):
                pts = [(nodes[n]["x"], nodes[n]["y"]) for m, n in placed.items() if r["metabolites"][m] * sign > 0]
                if pts:
                    cx = sum(p[0] for p in pts) / len(pts)
                    cy = sum(p[1] for p in pts) / len(pts)
                else:
                    cx, cy = mx + sign * 60 * scale, my
                d = math.dist((cx, cy), (mx, my)) or 1.0
                step = min(10 * scale, d / 3)
                markers[side] = (mx + (cx - mx) / d * step, my + (cy - my) / d * step)
            mid = new_id()
            nodes[mid] = dict(node_type="midmarker", x=mx, y=my)
            segments = {}
            marker_ids = {}
            for side, (x, y) in markers.items():
                marker_ids[side] = new_id()
                nodes[marker_ids[side]] = dict(node_type="multimarker", x=x, y=y)
                segments[new_id()] = dict(from_node_id=mid, to_node_id=marker_ids[side], b1=None, b2=None)
            # metabolites that are not drawn: secondary nodes around the multimarker of their side
            for side, sign in (("sub", -1), ("prod", 1)):
                ms = [m for m in missing if r["metabolites"][m] * sign > 0]
                x, y = markers[side]
                ax, ay = x - mx, y - my
                angle = math.atan2(ay, ax) if (ax or ay) else (0.0 if sign > 0 else math.pi)
                for k, m in enumerate(ms):
                    a = angle + (k - (len(ms) - 1) / 2) * 0.6
                    placed[m] = free_node(m, x + math.cos(a) * 40 * scale, y + math.sin(a) * 40 * scale)
            for m, node_id in placed.items():
                c = r["metabolites"][m]
                marker = marker_ids["prod" if c > 0 else "sub"]
                if c > 0:
                    segments[new_id()] = dict(from_node_id=marker, to_node_id=node_id, b1=None, b2=None)
                else:
                    segments[new_id()] = dict(from_node_id=node_id, to_node_id=marker, b1=None, b2=None)
            genes = [dict(bigg_id=g, name=name) for g, name in r.get("genes", [])]
            reactions[new_id()] = dict(
                name=r["name"] or rid, bigg_id=rid, reversibility=bool(r["reversible"]),
                label_x=mx + 12 * scale / 2.5, label_y=my + 30 * scale / 2.5,
                gene_reaction_rule=r.get("gpr", "") or "", genes=genes,
                metabolites=[dict(bigg_id=m, coefficient=c) for m, c in r["metabolites"].items()],
                segments=segments,
            )
    if subsystem_labels:
        for s in dm.subsystems.values():
            text_labels[new_id()] = dict(text=s.label, x=(s.x + s.w / 2 - text_width(s.label, 18) / 2) * scale,
                                         y=s.y * scale)
    xs = [n["x"] for n in nodes.values()] + [t["x"] for t in text_labels.values()]
    ys = [n["y"] for n in nodes.values()] + [t["y"] for t in text_labels.values()]
    margin = 200
    canvas = dict(x=min(xs) - margin, y=min(ys) - margin, width=max(xs) - min(xs) + 2 * margin,
                  height=max(ys) - min(ys) + 2 * margin)
    head = dict(map_name=map_name, map_id=re.sub(r"\W+", "_", map_name), map_description=description,
                homepage="https://escher.github.io", schema=SCHEMA)
    return [head, dict(reactions=reactions, nodes=nodes, text_labels=text_labels, canvas=canvas)]


def write_escher(dm: DraftMap, md: ModelData, path: str, map_name: str, **kwargs) -> dict:
    data = to_escher(dm, md, map_name, **kwargs)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    body = data[1]
    return dict(reactions=len(body["reactions"]), nodes=len(body["nodes"]), text_labels=len(body["text_labels"]),
                bigg_reactions=len({r["bigg_id"] for r in body["reactions"].values()}))
