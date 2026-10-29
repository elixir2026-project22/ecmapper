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
              protons: bool = False, description: str = "", curved: bool = True) -> list:
    """An Escher map of the model reactions of a draft map.

    Every model reaction of an enzyme node becomes an Escher reaction with its
    midmarker at the node and a multimarker towards its substrates and one
    towards its products. Drawn metabolites are the metabolite nodes of the
    draft map: main metabolites (primary nodes) and cofactors (secondary
    nodes) whose names match the model metabolites; one Escher node per drawn
    node and model metabolite (compartment). Metabolites of a reaction that
    are not drawn get secondary nodes next to its multimarkers; protons are
    not drawn unless `protons`. Inter-pathway links and subsystem boxes are
    not exported; the names of the subsystems are text labels if
    `subsystem_labels`. Metabolite nodes that continue on other maps
    (`stubs`, see :mod:`.split`) get a text label "→ map names" next to them
    (:func:`place_text_labels`).

    With `curved`, the reactions get straight axes and curved branches
    (:func:`curve_reactions`): the reactions of an enzyme node side by side
    and the metabolites that are not drawn next to the reaction axis.
    Otherwise the segments are straight and further reactions of an enzyme
    node are drawn below the first.
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

    movable = set()

    def free_node(bigg_id, x, y):
        node_id = new_id()
        movable.add(node_id)
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
            mx, my = enzyme.x * scale, (enzyme.y + (0 if curved else 34 * idx)) * scale
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
    references = {}  # text label ID -> the extent of the metabolite nodes it belongs to (x, top y, bottom y)
    for drawn_id, ids in offsets.items():
        drawn = dm.nodes[drawn_id]
        if getattr(drawn, "stubs", None):
            label_id = new_id()
            text_labels[label_id] = dict(text="→ " + ", ".join(drawn.stubs), x=drawn.x * scale,
                                         y=(drawn.y + 34 * (len(ids) - 1)) * scale + 66)
            references[label_id] = (drawn.x * scale, drawn.y * scale, (drawn.y + 34 * (len(ids) - 1)) * scale)
    if subsystem_labels:
        for s in dm.subsystems.values():
            text_labels[new_id()] = dict(text=s.label, x=(s.x + s.w / 2 - text_width(s.label, 18) / 2) * scale,
                                         y=s.y * scale)
    body = dict(reactions=reactions, nodes=nodes, text_labels=text_labels)
    if curved:
        curve_reactions(body, movable)
    place_text_labels(body, references)
    xs = [n["x"] for n in nodes.values()] + [t["x"] for t in text_labels.values()]
    xs += [t["x"] + ESCHER_FONTS["text"][1] * len(t["text"]) for t in text_labels.values()]
    ys = [n["y"] for n in nodes.values()] + [t["y"] for t in text_labels.values()]
    margin = 200
    canvas = dict(x=min(xs) - margin, y=min(ys) - margin, width=max(xs) - min(xs) + 2 * margin,
                  height=max(ys) - min(ys) + 2 * margin)
    head = dict(map_name=map_name, map_id=re.sub(r"\W+", "_", map_name), map_description=description,
                homepage="https://escher.github.io", schema=SCHEMA)
    return [head, dict(body, canvas=canvas)]


def _unit(x: float, y: float) -> tuple[float, float] | None:
    d = math.hypot(x, y)
    return (x / d, y / d) if d > 1e-9 else None


def _reaction_parts(nodes: dict, reaction: dict) -> tuple[str, dict[str, int], dict[str, list]] | None:
    """The midmarker of a reaction, its multimarkers with their side (-1 substrates, 1 products) and their segments.

    None if the reaction is not drawn as Escher's builder draws it (a
    midmarker, a multimarker for each side, segments from the multimarkers
    to the metabolites).
    """
    coefficient = {m["bigg_id"]: m["coefficient"] for m in reaction["metabolites"]}
    mid, markers = None, {}
    for sid, s in reaction["segments"].items():
        a, b = s["from_node_id"], s["to_node_id"]
        types = nodes[a]["node_type"], nodes[b]["node_type"]
        if types in (("midmarker", "multimarker"), ("multimarker", "midmarker")):
            m = a if types[1] == "multimarker" else b
            if mid not in (None, m):
                return None
            mid = m
            markers.setdefault(b if m == a else a, [])
        elif types in (("multimarker", "metabolite"), ("metabolite", "multimarker")):
            marker, met = (a, b) if types[0] == "multimarker" else (b, a)
            markers.setdefault(marker, []).append((sid, met))
        else:
            return None
    if mid is None or len(markers) != 2:
        return None
    sides = {}
    for marker, segments in markers.items():
        total = 0
        for _, met in segments:
            c = coefficient.get(nodes[met]["bigg_id"], 0)
            total += (c > 0) - (c < 0)
        sides[marker] = (total > 0) - (total < 0)
    (a, side_a), (b, side_b) = sides.items()
    if not markers[a] and side_b:  # a side without metabolites (drawn) is the other side
        sides[a] = -side_b
    elif not markers[b] and side_a:
        sides[b] = -side_a
    if sorted(sides.values()) != [-1, 1]:
        return None
    return mid, sides, markers


def curve_reactions(body: dict, movable=frozenset(), lane: float = 25.0, trunk: float = 0.3,
                    marker_distance: tuple[float, float] = (20.0, 80.0), max_lead: float = 100.0) -> dict:
    """Draw the reactions of an Escher map as Escher's builder does: a straight axis and curved branches.

    Post-processing of an Escher map (its second part, `body`); similar to
    BioEMMA (multimarkers on a straight reaction axis) and to the new
    reactions of Escher's builder (curved branches):

    - The axis of a reaction runs through its midmarker from its substrates
      to its products (the primary metabolites, all metabolites of a side
      if none is primary): along the difference of the directions from the
      midmarker to the centres of the products and of the substrates.
    - Reactions with the same midmarker position (the model reactions of an
      enzyme node) are drawn side by side in lanes `lane` apart,
      perpendicular to their common axis; their labels are stacked next to
      the lanes.
    - The multimarkers are on the axis, `trunk` of the way to the centre
      of their metabolites, within `marker_distance`.
    - The segments between a multimarker and its metabolites are cubic
      Bezier curves (b1, b2): they leave the multimarker along the axis
      and bend towards the metabolite; metabolites on the axis get straight
      segments.
    - Metabolite nodes in `movable` (secondary metabolites that are not
      drawn on the source map) are put next to the multimarker of their
      side, alternately left and right of the axis.

    Metabolite nodes in general are not moved. Reactions drawn in other
    ways are left as they are. Returns counts of the curved reactions and
    segments.
    """
    nodes = body["nodes"]
    parts = {rid: p for rid, r in body["reactions"].items() if (p := _reaction_parts(nodes, r))}

    def axis(mid, sides, markers):
        """The direction of the reaction (from the substrates to the products)."""
        directions = {}
        for side in (-1, 1):
            mets = [met for marker, s in sides.items() if s == side for _, met in markers[marker]]
            mets = [m for m in mets if nodes[m].get("node_is_primary")] or mets
            if not mets:
                directions[side] = None
                continue
            cx = sum(nodes[m]["x"] for m in mets) / len(mets) - nodes[mid]["x"]
            cy = sum(nodes[m]["y"] for m in mets) / len(mets) - nodes[mid]["y"]
            directions[side] = _unit(cx, cy)
        sub, prod = directions[-1], directions[1]
        if sub and prod:
            return _unit(prod[0] - sub[0], prod[1] - sub[1]) or (-prod[1], prod[0])
        if prod:
            return prod
        if sub:
            return -sub[0], -sub[1]
        return 0.0, 1.0

    groups = {}
    for rid, (mid, _, _) in parts.items():
        groups.setdefault((round(nodes[mid]["x"], 1), round(nodes[mid]["y"], 1)), []).append(rid)
    segments = 0
    for (x0, y0), rids in groups.items():
        axes = {rid: axis(*parts[rid]) for rid in rids}
        ux, uy = axes[rids[0]]
        sx = sum(a[0] * math.copysign(1, a[0] * ux + a[1] * uy) for a in axes.values())
        sy = sum(a[1] * math.copysign(1, a[0] * ux + a[1] * uy) for a in axes.values())
        ux, uy = _unit(sx, sy) or (ux, uy)
        nx, ny = -uy, ux  # normal, to the right or downwards
        if nx < -1e-9 or (abs(nx) <= 1e-9 and ny < 0):
            nx, ny = -nx, -ny
        k = len(rids)
        if k > 1:
            # lanes in the order of the metabolites they do not share (fewer crossings)
            mets = {rid: {met for segs in parts[rid][2].values() for _, met in segs} for rid in rids}
            shared = set.intersection(*mets.values())
            lateral = {}
            for rid in rids:
                own = mets[rid] - shared
                lateral[rid] = sum((nodes[m]["x"] - x0) * nx + (nodes[m]["y"] - y0) * ny for m in own) / max(len(own), 1)
            rids = sorted(rids, key=lateral.get)
        for i, rid in enumerate(rids):
            a = axes[rid]
            mid, sides, markers = parts[rid]
            sign = math.copysign(1, a[0] * ux + a[1] * uy)
            offset = (i - (k - 1) / 2) * lane
            mx, my = x0 + nx * offset, y0 + ny * offset
            nodes[mid]["x"], nodes[mid]["y"] = mx, my
            for marker, side in sides.items():
                wx, wy = side * sign * ux, side * sign * uy  # outwards along the axis
                mets = [met for _, met in markers[marker]]
                free = [m for m in mets if m in movable]
                fixed = [m for m in mets if m not in movable] or mets
                along = sum((nodes[m]["x"] - mx) * wx + (nodes[m]["y"] - my) * wy for m in fixed) / max(len(fixed), 1)
                d = min(max(trunk * along, marker_distance[0]), marker_distance[1])
                px, py = mx + wx * d, my + wy * d
                nodes[marker]["x"], nodes[marker]["y"] = px, py
                outer = 1 if offset > 0 else -1 if offset < 0 else 0
                for j, m in enumerate(free):
                    lateral = outer or (1 if j % 2 == 0 else -1)
                    step = j if outer else j // 2
                    dist = 45.0 + 30.0 * step
                    reach = abs(offset) + 55.0 + (25.0 * step if outer else 0.0)
                    nodes[m]["x"] = px + wx * dist + nx * lateral * reach
                    nodes[m]["y"] = py + wy * dist + ny * lateral * reach
                    nodes[m]["label_x"], nodes[m]["label_y"] = nodes[m]["x"] + 8, nodes[m]["y"] - 12
                for sid, met in markers[marker]:
                    vx, vy = nodes[met]["x"] - px, nodes[met]["y"] - py
                    length = math.hypot(vx, vy) or 1.0
                    lead = min(0.4 * length * max((vx * wx + vy * wy) / length, 0.25), max_lead)
                    b_marker = dict(x=px + wx * lead, y=py + wy * lead)
                    b_met = dict(x=nodes[met]["x"] + (b_marker["x"] - nodes[met]["x"]) * 0.35,
                                 y=nodes[met]["y"] + (b_marker["y"] - nodes[met]["y"]) * 0.35)
                    s = body["reactions"][rid]["segments"][sid]
                    s["b1"], s["b2"] = (b_marker, b_met) if s["from_node_id"] == marker else (b_met, b_marker)
                    segments += 1
        # labels: stacked right of the lanes (steep axes) or below them
        half = (k - 1) / 2 * lane
        for i, rid in enumerate(rids):
            r = body["reactions"][rid]
            width = 0.55 * 30 * len(r["bigg_id"])
            if abs(uy) >= abs(ux):
                r["label_x"] = x0 + abs(nx) * half + 15 + 30 * abs(ux / uy)
                r["label_y"] = y0 + 10 + (i - (k - 1) / 2) * 32
            else:
                r["label_x"] = x0 - width / 2
                r["label_y"] = y0 + abs(ny) * half + 15 + 30 * abs(uy / ux) + 22 + i * 32
    return dict(curved_reactions=len(parts), curved_segments=segments)


#: approximate label sizes in Escher (its CSS): font size, width per character
ESCHER_FONTS = {"metabolite": (20, 11.5), "reaction": (30, 17.0), "text": (50, 28.0)}


def _label_box(kind: str, text: str, x: float, y: float) -> tuple[float, float, float, float]:
    """The box of a label at (x, y) (start of the baseline, as Escher places labels)."""
    size, per_char = ESCHER_FONTS[kind]
    return x, y - 0.8 * size, x + per_char * len(text), y + 0.25 * size


def _segment_points(nodes: dict, segment: dict, step: float = 40.0) -> list[tuple[float, float]]:
    """Points about `step` apart along a segment (a cubic Bezier curve if it has control points)."""
    a, b = nodes[segment["from_node_id"]], nodes[segment["to_node_id"]]
    p0, p3 = (a["x"], a["y"]), (b["x"], b["y"])
    b1, b2 = segment.get("b1"), segment.get("b2")
    p1, p2 = ((b1["x"], b1["y"]), (b2["x"], b2["y"])) if b1 and b2 else (p0, p3)
    n = max(4, int((math.dist(p0, p1) + math.dist(p1, p2) + math.dist(p2, p3)) / step))
    points = []
    for i in range(n + 1):
        t = i / n
        c = ((1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3)
        points.append((c[0] * p0[0] + c[1] * p1[0] + c[2] * p2[0] + c[3] * p3[0],
                       c[0] * p0[1] + c[1] * p1[1] + c[2] * p2[1] + c[3] * p3[1]))
    return points


def place_text_labels(body: dict, anchors: dict[str, tuple[float, float, float]], cell: float = 100.0) -> None:
    """Put text labels next to their metabolite nodes where they cover the least.

    `anchors`: {text label ID: (x, top y, bottom y) of its metabolite nodes}.
    Each label goes right of, below, left of or above its nodes, whichever
    position covers the fewest points of segments, nodes and labels (other
    text labels count most); ties in this order.
    """
    nodes = body["nodes"]
    grid = {}  # (cell x, cell y) -> [(x, y, weight)]

    def add(x, y, weight):
        grid.setdefault((int(x // cell), int(y // cell)), []).append((x, y, weight))

    def add_box(box, weight, step=25.0):
        x0, y0, x1, y1 = box
        nx, ny = max(1, int((x1 - x0) / step)), max(1, int((y1 - y0) / step))
        for i in range(nx + 1):
            for j in range(ny + 1):
                add(x0 + (x1 - x0) * i / nx, y0 + (y1 - y0) * j / ny, weight)

    for r in body["reactions"].values():
        for segment in r["segments"].values():
            for x, y in _segment_points(nodes, segment):
                add(x, y, 1.0)
        if r.get("label_x") is not None:
            add_box(_label_box("reaction", r["bigg_id"], r["label_x"], r["label_y"]), 1.0)
    for n in nodes.values():
        if n["node_type"] == "metabolite":
            add(n["x"], n["y"], 4.0)
            add_box(_label_box("metabolite", n["bigg_id"], n["label_x"], n["label_y"]), 1.0)
    for label_id, t in body["text_labels"].items():
        if label_id not in anchors:
            add_box(_label_box("text", t["text"], t["x"], t["y"]), 10.0)

    def cost(box):
        x0, y0, x1, y1 = box
        total = 0.0
        for i in range(int(x0 // cell), int(x1 // cell) + 1):
            for j in range(int(y0 // cell), int(y1 // cell) + 1):
                total += sum(w for x, y, w in grid.get((i, j), ()) if x0 <= x <= x1 and y0 <= y <= y1)
        return total

    for label_id, (x, top, bottom) in anchors.items():
        t = body["text_labels"][label_id]
        width = ESCHER_FONTS["text"][1] * len(t["text"])
        candidates = [(x + 30, (top + bottom) / 2 + 22),  # right
                      (x - 20, bottom + 66),  # below
                      (x - 35 - width, (top + bottom) / 2 + 22),  # left
                      (x - 20, top - 50)]  # above
        best = min(candidates, key=lambda c: cost(_label_box("text", t["text"], *c)))
        t["x"], t["y"] = best
        add_box(_label_box("text", t["text"], *best), 10.0)


def write_escher(dm: DraftMap, md: ModelData, path: str, map_name: str, **kwargs) -> dict:
    data = to_escher(dm, md, map_name, **kwargs)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    body = data[1]
    return dict(reactions=len(body["reactions"]), nodes=len(body["nodes"]), text_labels=len(body["text_labels"]),
                bigg_reactions=len({r["bigg_id"] for r in body["reactions"].values()}))
