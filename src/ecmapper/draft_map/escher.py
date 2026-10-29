"""Export of a draft map as an Escher map (JSON, schema 1-0-0)."""

from __future__ import annotations

import json
import math
import re

from .core import LANE, PROTONS, DraftMap, ModelData, text_width

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

    The reactions get straight axes and curved branches
    (:func:`curve_reactions`): the reactions of an enzyme node side by side
    (`LANE` apart, as the layout leaves room for them), U-shaped reactions
    moved towards their metabolites where there is room, and the
    metabolites that are not drawn next to the reaction axis.
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
            mx, my = enzyme.x * scale, enzyme.y * scale
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
    curve_reactions(body, movable, lane=LANE * scale)
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


class _Points:
    """Weighted points with owners in a grid: what is near a place, leaving out some owners."""

    def __init__(self, cell: float = 100.0):
        self.cell, self.grid, self.owned = cell, {}, {}

    def add(self, owner: str, x: float, y: float, weight: float = 1.0) -> None:
        key, item = (int(x // self.cell), int(y // self.cell)), (x, y, weight, owner)
        self.grid.setdefault(key, []).append(item)
        self.owned.setdefault(owner, []).append((key, item))

    def remove(self, owner: str) -> None:
        for key, item in self.owned.pop(owner, []):
            self.grid[key].remove(item)

    def near(self, x: float, y: float, r: float, skip=frozenset()) -> float:
        total = 0.0
        for i in range(int((x - r) // self.cell), int((x + r) // self.cell) + 1):
            for j in range(int((y - r) // self.cell), int((y + r) // self.cell) + 1):
                for px, py, weight, owner in self.grid.get((i, j), ()):
                    if owner not in skip and (px - x) ** 2 + (py - y) ** 2 <= r * r:
                        total += weight
        return total


def _line(p, q, step: float = 30.0, skip: float = 0.0) -> list[tuple[float, float]]:
    """Points about `step` apart on the line from p to q, leaving out `skip` at both ends."""
    length = math.dist(p, q)
    n = max(1, int(length / step))
    return [(p[0] + (q[0] - p[0]) * i / n, p[1] + (q[1] - p[1]) * i / n) for i in range(n + 1)
            if skip <= length * i / n <= length - skip]


def _axis(nodes: dict, mid: str, sides: dict, markers: dict) -> tuple[float, float]:
    """The direction of a reaction (from its substrates to its products): along the difference of the directions
    from the midmarker to the centres of its products and of its substrates (primary ones if there are)."""
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


def _normal(ux: float, uy: float) -> tuple[float, float]:
    """The normal of a direction that points to the right or downwards."""
    nx, ny = -uy, ux
    return (-nx, -ny) if nx < -1e-9 or (abs(nx) <= 1e-9 and ny < 0) else (nx, ny)


def _lane_labels(x: float, y: float, ux: float, uy: float, ids: list[str], lane: float,
                 side: int = 1) -> list[tuple[float, float]]:
    """Label positions of reactions drawn side by side at (x, y) along (ux, uy).

    Stacked right of the lanes (steep axes) or below them (flat axes), or
    left of or above them (`side` -1).
    """
    k = len(ids)
    gap = (k - 1) / 2 * lane + 15
    out = []
    for i, text in enumerate(ids):
        width = ESCHER_FONTS["reaction"][1] * len(text)
        if abs(uy) >= abs(ux):
            dx = gap + 30 * abs(ux / uy)
            out.append((x + dx if side > 0 else x - dx - width, y + 10 + (i - (k - 1) / 2) * 32))
        else:
            dy = gap + 30 * abs(uy / ux)
            out.append((x - width / 2, y + dy + 22 + i * 32 if side > 0 else y - dy - 8 - (k - 1 - i) * 32))
    return out


def _footprint(x: float, y: float, ux: float, uy: float, ids: list[str], lane: float, side: int,
               reach: float) -> list[tuple[float, float]]:
    """Points of the area that reactions drawn side by side at (x, y) along (ux, uy) take.

    Their lanes, multimarkers and the secondary metabolites beside them (on
    `side` of the normal), up to `reach` along the axis, and their labels.
    """
    nx, ny = _normal(ux, uy)
    half = (len(ids) - 1) / 2 * lane
    near, outside = -(half + 20.0), half + 70.0
    n = max(1, int((outside - near) / 30))
    across = [side * (near + (outside - near) * i / n) for i in range(n + 1)]
    points = [(x + ux * a + nx * b, y + uy * a + ny * b) for a in range(-int(reach), int(reach) + 1, 30) for b in across]
    for (lx, ly), text in zip(_lane_labels(x, y, ux, uy, ids, lane, side), ids):
        x0, y0, x1, y1 = _label_box("reaction", text, lx, ly)
        points += [(x0 + (x1 - x0) * i / 4, (y0 + y1) / 2) for i in range(5)]
    return points


def _pull_reactions(body: dict, parts: dict, groups: list, lane: float, movable: set, offset: float = 90.0,
                    reach: float = 150.0, margin: float = 15.0) -> int:
    """Move U-shaped reactions towards their metabolites, into free space.

    A group of reactions drawn at one place (an enzyme node) is U-shaped if
    its substrates and products lie on the same side of it (less than 100°
    apart, seen from it). It moves towards the line between the centres of
    its substrates and of its products, to `offset` from that line plus half
    the width of its lanes, or 80%, 60%, ... of the way there: the farthest
    place where its footprint (its lanes with multimarkers and secondary
    metabolites up to `reach` along the line, and its labels) holds no node
    of other reactions, no more of their lines than where it is, and where
    its lines pass no more nodes. The secondary metabolites that only a
    moved group uses are added to `movable` (to be put next to its
    multimarkers, on the side away from the line, like its labels). Changes
    `groups` ([x, y, reaction IDs, direction moved away from: None]) and the
    midmarkers; returns the number of moved reactions.
    """
    nodes, reactions = body["nodes"], body["reactions"]
    ends = {rid: [met for segs in markers.values() for _, met in segs] for rid, (_, _, markers) in parts.items()}
    users = {}  # metabolite node -> reactions using it
    for rid, mets in ends.items():
        for met in mets:
            users.setdefault(met, set()).add(rid)
    obstacles, lines = _Points(), _Points()  # nodes; lines of reactions from the midmarker to their metabolites

    def add_lines(rid, x, y):
        for met in ends[rid]:
            for p in _line((x, y), (nodes[met]["x"], nodes[met]["y"]), skip=25.0):
                lines.add(rid, *p)

    for x, y, rids, _ in groups:
        for rid in rids:
            obstacles.add(rid, x, y)
            add_lines(rid, x, y)
        # their labels (where they go if the group stays)
        ids = [reactions[rid]["bigg_id"] for rid in rids]
        for (lx, ly), text in zip(_lane_labels(x, y, *_axis(nodes, *parts[rids[0]]), ids, lane), ids):
            x0_, y0_, x1_, y1_ = _label_box("reaction", text, lx, ly)
            for i in range(5):
                obstacles.add(rids[0], x0_ + (x1_ - x0_) * i / 4, (y0_ + y1_) / 2)
    for i, n in nodes.items():
        if n["node_type"] == "metabolite":
            obstacles.add(i, n["x"], n["y"])
    candidates = []
    for g in groups:
        x0, y0, rids, _ = g
        _, sides, markers = parts[rids[0]]
        centres = []
        for side in (-1, 1):
            mets = [met for marker, sd in sides.items() if sd == side for _, met in markers[marker]]
            mets = [m for m in mets if nodes[m].get("node_is_primary")] or mets
            if mets:
                centres.append((sum(nodes[m]["x"] for m in mets) / len(mets), sum(nodes[m]["y"] for m in mets) / len(mets)))
        if len(centres) < 2:
            continue
        (sx, sy), (px, py) = centres
        a, b, axis = _unit(sx - x0, sy - y0), _unit(px - x0, py - y0), _unit(px - sx, py - sy)
        if not a or not b or not axis or a[0] * b[0] + a[1] * b[1] < math.cos(math.radians(100)):
            continue
        t = min(max(((x0 - sx) * axis[0] + (y0 - sy) * axis[1]) / math.hypot(px - sx, py - sy), 0.0), 1.0)
        fx, fy = sx + t * (px - sx), sy + t * (py - sy)
        dist = math.hypot(x0 - fx, y0 - fy)
        target = offset + (len(rids) - 1) / 2 * lane
        if dist > target + 20:
            candidates.append((dist, g, (fx - x0) / dist * (dist - target), (fy - y0) / dist * (dist - target), axis))
    moved = 0
    for _, g, dx, dy, (ux, uy) in sorted(candidates, key=lambda c: -c[0]):
        x0, y0, rids, _ = g
        group = set(rids)
        mets = {m for rid in rids for m in ends[rid]}
        own = {m for m in mets if users[m] <= group and not nodes[m].get("node_is_primary")}
        skip = group | mets
        nx, ny = _normal(ux, uy)
        away = _unit(-dx, -dy)
        far = 1 if away[0] * nx + away[1] * ny >= 0 else -1
        ids = [reactions[rid]["bigg_id"] for rid in rids]
        # the farthest place with no nodes in the footprint, no more lines in it and of its lines through nodes
        f, base = 0.0, None
        for fraction in (0.0, 1.0, 0.8, 0.6, 0.4, 0.2):
            x, y = x0 + dx * fraction, y0 + dy * fraction
            points = _footprint(x, y, ux, uy, ids, lane, far, reach)
            hard = sum(obstacles.near(px, py, margin, skip) for px, py in points)
            soft = sum(lines.near(px, py, margin, skip) for px, py in points)
            for m in mets - own:
                soft += 3 * sum(obstacles.near(px, py, 20.0, skip)
                                for px, py in _line((x, y), (nodes[m]["x"], nodes[m]["y"]), skip=45.0))
            if base is None:
                base = soft
            elif hard == 0 and soft <= base:
                f = fraction
                break
        if not f:
            continue
        x, y = x0 + dx * f, y0 + dy * f
        g[0], g[1], g[3] = x, y, away
        for rid in rids:
            nodes[parts[rid][0]]["x"], nodes[parts[rid][0]]["y"] = x, y
            obstacles.remove(rid)
            lines.remove(rid)
            add_lines(rid, x, y)
        for m in own:
            obstacles.remove(m)
        for px, py in _footprint(x, y, ux, uy, ids, lane, far, reach):  # later groups keep out of it
            obstacles.add(rids[0], px, py)
        movable.update(own)
        moved += len(rids)
    return moved


def curve_reactions(body: dict, movable=frozenset(), lane: float = 60.0, trunk: float = 0.3,
                    marker_distance: tuple[float, float] = (20.0, 80.0), max_lead: float = 100.0,
                    pull: bool = True) -> dict:
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
    - With `pull`, U-shaped reactions (substrates and products on the same
      side) move towards their metabolites where there is room
      (:func:`_pull_reactions`).
    - The multimarkers are on the axis, `trunk` of the way to the centre
      of their metabolites, within `marker_distance`.
    - The segments between a multimarker and its metabolites are cubic
      Bezier curves (b1, b2): they leave the multimarker along the axis
      and bend towards the metabolite; metabolites on the axis get straight
      segments.
    - Metabolite nodes in `movable` (secondary metabolites that are not
      drawn on the source map) are put next to the multimarker of their
      side, alternately left and right of the axis.

    Other metabolite nodes are not moved. Reactions drawn in other ways are
    left as they are. Returns counts of the curved, the moved (`pull`)
    reactions and the curved segments.
    """
    nodes = body["nodes"]
    parts = {rid: p for rid, r in body["reactions"].items() if (p := _reaction_parts(nodes, r))}

    places = {}
    for rid, (mid, _, _) in parts.items():
        places.setdefault((round(nodes[mid]["x"], 1), round(nodes[mid]["y"], 1)), []).append(rid)
    groups = [[x, y, rids, None] for (x, y), rids in places.items()]
    movable = set(movable)
    moved = _pull_reactions(body, parts, groups, lane, movable) if pull else 0
    segments = 0
    for x0, y0, rids, away in groups:
        axes = {rid: _axis(nodes, *parts[rid]) for rid in rids}
        ux, uy = axes[rids[0]]
        sx = sum(a[0] * math.copysign(1, a[0] * ux + a[1] * uy) for a in axes.values())
        sy = sum(a[1] * math.copysign(1, a[0] * ux + a[1] * uy) for a in axes.values())
        ux, uy = _unit(sx, sy) or (ux, uy)
        nx, ny = _normal(ux, uy)
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
        far = (1 if away[0] * nx + away[1] * ny >= 0 else -1) if away else 1
        half = (k - 1) / 2 * lane
        markers_at = []  # (offset, reaction, multimarker, its direction outwards, its position, movable metabolites)
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
                markers_at.append((offset, rid, marker, (wx, wy), (px, py), free))
        # movable metabolites: alternately left and right of a single reaction; for lanes and moved reactions in a
        # column beside the lanes (the outer lanes' nearest to their multimarkers), away from the line moved towards
        column = {}
        outer = {e[0]: far if away else (1 if e[0] >= 0 else -1) if k > 1 else 0 for e in markers_at}
        for offset, rid, marker, (wx, wy), (px, py), free in sorted(markers_at, key=lambda e: -outer[e[0]] * e[0]):
            for j, m in enumerate(free):
                if outer[offset]:
                    key = (outer[offset], round(wx * ux + wy * uy))
                    q = column.get(key, 0)
                    column[key] = q + 1
                    dist, lat = 45.0 + 32.0 * q, outer[offset] * (half + 55.0) - offset
                else:
                    dist, lat = 45.0 + 30.0 * (j // 2), (55.0 if j % 2 == 0 else -55.0)
                nodes[m]["x"] = px + wx * dist + nx * lat
                nodes[m]["y"] = py + wy * dist + ny * lat
                nodes[m]["label_x"], nodes[m]["label_y"] = nodes[m]["x"] + 8, nodes[m]["y"] - 12
        for offset, rid, marker, (wx, wy), (px, py), free in markers_at:
            for sid, met in parts[rid][2][marker]:
                vx, vy = nodes[met]["x"] - px, nodes[met]["y"] - py
                length = math.hypot(vx, vy) or 1.0
                lead = min(0.4 * length * max((vx * wx + vy * wy) / length, 0.25), max_lead)
                b_marker = dict(x=px + wx * lead, y=py + wy * lead)
                b_met = dict(x=nodes[met]["x"] + (b_marker["x"] - nodes[met]["x"]) * 0.35,
                             y=nodes[met]["y"] + (b_marker["y"] - nodes[met]["y"]) * 0.35)
                s = body["reactions"][rid]["segments"][sid]
                s["b1"], s["b2"] = (b_marker, b_met) if s["from_node_id"] == marker else (b_met, b_marker)
                segments += 1
        labels = _lane_labels(x0, y0, ux, uy, [body["reactions"][rid]["bigg_id"] for rid in rids], lane,
                              far if away else 1)
        for rid, (lx, ly) in zip(rids, labels):
            body["reactions"][rid]["label_x"], body["reactions"][rid]["label_y"] = lx, ly
    return dict(curved_reactions=len(parts), moved_reactions=moved, curved_segments=segments)


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
