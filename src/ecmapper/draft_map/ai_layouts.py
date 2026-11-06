"""AI layouts of whole maps (--split): Escher maps laid out by an AI with a toolkit (the GUI, :mod:`.gui.chat`).

The GUI exports a map as the pipeline draws it (:func:`export_draft`; the
draft objects of the Escher IDs are recorded and keep their IDs in later
exports), an AI lays it out with a toolkit (e.g. escher_handlayout), and the
result is the map's layout (:class:`MapLayout`, a folder next to the
assignment table, sheet "AI layouts"). The pipeline writes its Escher map
and its SVG (exported by Escher in the GUI) as the files of the map; the
positions of the draft nodes in it (:func:`map_positions`) give the node and
box positions of the map's pickle, CNApy maps and reaction positions
(:func:`apply`). A layout holds as long as the map has the same reactions
and subsystems (:func:`check`); the "→ map" labels of metabolites that
continue on other maps follow the assignment (:func:`refresh_exit_labels`).
"""

from __future__ import annotations

import copy
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .core import DraftMap, ModelData
from .escher import place_text_labels, stub_text, to_escher
from .layout import PADDING, TITLE_GAP

LAYOUT, FINAL, SVG = "layout.json", "final.escher.json", "final.svg"
#: the exported map in the work folder of an AI job
DRAFT = "draft.escher.json"


def map_reactions(m: DraftMap) -> dict[str, int]:
    """{reaction ID: how often it is drawn} of a map."""
    return dict(sorted(Counter(r for n in m.nodes.values() if n.kind == "enzyme" for r in n.reactions).items()))


def map_stubs(m: DraftMap) -> dict[str, list[str]]:
    """{metabolite node: the other maps on which it continues}."""
    return {i: list(n.stubs) for i, n in m.nodes.items() if getattr(n, "stubs", None)}


@dataclass(eq=False)  # compared by identity (states of the GUI)
class MapLayout:
    """The AI layout of a map: the Escher map ([head, body]), its SVG export and what it was made from.

    `positions`: {draft node: (x, y)} in map units; `ids`: {Escher ID: draft
    object} (see :func:`.escher.to_escher`; None for layouts made elsewhere
    before they are matched); `stubs`: the "→ map" labels it has.
    """

    map: str
    scale: float
    reactions: dict[str, int]
    subsystems: list[str]
    positions: dict[str, tuple[float, float]]
    escher: list
    svg: str | None = None
    ids: dict[str, list] | None = None
    stubs: dict[str, list[str]] = field(default_factory=dict)
    source: dict = field(default_factory=dict)

    META = ("map", "scale", "reactions", "subsystems", "positions", "ids", "stubs", "source")

    def save(self, folder) -> Path:
        """Write layout.json, final.escher.json and final.svg (if any) to `folder`."""
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        (folder / LAYOUT).write_text(json.dumps({k: getattr(self, k) for k in self.META}, indent=1), encoding="utf-8")
        (folder / FINAL).write_text(json.dumps(self.escher), encoding="utf-8")
        svg = folder / SVG
        if self.svg is not None:
            svg.write_text(self.svg, encoding="utf-8")
        elif svg.exists():
            svg.unlink()
        return folder

    @classmethod
    def load(cls, folder) -> MapLayout:
        folder = Path(folder)
        meta = json.loads((folder / LAYOUT).read_text(encoding="utf-8"))
        meta["positions"] = {i: tuple(p) for i, p in meta["positions"].items()}
        escher = json.loads((folder / FINAL).read_text(encoding="utf-8"))
        svg = (folder / SVG).read_text(encoding="utf-8") if (folder / SVG).is_file() else None
        return cls(escher=escher, svg=svg, **{k: meta[k] for k in cls.META if k in meta})

    def map_json(self, name: str, description: str) -> list:
        """The Escher map with the name `name`."""
        head, body = copy.deepcopy(self.escher)
        head = dict(head, map_name=name, map_id=re.sub(r"\W+", "_", name), map_description=description)
        return [head, body]

    def extent(self) -> tuple[float, float, float, float]:
        """The canvas of the Escher map in map units: x, y, width, height."""
        c = self.escher[1]["canvas"]
        return c["x"] / self.scale, c["y"] / self.scale, c["width"] / self.scale, c["height"] / self.scale


def load_layouts(table, base, log=print) -> dict[str, MapLayout]:
    """The layouts of the sheet "AI layouts" of `table` (folders relative to `base`, the folder of the table)."""
    result = {}
    for name, folder in (table.layouts if table is not None else {}).items():
        path = Path(folder) if Path(folder).is_absolute() else Path(base) / folder
        try:
            result[name] = MapLayout.load(path)
        except (OSError, ValueError, KeyError, TypeError) as e:
            log(f"  warning: the AI layout of map {name!r} could not be read from {path}: {e}")
    return result


def export_draft(m: DraftMap, md: ModelData, name: str, scale: float = 2.0, reuse: dict | None = None,
                 description: str = "Draft map from the BRENDA map, made with ecmapper.draft_map") -> tuple[list, dict]:
    """The Escher map of map `m` (absolute positions) as the pipeline draws it, with the names of the subsystems, and
    the draft objects of its Escher IDs (`reuse`: those of an earlier export, whose objects keep their IDs)."""
    ids = {}
    data = to_escher(m, md, name, scale=scale, subsystem_labels=True, description=description, ids=ids, reuse=reuse)
    return data, ids


def _keepers(run_log: dict | None) -> dict[str, str]:
    """{Escher node removed by a merge of the toolkit: the node it ended up in} (run_log.json of escher_handlayout:
    "merges", [BiGG ID, keeper, removed])."""
    keeper = {str(entry[2]): str(entry[1]) for entry in (run_log or {}).get("merges", [])}

    def last(i):
        seen = set()
        while i in keeper and i not in seen:
            seen.add(i)
            i = keeper[i]
        return i

    return {i: last(i) for i in keeper}


def _mean(points):
    return sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points)


def map_positions(final: list, ids: dict, m: DraftMap, scale: float, run_log: dict | None = None
                  ) -> dict[str, tuple[float, float]]:
    """The positions (map units) of the draft nodes of map `m` in the laid-out Escher map `final`.

    A metabolite node is at (the mean of) its Escher nodes, an enzyme node at
    the mean of the midmarkers of its reactions. Escher nodes merged away by
    the toolkit lead to the node they were merged into (`run_log`); others that
    are gone to the nearest Escher node of the same metabolite (next to the
    node's enzymes).
    """
    nodes = final[1]["nodes"]
    keeper = _keepers(run_log)
    found, gone = defaultdict(list), defaultdict(list)
    for i, key in ids.items():
        if key[0] not in ("node", "mid"):
            continue
        j = i if i in nodes else keeper.get(i)
        if j in nodes:
            found[key[1]].append((nodes[j]["x"] / scale, nodes[j]["y"] / scale))
        elif key[0] == "node":
            gone[key[1]].append(key[2])
    positions = {n: _mean(points) for n, points in found.items()}
    if gone:
        of_metabolite = defaultdict(list)
        for n in nodes.values():
            if n["node_type"] == "metabolite":
                of_metabolite[n["bigg_id"]].append((n["x"] / scale, n["y"] / scale))
        enzymes = defaultdict(list)
        for e in m.edges:
            enzymes[e.met].append(e.enzyme)
        for node, biggs in gone.items():
            near = [positions[e] for e in enzymes.get(node, []) if e in positions]
            candidates = [p for b in biggs for p in of_metabolite.get(b, [])]
            if node not in positions and near and candidates:
                centre = _mean(near)
                positions[node] = min(candidates, key=lambda p: math.dist(p, centre))
    return positions


def _reaction_nodes(body: dict, reaction: dict) -> tuple[str | None, dict[str, list[str]]]:
    """The midmarker and {BiGG ID: metabolite nodes} of an Escher reaction."""
    mid, mets = None, defaultdict(list)
    for s in reaction["segments"].values():
        for i in (s["from_node_id"], s["to_node_id"]):
            n = body["nodes"].get(i)
            if n is None:
                continue
            if n["node_type"] == "midmarker":
                mid = i
            elif n["node_type"] == "metabolite" and i not in mets[n["bigg_id"]]:
                mets[n["bigg_id"]].append(i)
    return mid, mets


def match_ids(draft: list, draft_ids: dict, final: list) -> dict[str, list]:
    """The draft objects of the Escher IDs of `final`, a layout of the map `draft` (exported with `draft_ids`) made
    elsewhere (other Escher IDs): reactions by BiGG ID (several of one in the order of the maps), metabolite nodes by
    the reactions they belong to, "→ map" labels by the nearest of them."""
    d, f = draft[1], final[1]
    by_bigg = defaultdict(lambda: ([], []))
    for i, r in d["reactions"].items():
        by_bigg[r["bigg_id"]][0].append(i)
    for i, r in f["reactions"].items():
        by_bigg[r["bigg_id"]][1].append(i)
    result = {}
    votes = defaultdict(Counter)  # final metabolite node -> draft objects
    for draft_reactions, final_reactions in by_bigg.values():
        for dr, fr in zip(draft_reactions, final_reactions):
            result[fr] = draft_ids[dr]
            d_mid, d_mets = _reaction_nodes(d, d["reactions"][dr])
            f_mid, f_mets = _reaction_nodes(f, f["reactions"][fr])
            if d_mid is not None and f_mid is not None:
                result[f_mid] = draft_ids[d_mid]
            for bigg, f_nodes in f_mets.items():
                for d_node in d_mets.get(bigg, []):
                    key = tuple(draft_ids[d_node])
                    if key[0] == "node":
                        for i in f_nodes:
                            votes[i][key] += 1
    for i, counter in votes.items():
        result[i] = list(counter.most_common(1)[0][0])
    # "→ map" labels: the label of the final map nearest to the metabolite node it belongs to
    at = defaultdict(list)
    for i, key in result.items():
        if key[0] == "node":
            at[key[1]].append((f["nodes"][i]["x"], f["nodes"][i]["y"]))
    arrows = {i: t for i, t in f["text_labels"].items() if t["text"].startswith("→")}
    for i, key in draft_ids.items():
        if key[0] == "stub" and at.get(key[1]) and arrows:
            centre = _mean(at[key[1]])
            j = min(arrows, key=lambda a: math.dist((arrows[a]["x"], arrows[a]["y"]), centre))
            if math.dist((arrows[j]["x"], arrows[j]["y"]), centre) < 800:
                result[j] = list(key)
                del arrows[j]
    return result


def check(layout: MapLayout, m: DraftMap, md: ModelData) -> list[str]:
    """Why the layout no longer fits map `m` (empty: it fits): other reactions or subsystems, other stoichiometry."""

    def listed(items):
        items = sorted(items)
        return ", ".join(items[:10]) + (f" and {len(items) - 10} more" if len(items) > 10 else "")

    problems = []
    now = Counter(map_reactions(m))
    drawn = Counter(r["bigg_id"] for r in layout.escher[1]["reactions"].values())
    if drawn != now:
        missing, extra = now - drawn, drawn - now
        if missing:
            problems.append(f"reactions of the map not in the layout: {listed(missing)}")
        if extra:
            problems.append(f"reactions of the layout not on the map: {listed(extra)}")
    other = []
    for r in layout.escher[1]["reactions"].values():
        rid = r["bigg_id"]
        if rid in md.reactions and rid in now:
            given = {x["bigg_id"]: float(x["coefficient"]) for x in r["metabolites"]}
            model = {k: float(v) for k, v in md.reactions[rid]["metabolites"].items()}
            if given != model:
                other.append(rid)
    if other:
        problems.append(f"reactions with another stoichiometry than in the model: {listed(set(other))}")
    if set(layout.subsystems) != set(m.subsystems):
        added, removed = set(m.subsystems) - set(layout.subsystems), set(layout.subsystems) - set(m.subsystems)
        if added:
            problems.append(f"subsystems added to the map: {listed(added)}")
        if removed:
            problems.append(f"subsystems taken off the map: {listed(removed)}")
    return problems


def make_layout(m: DraftMap, md: ModelData, name: str, final: list, ids: dict, scale: float, svg: str | None = None,
                run_log: dict | None = None, source: dict | None = None) -> MapLayout:
    """The layout of map `m` (exported with `ids`) laid out as `final` (a copy, with the "→ map" labels of `m`)."""
    labelled = {key[1]: [] for i, key in ids.items() if key[0] == "stub" and i in final[1]["text_labels"]}
    layout = MapLayout(map=name, scale=scale, reactions=map_reactions(m), subsystems=sorted(m.subsystems),
                       positions=map_positions(final, ids, m, scale, run_log), escher=copy.deepcopy(final), svg=svg,
                       ids=dict(ids), stubs=labelled, source=dict(source or {}))
    refresh_exit_labels(layout, m)
    return layout


def refresh_exit_labels(layout: MapLayout, m: DraftMap) -> bool:
    """Make the "→ map" labels of the layout those of map `m` (the assignment of other maps may have changed).

    Labels of metabolites that no longer continue on other maps are removed,
    those of metabolites that do now (`stubs` of the layout: those that had a
    label) are added next to their nodes. Returns whether the labels changed
    (the SVG export is then outdated).
    """
    now = map_stubs(m)
    if layout.ids is None:
        return False
    body = layout.escher[1]
    labels = body["text_labels"]
    labelled, changed = set(), False
    for i, key in list(layout.ids.items()):
        if key[0] != "stub" or i not in labels:
            continue
        if key[1] in now:
            text = stub_text(now[key[1]])
            changed |= labels[i]["text"] != text
            labels[i]["text"] = text
            labelled.add(key[1])
        else:
            del labels[i]
            changed = True
    at = defaultdict(list)
    for i, key in layout.ids.items():
        if key[0] == "node" and i in body["nodes"]:
            at[key[1]].append(body["nodes"][i])
    used = list(body["nodes"]) + list(body["reactions"]) + list(labels)
    used += [s for r in body["reactions"].values() for s in r["segments"]]
    n = max((int(i) for i in used if str(i).isdigit()), default=0)
    anchors = {}
    for node, maps in now.items():
        if node in labelled or node in layout.stubs or not at.get(node):
            continue
        n += 1
        label = str(n)
        xs, ys = [p["x"] for p in at[node]], [p["y"] for p in at[node]]
        labels[label] = dict(text=stub_text(maps), x=xs[0], y=max(ys) + 66)
        layout.ids[label] = ["stub", node]
        anchors[label] = (sum(xs) / len(xs), min(ys), max(ys))
    if anchors:
        place_text_labels(body, anchors)
    layout.stubs = now
    return changed or bool(anchors)


def apply(m: DraftMap, layout: MapLayout) -> int:
    """Put the nodes of map `m` where the layout has them and the subsystem boxes around their nodes (absolute).

    Nodes the layout does not place keep their offset from a placed
    neighbour. Boxes may overlap (subsystems of a layout may be
    interleaved). Returns the number of nodes not placed by the layout.
    """
    positions = layout.positions
    before = {i: (n.x, n.y) for i, n in m.nodes.items()}
    neighbours = defaultdict(list)
    for e in m.edges:
        neighbours[e.met].append(e.enzyme)
        neighbours[e.enzyme].append(e.met)
    for i, n in m.nodes.items():
        if i in positions:
            n.x, n.y = positions[i]
    unplaced = 0
    for i, n in m.nodes.items():
        if i in positions:
            continue
        unplaced += 1
        ref = next((j for j in neighbours.get(i, []) if j in positions), None)
        if ref is not None:
            n.x = positions[ref][0] + before[i][0] - before[ref][0]
            n.y = positions[ref][1] + before[i][1] - before[ref][1]
    for s in m.subsystems.values():
        members = m.members(s.id)
        if not members:
            continue
        x0, x1 = min(n.x - n.w / 2 for n in members), max(n.x + n.w / 2 for n in members)
        y0, y1 = min(n.y - n.h / 2 for n in members), max(n.y + n.h / 2 for n in members)
        top = s.title_h / 2 + TITLE_GAP
        s.w = max(x1 - x0 + 2 * PADDING, s.title_w + 2 * PADDING)
        s.h = y1 - y0 + top + PADDING
        s.x, s.y = (x0 + x1) / 2 - s.w / 2, y0 - top
        s.flip_x = s.flip_y = False
    return unplaced
