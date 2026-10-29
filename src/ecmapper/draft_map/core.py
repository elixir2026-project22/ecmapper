"""Draft maps from the BRENDA map: data model, pruning and completion with model reactions."""

from __future__ import annotations

import copy
import json
import math
import pickle
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

#: short labels of currency metabolites for added cofactor nodes (BiGG IDs without compartment)
SHORT_NAMES = {
    "atp": "ATP", "adp": "ADP", "amp": "AMP", "gtp": "GTP", "gdp": "GDP", "gmp": "GMP", "ctp": "CTP",
    "cdp": "CDP", "cmp": "CMP", "utp": "UTP", "udp": "UDP", "ump": "UMP", "itp": "ITP", "idp": "IDP",
    "nad": "NAD+", "nadh": "NADH", "nadp": "NADP+", "nadph": "NADPH", "fad": "FAD", "fadh2": "FADH2",
    "h2o": "H2O", "co2": "CO2", "pi": "Pi", "ppi": "PPi", "coa": "CoA", "nh4": "NH4+", "o2": "O2",
    "h2o2": "H2O2", "q8": "ubiquinone", "q8h2": "ubiquinol", "mqn8": "menaquinone", "mql8": "menaquinol",
    "2dmmq8": "DMK", "2dmmql8": "DMKH2", "hco3": "HCO3-", "glu__L": "Glu", "akg": "2-oxoglutarate",
    "gln__L": "Gln", "accoa": "acetyl-CoA", "amet": "SAM", "ahcys": "SAH", "fdxrd": "Fd_red",
    "fdxox": "Fd_ox", "flxr": "Fld_red", "flxso": "Fld_ox", "trdrd": "Trx_red", "trdox": "Trx_ox",
    "grxrd": "Grx_red", "grxox": "Grx_ox", "gthrd": "GSH", "gthox": "GSSG", "h2s": "H2S", "so3": "sulfite",
    "acp": "ACP", "thf": "THF", "mlthf": "5,10-methylene-THF", "10fthf": "10-formyl-THF", "5mthf": "5-methyl-THF",
}
PROTONS = {"h"}


def text_width(text: str, font_size: float = 12.0) -> float:
    """Approximate width of a label in map units (Arial), as the BRENDA map sizes its nodes."""
    return 0.54 * font_size * len(text) + 10.0


@dataclass
class Node:
    id: str
    kind: str  # enzyme, main, cofactor
    subsystem: str
    x: float
    y: float
    w: float
    h: float
    label: str
    title: str = ""
    svg_class: str = ""  # nodeEN, nodeHM, nodeSM, nodeCM
    ligand_id: int | None = None
    ec_label: str = ""
    reactions: list[str] = field(default_factory=list)  # model reactions of enzyme nodes
    status: str = "map"  # enzyme nodes: on-map, vicinity, added; metabolite nodes: map, added
    model_mets: set[str] = field(default_factory=set)  # metabolite nodes: model metabolites (without compartment)

    @property
    def box(self) -> tuple[float, float, float, float]:
        return self.x - self.w / 2, self.y - self.h / 2, self.x + self.w / 2, self.y + self.h / 2


@dataclass
class Edge:
    """A link between an enzyme node and a metabolite node."""

    enzyme: str
    met: str
    kind: str  # main or cofactor
    arrow_enzyme: bool
    arrow_met: bool
    model: bool = False  # direction taken from the model


@dataclass
class MetLink:
    """A link between two drawn nodes of the same metabolite (in different subsystems)."""

    a: str
    b: str
    origin: str = "map"  # map, bridged (through removed nodes) or added


@dataclass
class Subsystem:
    id: str
    label: str
    color_group: str
    fill: str
    x: float  # box, top left
    y: float
    w: float
    h: float
    title_w: float = 0.0
    title_h: float = 31.0
    flip_x: bool = False
    flip_y: bool = False
    aside: str = ""  # set aside by drop_redundant_subsystems: the rule


class DraftMap:
    """A (part of the) BRENDA map with the model reactions of its enzyme nodes."""

    def __init__(self):
        self.subsystems: dict[str, Subsystem] = {}
        self.nodes: dict[str, Node] = {}
        self.edges: list[Edge] = []
        self.links: list[MetLink] = []
        self.log: list[str] = []
        self.stats: dict = {}

    # --- loading ------------------------------------------------------------------------------

    @classmethod
    def from_brenda(cls, json_path: str, svg_path: str | None = None) -> DraftMap:
        with open(json_path, encoding="utf-8") as f:
            data = json.load(f)
        fills = {}
        if svg_path:
            text = open(svg_path, encoding="utf-8").read()
            for m in re.finditer(r'<polygon id="pw_([^"]+)"[^>]*style="([^"]*)"', text):
                fill = re.search(r"fill:\s*([^;]+)", m.group(2))
                fills[m.group(1)] = fill.group(1).strip() if fill else "#dedede"
        dm = cls()
        for s in data["subsystems"]:
            b, t = s["bbox"], s.get("title_node") or {}
            dm.subsystems[s["id"]] = Subsystem(
                s["id"], s.get("label_text") or (t.get("title") or s["id"]), s.get("color_group", ""),
                fills.get(s["id"], "#dedede"), b["x"], b["y"], b["width"], b["height"],
                t.get("width", text_width(s.get("label_text", ""), 18)), t.get("height", 31.0),
            )
        for m in data["metabolite_nodes"]:
            kind = "cofactor" if m["node_class"] == "nodeCM" else "main"
            dm.nodes[m["node_id"]] = Node(
                m["node_id"], kind, m["subsystem"], m["x"], m["y"], m["width"], m["height"],
                m.get("display_name") or m.get("brenda_name") or "", (m.get("brenda_name") or "").strip(),
                m["node_class"], m.get("brenda_ligand_id"),
            )
        for r in data["reactions"]:
            lay = r["layout"]
            dm.nodes[r["node_id"]] = Node(
                r["node_id"], "enzyme", r["subsystem"], lay["x"], lay["y"], lay["width"], lay["height"],
                r.get("ec_label") or "", r.get("enzyme_name") or "", "nodeEN", ec_label=r.get("ec_label") or "",
            )
            for kind, items in (
                ("main", r["main_metabolites"]),
                ("consumed", r["cofactors"]["consumed"]),
                ("produced", r["cofactors"]["produced"]),
                ("undirected", r["cofactors"]["undirected"]),
            ):
                for m in items:
                    if m["node_id"] not in dm.nodes:
                        continue
                    if kind == "main":
                        link = m.get("link") or {}
                        start, end = link.get("start"), link.get("end")
                        if start and end:
                            # which end of the line touches the enzyme
                            ds = math.dist(start, (lay["x"], lay["y"]))
                            de = math.dist(end, (lay["x"], lay["y"]))
                            at_enzyme, at_met = (
                                (link.get("arrow_at_start"), link.get("arrow_at_end"))
                                if ds < de else (link.get("arrow_at_end"), link.get("arrow_at_start"))
                            )
                        else:
                            at_enzyme = at_met = True
                        dm.edges.append(Edge(r["node_id"], m["node_id"], "main", bool(at_enzyme), bool(at_met)))
                    else:
                        dm.edges.append(
                            Edge(r["node_id"], m["node_id"], "cofactor", kind == "consumed", kind == "produced")
                        )
        for link in data["metabolite_links"]:
            a, b = link["start_node_id"], link["end_node_id"]
            if a != b and a in dm.nodes and b in dm.nodes:
                dm.links.append(MetLink(a, b))
        # metabolite nodes drawn in another subsystem than all their enzymes belong to those
        enzymes_of = defaultdict(set)
        for e in dm.edges:
            enzymes_of[e.met].add(dm.nodes[e.enzyme].subsystem)
        for met, subs in enzymes_of.items():
            if len(subs) == 1 and dm.nodes[met].subsystem not in subs:
                dm.nodes[met].subsystem = next(iter(subs))
        return dm

    # --- helpers ------------------------------------------------------------------------------

    def members(self, subsystem: str) -> list[Node]:
        return [n for n in self.nodes.values() if n.subsystem == subsystem]

    def edges_of(self) -> dict[str, list[Edge]]:
        result = defaultdict(list)
        for e in self.edges:
            result[e.enzyme].append(e)
            result[e.met].append(e)
        return result

    def remove_nodes(self, ids: set[str]) -> None:
        for i in ids:
            self.nodes.pop(i, None)
        self.edges = [e for e in self.edges if e.enzyme in self.nodes and e.met in self.nodes]

    def summary(self) -> dict:
        kinds = Counter(n.kind for n in self.nodes.values())
        status = Counter(n.status for n in self.nodes.values() if n.kind == "enzyme")
        reactions = {r for n in self.nodes.values() for r in n.reactions}
        inter = sum(1 for l in self.links if self.nodes[l.a].subsystem != self.nodes[l.b].subsystem)
        return dict(
            subsystems=len(self.subsystems), enzyme_nodes=kinds["enzyme"], main_metabolite_nodes=kinds["main"],
            cofactor_nodes=kinds["cofactor"], edges=len(self.edges), inter_pathway_links=inter,
            model_reactions=len(reactions), **{f"enzyme_nodes_{k}": v for k, v in status.items()},
        )


def base_id(met_id: str, compartment: str | None) -> str:
    """Metabolite ID without compartment suffix (glc__D_c -> glc__D)."""
    if compartment and met_id.endswith("_" + compartment):
        return met_id[: -len(compartment) - 1]
    return re.sub(r"_[a-z][a-z0-9]?$", "", met_id)


class ModelData:
    """What :func:`.prepare.prepare` found about a model and its matches to the map (the dict or a pickle of it)."""

    def __init__(self, data: dict | str):
        if not isinstance(data, dict):
            with open(data, "rb") as f:
                data = pickle.load(f)
        self.__dict__.update(data)
        self.base = {m: base_id(m, v["compartment"]) for m, v in self.metabolites.items()}

    def is_main(self, met_id: str) -> bool:
        m = self.metabolites[met_id]
        return bool(m["keys"]) and not m["currency"] and self.base[met_id] not in PROTONS

    def accepted_nodes(self) -> dict[str, list[tuple[str, bool]]]:
        """{node ID: [(reaction ID, found only in the vicinity)]} of the accepted matches."""
        result = defaultdict(list)
        for rid, m in self.matches.items():
            for c in m["candidates"]:
                if c["accepted"]:
                    result[c["node_id"]].append((rid, "vicinity" in c["via"]))
        return result


# --- steps ----------------------------------------------------------------------------------


def drop_empty_subsystems(dm: DraftMap, md: ModelData) -> None:
    """Step 1: keep the subsystems with at least one enzyme node with model reactions."""
    accepted = md.accepted_nodes()
    keep = {dm.nodes[n].subsystem for n in accepted if n in dm.nodes}
    before = len(dm.subsystems)
    dm.subsystems = {k: v for k, v in dm.subsystems.items() if k in keep}
    dm.remove_nodes({n.id for n in dm.nodes.values() if n.subsystem not in keep})
    dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]
    for node_id, items in accepted.items():
        node = dm.nodes[node_id]
        node.reactions = sorted({r for r, _ in items})
        node.status = "vicinity" if all(v for _, v in items) else "on-map"
    dm.log.append(f"step 1: {len(keep)} of {before} subsystems have model reactions")


def drop_redundant_subsystems(dm: DraftMap, drop_subsets: bool = True, coverage: float = 0.9,
                              min_matched: float = 0.0) -> dict:
    """Step 1, continued: drop or set aside subsystems whose model reactions are drawn elsewhere.

    Model reactions of a subsystem: those of its enzyme nodes (after step 1).

    - Rule A (`drop_subsets`): a subsystem whose model reactions are a proper
      subset of those of another subsystem is removed.
    - Rule B/C (`coverage`): going through the remaining subsystems from the
      smallest (fewest model reactions), a subsystem of which at least the
      fraction `coverage` of the model reactions are in the other remaining
      subsystems is set aside (B: all of them, C: the fraction). ``1.0``:
      only B; above 1: off.
    - Rule D (`min_matched`): a subsystem in which less than this fraction of
      the enzyme nodes have model reactions is set aside.

    Set-aside subsystems stay in the map with ``aside`` set to the rule; they
    are drawn on a separate map and are not used by the completion (step 3).
    Returns the dropped and set-aside subsystems.
    """
    reactions = {k: frozenset(r for n in dm.members(k) for r in n.reactions) for k in dm.subsystems}
    dropped = []
    if drop_subsets:
        for k, rs in sorted(reactions.items()):
            for other, ors in sorted(reactions.items()):
                if other != k and rs < ors:
                    dropped.append(dict(subsystem=k, reactions=sorted(rs), subset_of=other))
                    break
        gone = {d["subsystem"] for d in dropped}
        dm.subsystems = {k: v for k, v in dm.subsystems.items() if k not in gone}
        dm.remove_nodes({n.id for n in dm.nodes.values() if n.subsystem in gone})
        dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]
    aside = []
    remaining = {k: reactions[k] for k in dm.subsystems}
    if coverage <= 1.0:
        for k in sorted(remaining, key=lambda k: (len(reactions[k]), k)):
            others = set().union(*(rs for o, rs in remaining.items() if o != k))
            covered = len(reactions[k] & others) / len(reactions[k])
            if covered >= coverage:
                del remaining[k]
                rule = "B" if covered == 1.0 else "C"
                dm.subsystems[k].aside = rule
                aside.append(dict(subsystem=k, rule=rule, reactions=len(reactions[k]), covered=covered))
    if min_matched > 0:
        for k in sorted(remaining):
            enzymes = [n for n in dm.members(k) if n.kind == "enzyme"]
            matched = sum(bool(n.reactions) for n in enzymes) / len(enzymes)
            if matched < min_matched:
                dm.subsystems[k].aside = "D"
                aside.append(dict(subsystem=k, rule="D", reactions=len(reactions[k]), matched=matched))
    rules = Counter(a["rule"] for a in aside)
    dm.log.append(
        f"step 1: removed {len(dropped)} subsystems whose model reactions are a proper subset of another's; set "
        f"aside {len(aside)} (rule B: {rules['B']}, C: {rules['C']}, D: {rules['D']})"
    )
    return dict(dropped=dropped, set_aside=aside)


def anchor_metabolites(dm: DraftMap, md: ModelData) -> None:
    """Identify the model metabolites of the drawn metabolite nodes of enzyme nodes with model reactions.

    A drawn metabolite of an enzyme node is a model metabolite of one of its
    model reactions if they share a key (KEGG compound or normalised name), as
    in the matching. Also sets the arrows of the main metabolites from the
    directions of the model reactions.
    """
    for e in dm.edges:
        enzyme, met = dm.nodes[e.enzyme], dm.nodes[e.met]
        if not enzyme.reactions:
            continue
        keys = md.node_met_keys.get((enzyme.id, met.id), frozenset())
        if not keys:
            continue
        sides = set()
        for rid in enzyme.reactions:
            r = md.reactions[rid]
            for mid, coefficient in r["metabolites"].items():
                mk = md.metabolites[mid]["keys"]
                if mk & keys:
                    if e.kind == "main" and md.is_main(mid):
                        met.model_mets.add(md.base[mid])
                    sides.add("both" if r["reversible"] else ("product" if coefficient > 0 else "substrate"))
        if e.kind == "main" and sides:
            if len(sides) == 1 and "both" not in sides:
                e.arrow_met, e.arrow_enzyme = sides == {"product"}, sides == {"substrate"}
            else:
                e.arrow_met = e.arrow_enzyme = True
            e.model = True


def prune_unassociated(dm: DraftMap) -> None:
    """Step 2: remove enzyme nodes without model reactions and the metabolites only they have.

    Links between drawn metabolites whose nodes are removed are bridged: if
    the remaining nodes of a metabolite were connected through removed nodes,
    they are connected directly.
    """
    enzymes = {n.id for n in dm.nodes.values() if n.kind == "enzyme" and not n.reactions}
    before = sum(1 for n in dm.nodes.values() if n.kind == "enzyme")
    dm.remove_nodes(enzymes)
    attached = {e.met for e in dm.edges}
    linked_ok = attached
    mets = {n.id for n in dm.nodes.values() if n.kind != "enzyme" and n.id not in attached}
    # bridge links through removed metabolite nodes
    adjacency = defaultdict(set)
    for l in dm.links:
        adjacency[l.a].add(l.b)
        adjacency[l.b].add(l.a)
    new_links, seen = [], set()
    for l in dm.links:
        if l.a in linked_ok and l.b in linked_ok:
            key = frozenset((l.a, l.b))
            if key not in seen:
                seen.add(key)
                new_links.append(l)
    for start in sorted(linked_ok & set(adjacency)):
        # search through removed nodes only
        stack, visited = [start], {start}
        while stack:
            current = stack.pop()
            for nxt in adjacency[current]:
                if nxt in visited:
                    continue
                visited.add(nxt)
                if nxt in linked_ok:
                    if current != start:
                        key = frozenset((start, nxt))
                        if key not in seen:
                            seen.add(key)
                            new_links.append(MetLink(start, nxt, "bridged"))
                else:
                    stack.append(nxt)
    dm.remove_nodes(mets)
    dm.links = [l for l in new_links if l.a in dm.nodes and l.b in dm.nodes]
    after = sum(1 for n in dm.nodes.values() if n.kind == "enzyme")
    dm.log.append(
        f"step 2: removed {before - after} of {before} enzyme nodes without model reactions and "
        f"{len(mets)} metabolite nodes; {sum(l.origin == 'bridged' for l in dm.links)} links bridged"
    )


# --- completion with model reactions -----------------------------------------------------------


def _overlaps(box, boxes, margin):
    x0, y0, x1, y1 = box
    return any(x0 < b[2] + margin and b[0] - margin < x1 and y0 < b[3] + margin and b[1] - margin < y1 for b in boxes)


def place_free(x, y, w, h, boxes, margin=12.0, step=45.0, prefer=None):
    """A position near (x, y) where a w x h box does not overlap `boxes` (spiral search)."""
    candidates = [(x, y)]
    directions = 16
    for ring in range(1, 40):
        r = ring * step
        ring_candidates = []
        for k in range(directions):
            a = 2 * math.pi * k / directions
            ring_candidates.append((x + r * math.cos(a) * 1.6, y + r * math.sin(a)))
        if prefer is not None:
            ring_candidates.sort(key=lambda p: math.dist(p, prefer))
        candidates.extend(ring_candidates)
    for cx, cy in candidates:
        box = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        if not _overlaps(box, boxes, margin):
            return cx, cy
    return x, y


#: carriers and redox partners that are not main metabolites of added reactions (BiGG IDs, no compartment)
CARRIERS = {
    "ACP", "apoACP", "trdrd", "trdox", "grxrd", "grxox", "fdxrd", "fdxox", "flxr", "flxso", "gthrd", "gthox",
    "q8", "q8h2", "mqn8", "mql8", "2dmmq8", "2dmmql8", "h2o", "co2", "o2", "nh4", "pi", "ppi", "hco3", "h2o2",
}

DRAWN = ("map", "linked")


def _stem(identifier: str) -> str:
    """ID without digits: variants of reactions and metabolites with other chain lengths get the same stem."""
    return re.sub(r"\d+", "", identifier)  # metabolite nodes of metabolites drawn on the map (not new ones)


class Completion:
    """Step 3: add model reactions that are not on the map but connect metabolites drawn on it.

    Candidates are the model reactions that are not on the map, no boundary
    reactions or transporters, and whose main metabolites (no currency
    metabolites, carriers such as ACP, protons) are in one compartment. A
    candidate is added if at least `min_anchored` of its main metabolites are
    already drawn, at least one of them on the map itself (not only on a node
    added for an earlier reaction): the reaction closes a gap of the map or
    extends it by one step.

    It is drawn in the subsystem in which most of these metabolites are drawn
    (ties: the subsystem of the best related node, the subsystem in which the
    model reactions of the same model subsystem are drawn most often, the
    subsystem whose original position is nearest). Main metabolites drawn only
    in other subsystems get a new node there, linked to their nearest drawn
    node; main metabolites not drawn at all get a new node. Reactions that
    connect the same drawn nodes in the same way and are variants of each
    other (IDs of the reactions or of their new metabolites equal apart from
    digits, e.g. other chain lengths) share one enzyme node. Repeated until no reaction is added.
    """

    def __init__(self, dm: DraftMap, md: ModelData, min_anchored: int = 2, cofactors: bool = True,
                 origin_positions: dict[str, tuple[float, float, str]] | None = None):
        self.dm, self.md = dm, md
        self.min_anchored = min_anchored
        self.cofactors = cofactors
        self.counter = 0
        self.origin_positions = origin_positions or {}
        self.added: list[dict] = []
        self.groups: dict[tuple, list[Node]] = {}
        self.centres = {}
        for s in dm.subsystems.values():
            members = dm.members(s.id)
            self.centres[s.id] = (sum(n.x for n in members) / len(members), sum(n.y for n in members) / len(members))

    def is_main(self, met_id: str) -> bool:
        return self.md.is_main(met_id) and self.md.base[met_id] not in CARRIERS

    def anchors(self) -> dict[str, list[str]]:
        result = defaultdict(list)
        for n in self.dm.nodes.values():
            if n.kind == "main" and not self.dm.subsystems[n.subsystem].aside:
                for m in n.model_mets:
                    result[m].append(n.id)
        return result

    def candidates(self) -> list[str]:
        # reactions drawn only in set-aside subsystems are candidates, too
        on_map = {r for n in self.dm.nodes.values() if not self.dm.subsystems[n.subsystem].aside for r in n.reactions}
        result = []
        for rid in self.md.matches:
            r = self.md.reactions[rid]
            if rid in on_map or r["boundary"] or r["transport"]:
                continue
            mains = [m for m in r["metabolites"] if self.is_main(m)]
            if len(mains) < 2 or len({self.md.metabolites[m]["compartment"] for m in mains}) > 1:
                continue
            result.append(rid)
        return sorted(result)

    def main_mets(self, rid: str) -> dict[str, float]:
        return {self.md.base[m]: c for m, c in self.md.reactions[rid]["metabolites"].items() if self.is_main(m)}

    def model_subsystem_votes(self) -> dict[str, Counter]:
        votes = defaultdict(Counter)
        for n in self.dm.nodes.values():
            for rid in n.reactions:
                votes[self.md.reactions[rid]["subsystem"]][n.subsystem] += 1
        return votes

    def run(self) -> list[dict]:
        candidates = self.candidates()
        votes = self.model_subsystem_votes()
        done = set()
        while True:
            added_now = 0
            for rid in candidates:
                if rid in done:
                    continue
                anchors = self.anchors()
                mains = self.main_mets(rid)
                anchored = {m: anchors[m] for m in mains if anchors.get(m)}
                drawn = [m for m, nodes in anchored.items() if any(self.dm.nodes[n].status == "map" for n in nodes)]
                if len(anchored) < self.min_anchored or not drawn:
                    continue
                self.add(rid, mains, anchored, votes)
                done.add(rid)
                added_now += 1
            if not added_now:
                break
        nodes = {a["node"] for a in self.added}
        self.dm.log.append(
            f"step 3: added {len(self.added)} of {len(candidates)} candidate model reactions (not on the map, no "
            f"boundary reactions, transporters) on {len(nodes)} new enzyme nodes, "
            f"{sum(a['new_metabolites'] for a in self.added)} nodes of metabolites not on the map, "
            f"{sum(a['linked_metabolites'] for a in self.added)} linked nodes of metabolites drawn in other subsystems"
        )
        return self.added

    def new_id(self, prefix: str) -> str:
        self.counter += 1
        return f"{prefix}{self.counter}"

    def target_subsystem(self, rid, anchored, votes) -> str:
        dm, md = self.dm, self.md
        per_sub = defaultdict(lambda: [0, 0])
        for m, nodes in anchored.items():
            subs = {dm.nodes[n].subsystem for n in nodes}
            drawn = {dm.nodes[n].subsystem for n in nodes if dm.nodes[n].status in DRAWN}
            for s in subs:
                per_sub[s][0] += 1
                per_sub[s][1] += s in drawn
        related = {
            self.origin_positions[c["node_id"]][2] for c in md.matches[rid]["candidates"]
            if set(c["via"]) & {"reaction-id", "ec"} and c["node_id"] in self.origin_positions
        }
        vote = votes.get(md.reactions[rid]["subsystem"], Counter()) if md.reactions[rid]["subsystem"] else Counter()
        pts = [(dm.nodes[n].x, dm.nodes[n].y) for nodes in anchored.values() for n in nodes]
        cx, cy = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
        return min(
            per_sub,
            key=lambda s: (-per_sub[s][0], -per_sub[s][1], s not in related, -vote.get(s, 0),
                           math.dist(self.centres[s], (cx, cy)), s),
        )

    def add(self, rid, mains, anchored, votes) -> None:
        dm, md = self.dm, self.md
        reaction = md.reactions[rid]
        target = self.target_subsystem(rid, anchored, votes)
        # the node of each metabolite drawn in the target subsystem, the one closest to the others
        local = {m: [n for n in nodes if dm.nodes[n].subsystem == target] for m, nodes in anchored.items()}
        local = {m: nodes for m, nodes in local.items() if nodes}
        chosen = {m: nodes[0] for m, nodes in local.items() if len(nodes) == 1}
        for m, nodes in local.items():
            if len(nodes) > 1:
                ref = [dm.nodes[n] for n in chosen.values()] or [dm.nodes[n] for n in nodes]
                cx = sum(n.x for n in ref) / len(ref)
                cy = sum(n.y for n in ref) / len(ref)
                chosen[m] = min(nodes, key=lambda n: math.dist((dm.nodes[n].x, dm.nodes[n].y), (cx, cy)))
        # reactions connecting the same nodes in the same way share an enzyme node
        sub_ids = frozenset(chosen[m] for m, c in mains.items() if m in chosen and c < 0)
        prod_ids = frozenset(chosen[m] for m, c in mains.items() if m in chosen and c > 0)
        n_sub = sum(1 for m, c in mains.items() if m not in chosen and c < 0)
        n_prod = sum(1 for m, c in mains.items() if m not in chosen and c > 0)
        key = (target, sub_ids, prod_ids, n_sub, n_prod)
        reverse = (target, prod_ids, sub_ids, n_prod, n_sub)
        new_stems = sorted(_stem(m) for m in mains if m not in chosen)
        for k, flipped in ((key, False), (reverse, True)):
            if flipped and not reaction["reversible"]:
                continue
            for group in self.groups.get(k, []):
                # variants (e.g. other chain lengths): same IDs apart from digits
                if _stem(group.reactions[0]) == _stem(rid) or (new_stems and group.new_stems == new_stems):
                    self.merge(group, rid, mains, chosen, flipped=flipped)
                    return
        members = dm.members(target)
        boxes = [n.box for n in members]
        position = None
        for c in md.matches[rid]["candidates"]:
            o = self.origin_positions.get(c["node_id"])
            if set(c["via"]) & {"reaction-id", "ec"} and o and o[2] == target:
                position = o[:2]  # the related node of the map, removed in step 2
                break
        if position is None:
            pts = [(dm.nodes[n].x, dm.nodes[n].y) for n in chosen.values()]
            if len(pts) >= 2:
                position = (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))
            else:
                cx, cy = self.centres[target]
                px, py = pts[0]
                d = math.dist((px, py), (cx, cy)) or 1.0
                position = (px + 120 * (px - cx) / d, py + 90 * (py - cy) / d)
        enzyme = Node(
            self.new_id("added_r"), "enzyme", target, 0, 0, text_width(rid) + 6, 30.0, rid,
            reaction["name"], "nodeEN", ec_label=";".join(reaction["ecs"][:3]), reactions=[rid], status="added",
        )
        enzyme.x, enzyme.y = place_free(position[0], position[1], enzyme.w, enzyme.h, boxes)
        dm.nodes[enzyme.id] = enzyme
        boxes.append(enzyme.box)
        enzyme.new_nodes = {"sub": [], "prod": []}
        reversible = reaction["reversible"]
        new_mets = linked = 0
        anchors = self.anchors()
        for m, coefficient in mains.items():
            if m in chosen:
                met_id = chosen[m]
            else:
                elsewhere = anchors.get(m, [])
                mid = next(k for k in reaction["metabolites"] if md.base[k] == m)
                source = None
                if elsewhere:
                    # the drawn node nearest to the target subsystem
                    tc = self.centres[target]
                    source = min((dm.nodes[n] for n in elsewhere),
                                 key=lambda n: (n.status not in DRAWN, math.dist((n.x, n.y), tc)))
                    label = getattr(source, "base_label", source.label)
                    title, ligand, css = source.title, source.ligand_id, source.svg_class
                else:
                    label = title = md.metabolites[mid]["name"] or m
                    ligand, css = None, "nodeHM"
                met = Node(self.new_id("added_m"), "main", target, 0, 0, text_width(label), 30.0, label, title, css,
                           ligand, status="linked" if source else "new", model_mets={m})
                met.base_label = label
                met.model_met_ids = {mid}
                prefer = (source.x, source.y) if source else None
                met.x, met.y = place_free(
                    enzyme.x + (140 if coefficient > 0 else -140), enzyme.y + (60 if coefficient > 0 else -60),
                    met.w, met.h, boxes, step=40, prefer=prefer,
                )
                dm.nodes[met.id] = met
                boxes.append(met.box)
                met_id = met.id
                enzyme.new_nodes["prod" if coefficient > 0 else "sub"].append(met.id)
                if source:
                    dm.links.append(MetLink(met.id, source.id, "added"))
                    linked += 1
                else:
                    new_mets += 1
            dm.edges.append(Edge(enzyme.id, met_id, "main", reversible or coefficient < 0,
                                 reversible or coefficient > 0, model=True))
        if self.cofactors:
            for mid, coefficient in reaction["metabolites"].items():
                b = md.base[mid]
                if self.is_main(mid) or b in PROTONS:
                    continue
                label = SHORT_NAMES.get(b, md.metabolites[mid]["name"] or b)
                cof = Node(self.new_id("added_c"), "cofactor", target, 0, 0, text_width(label), 30.0, label,
                           md.metabolites[mid]["name"], "nodeCM", status="added")
                cof.model_met_ids = {mid}
                cof.x, cof.y = place_free(
                    enzyme.x + (70 if coefficient > 0 else -70), enzyme.y + 45, cof.w, cof.h, boxes, margin=6, step=30
                )
                dm.nodes[cof.id] = cof
                boxes.append(cof.box)
                dm.edges.append(Edge(enzyme.id, cof.id, "cofactor", coefficient < 0 or reversible,
                                     coefficient > 0 or reversible, model=True))
        enzyme.new_stems = new_stems
        self.groups.setdefault(key, []).append(enzyme)
        self.added.append(dict(
            reaction=rid, name=reaction["name"], status=md.matches[rid]["status"], subsystem=target,
            model_subsystem=reaction["subsystem"], anchored=sorted(anchored), main_metabolites=sorted(mains),
            new_metabolites=new_mets, linked_metabolites=linked, node=enzyme.id,
        ))

    def merge(self, enzyme: Node, rid, mains, chosen, flipped: bool) -> None:
        """Add a reaction to the enzyme node of a reaction that connects the same nodes."""
        md = self.md
        enzyme.reactions.append(rid)
        enzyme.label = f"{enzyme.reactions[0]} (+{len(enzyme.reactions) - 1})"
        enzyme.w = text_width(enzyme.label) + 6
        sides = {"sub": [], "prod": []}
        for m, c in mains.items():
            if m not in chosen:
                side = "prod" if (c > 0) != flipped else "sub"
                sides[side].append(m)
        for side, mets in sides.items():
            for m, node_id in zip(mets, enzyme.new_nodes[side]):
                node = self.dm.nodes[node_id]
                mid = next(x for x in md.reactions[rid]["metabolites"] if md.base[x] == m)
                node.model_met_ids = getattr(node, "model_met_ids", set()) | {mid}
                if m not in node.model_mets:
                    node.model_mets.add(m)
                    k = len(node.model_mets) - 1
                    node.label = f"{node.base_label} (+{k})"
                    node.w = text_width(node.label)
                    mid = next(x for x in md.reactions[rid]["metabolites"] if md.base[x] == m)
                    node.title = f"{node.title}; {md.metabolites[mid]['name']}"
        self.added.append(dict(
            reaction=rid, name=md.reactions[rid]["name"], status=md.matches[rid]["status"], subsystem=enzyme.subsystem,
            model_subsystem=md.reactions[rid]["subsystem"], anchored=sorted(m for m in mains if m in chosen),
            main_metabolites=sorted(mains), new_metabolites=0, linked_metabolites=0, node=enzyme.id,
        ))


# --- hub metabolites ----------------------------------------------------------------------------------


def metabolite_key(node: Node):
    """Identity of a drawn metabolite across subsystems: BRENDA ligand ID, else the label."""
    return node.ligand_id if node.ligand_id is not None else node.label


def hub_table(dm: DraftMap) -> list[dict]:
    """Metabolites with inter-pathway links: in how many subsystems they are drawn, how many links they have."""
    subsystems = defaultdict(set)
    labels = {}
    for n in dm.nodes.values():
        if n.kind == "main":
            k = metabolite_key(n)
            subsystems[k].add(n.subsystem)
            labels.setdefault(k, n.title or n.label)
    links = Counter()
    for l in dm.links:
        a, b = dm.nodes[l.a], dm.nodes[l.b]
        if a.subsystem != b.subsystem:
            links[metabolite_key(a)] += 1
    rows = [dict(key=k, name=labels.get(k, str(k)), subsystems=len(subsystems[k]), links=v) for k, v in links.items()]
    rows.sort(key=lambda r: (-r["subsystems"], -r["links"], str(r["name"])))
    return rows


def disconnect_hubs(dm: DraftMap, min_subsystems: int, table: list[dict] | None = None) -> list[dict]:
    """Remove the inter-pathway links of metabolites drawn in at least `min_subsystems` subsystems.

    Their nodes stay in each subsystem (like currency metabolites, they are
    drawn where they are used) and are marked as hubs. `table`: the hub
    table to use (default: :func:`hub_table` of `dm`), e.g. that of the whole
    map for a part of it.
    """
    hubs = [r for r in (hub_table(dm) if table is None else table) if r["subsystems"] >= min_subsystems]
    keys = {r["key"] for r in hubs}
    before = len(dm.links)
    dm.links = [l for l in dm.links if metabolite_key(dm.nodes[l.a]) not in keys]
    for n in dm.nodes.values():
        if n.kind == "main" and metabolite_key(n) in keys:
            n.hub = True
    dm.log.append(
        f"hubs: removed {before - len(dm.links)} inter-pathway links of {len(hubs)} metabolites drawn in "
        f">= {min_subsystems} subsystems"
    )
    return hubs


def localise_edges(dm: DraftMap, subsystems: set[str] | None = None) -> int:
    """Draw the metabolites of edges between subsystems in the subsystem of the enzyme node.

    The BRENDA map has a few edges from an enzyme node to a metabolite node of
    another subsystem. If the enzyme node has the same metabolite in its own
    subsystem, the edge is removed; otherwise the metabolite gets a copy next
    to the enzyme node. Either way, the two nodes of the metabolite are linked
    by an inter-pathway link; an original that is no longer used is removed.
    `subsystems`: only the edges from or to these subsystems. Returns the
    number of edges removed or moved.
    """
    local = defaultdict(dict)  # enzyme node -> {metabolite key: metabolite node in its subsystem}
    for e in dm.edges:
        enzyme, met = dm.nodes[e.enzyme], dm.nodes[e.met]
        if enzyme.subsystem == met.subsystem:
            local[e.enzyme].setdefault(metabolite_key(met), met.id)
    copies: dict[tuple[str, str], str] = {}
    edges, links, moved = [], set(), 0
    for e in dm.edges:
        enzyme, met = dm.nodes[e.enzyme], dm.nodes[e.met]
        if enzyme.subsystem == met.subsystem or (
            subsystems is not None and enzyme.subsystem not in subsystems and met.subsystem not in subsystems
        ):
            edges.append(e)
            continue
        moved += 1
        own = local[e.enzyme].get(metabolite_key(met))
        if own is None:
            key = (met.id, enzyme.subsystem)
            if key not in copies:
                node = copy.copy(met)  # with the attributes set during completion (base_label, model_met_ids)
                node.id = f"{met.id}_in_{enzyme.subsystem}"
                node.subsystem = enzyme.subsystem
                node.model_mets = set(met.model_mets)
                boxes = [n.box for n in dm.nodes.values() if n.subsystem == enzyme.subsystem]
                node.x, node.y = place_free(enzyme.x, enzyme.y, node.w, node.h, boxes, prefer=(met.x, met.y))
                dm.nodes[node.id] = node
                copies[key] = node.id
            own = copies[key]
            e.met = own
            edges.append(e)
        links.add((own, met.id))
    dm.edges = edges
    existing = {(l.a, l.b) for l in dm.links} | {(l.b, l.a) for l in dm.links}
    dm.links.extend(MetLink(a, b, "added") for a, b in sorted(links) if (a, b) not in existing)
    used = {e.met for e in dm.edges}
    dm.remove_nodes({b for _, b in links if b not in used})
    dm.links = [l for l in dm.links if l.a in dm.nodes and l.b in dm.nodes]
    dm.log.append(f"{moved} edges between subsystems: metabolite drawn in the subsystem of the enzyme node")
    return moved
