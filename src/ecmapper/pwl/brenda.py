"""PWL 2.0 descriptions of the subsystems of a BRENDA map export (``BRENDA_map_reactions.json``).

Layout part (unchanged from PWL 1.1): the heading (``flow``) comes from the
drawn geometry, turn letters record the drawn steps, rings are cycles of the
drawing chosen by enclosed area, stars are hubs, fans (``||``) are reactions
between the same two drawn nodes.

Biochemical part (new in 2.0):

- Every compound has one token on the whole map and is declared with its
  BRENDA ligand IDs (:class:`~ecmapper.pwl.compounds.CompoundRegistry`).
- The sides (substrates, products) and the direction of a reaction come from
  the BRENDA equation of its BRENDA reaction IDs in BKMS-react and from the
  drawn cofactor arrows, see :func:`orient`.
- Group nodes ``(A + B)`` are compounds on the same side of a reaction;
  reactions whose sides are not known are written as ``react [...] {A, B, C};``.
- Motifs (cycles, substrate cycles, shunts, repeats, carrier cycles,
  transport) are detected on the result, see :mod:`ecmapper.pwl.motifs`.
"""

from __future__ import annotations

import collections
import itertools
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from .compounds import (
    CompoundInfo,
    CompoundRegistry,
    compartment_label,
    equation_signature,
    match_equation,
    split_location,
)
from .motifs import detect_motifs, format_motif
from .syntax import VERSION, parse, quote

# ---------------------------------------------------------------------------------- cofactor couples

CAND = {
    "ATP": ["ADP", "AMP"], "ADP": ["ATP", "AMP"], "AMP": ["ADP", "ATP"], "GTP": ["GDP", "GMP"], "GDP": ["GTP"],
    "CTP": ["CDP", "CMP"], "UTP": ["UDP", "UMP"], "NAD+": ["NADH"], "NADH": ["NAD+"], "NADP+": ["NADPH"],
    "NADPH": ["NADP+"], "NAD(P)+": ["NAD(P)H"], "NAD(P)H": ["NAD(P)+"], "FAD": ["FADH2"], "FADH2": ["FAD"],
    "FMN": ["FMNH2"], "FMNH2": ["FMN"], "Fd_ox": ["Fd_red"], "Fd_red": ["Fd_ox"],
    "S-adenosyl-L-methionine": ["S-adenosyl-L-homocysteine"],
    "PPi": ["Pi"], "CoA": ["Ac-CoA"], "Ac-CoA": ["CoA"], "malonyl-CoA": ["CoA"], "succinyl-CoA": ["CoA"],
    "Glu": ["2OG", "2-oxoglutarate"], "2OG": ["Glu"], "2-oxoglutarate": ["Glu"], "Gln": ["Glu"],
    "Ala": ["pyruvate"], "pyruvate": ["Ala"], "acceptor": ["reduced acceptor"],
    "reduced acceptor": ["acceptor", "oxidized acceptor"],
    "a quinone": ["reduced quinone"], "reduced quinone": ["a quinone"],
}


def partners(name: str) -> list[str]:
    c = list(CAND.get(name, []))
    if "reduced" in name:
        c.append(name.replace("reduced", "oxidized"))
    if "oxidized" in name:
        c.append(name.replace("oxidized", "reduced"))
    for pre in ("UDP-", "GDP-", "CDP-", "dTDP-", "ADP-"):
        if name.startswith(pre):
            base = pre[:-1]
            c += [base, "CMP"] if pre == "CDP-" else ([base, "AMP"] if pre == "ADP-" else [base])
    return c


# ---------------------------------------------------------------------------------- geometry helpers

HEAD = {"down": (0, 1), "up": (0, -1), "right": (1, 0), "left": (-1, 0)}


def cw(h):
    return (-h[1], h[0])


def ccw(h):
    return (h[1], -h[0])


def dir4(dx, dy):
    if abs(dx) >= abs(dy):
        return (1, 0) if dx >= 0 else (-1, 0)
    return (0, 1) if dy >= 0 else (0, -1)


def rel_letter(h, s):
    if s == h:
        return "F"
    if s == cw(h):
        return "R"
    if s == ccw(h):
        return "L"
    return "B"


def poly_area(pts):
    a = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
    return a / 2.0


def hull_area(pts):
    p = sorted(set(pts))
    if len(p) < 3:
        return 0.0

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lo, up = [], []
    for x in p:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], x) <= 0:
            lo.pop()
        lo.append(x)
    for x in reversed(p):
        while len(up) >= 2 and cross(up[-2], up[-1], x) <= 0:
            up.pop()
        up.append(x)
    return abs(poly_area(lo[:-1] + up[:-1]))


# ---------------------------------------------------------------------------------- orientation


@dataclass
class Orientation:
    """Sides and direction of a drawn reaction.

    `side` maps each main node ID to "S" (substrate side) or "P" (product
    side) in the drawn direction; it is empty if the sides are not known.
    `oriented` tells whether the direction is known; if not, "S" and "P" only
    separate the two sides. `evidence` is ``cof+eq`` (drawn cofactor arrows
    agree with the BRENDA equation), ``cof-eq`` (the cofactor arrows orient the
    reaction against the BRENDA equation; the drawn direction is used), ``eq``
    (BRENDA equation only, i.e. the direction of the natural substrates) or
    ``conflict`` (cofactor arrows contradict each other; no direction).
    """

    side: dict[str, str] = field(default_factory=dict)
    oriented: bool = False
    evidence: str | None = None
    reaction_id: str | None = None
    equation: str | None = None
    reversible: str | None = None  # evidence of reversibility: "eq" or "map"
    match: str | None = None  # how main metabolites without the same name were matched: kegg, single, similar
    coef: dict = field(default_factory=dict)  # (role, node or cofactor index) -> coefficient
    location: dict = field(default_factory=dict)  # (role, node or cofactor index) -> compartment

    @property
    def sided(self) -> bool:
        return bool(self.side)


def _participants(registry: CompoundRegistry, r: dict):
    """(role, key, compound) of the main metabolites (key: node ID) and cofactors (key: index) of a reaction."""
    out = []
    seen = set()
    for m in r["main_metabolites"]:
        if m["node_id"] in seen:
            continue
        seen.add(m["node_id"])
        out.append(("main", m["node_id"], registry.compound(m)))
    for role in ("consumed", "produced"):
        for k, m in enumerate(r["cofactors"][role]):
            out.append((role, k, registry.compound(m)))
    return out


def orient(registry: CompoundRegistry, r: dict, bkms=None) -> Orientation:
    """Sides and direction of the drawn reaction `r` from its BRENDA equations and its cofactor arrows.

    The equations of all BRENDA reaction IDs of the node are matched with the
    drawn metabolites (by name keys, see
    :func:`~ecmapper.pwl.compounds.match_equation`); the best one is used. A
    consumed cofactor on the left side or a produced one on the right side
    votes for the direction of the equation, the opposite for the reverse
    direction. Main metabolites take the side of their term; with two main
    metabolites of which only one is found, the other takes the opposite side.
    """
    parts = _participants(registry, r)
    mains = [p for p in parts if p[0] == "main"]
    if bkms is None or not r.get("brenda_reaction_ids"):
        return Orientation()
    best = None
    main_idx = [i for i, p in enumerate(parts) if p[0] == "main"]
    for br in r["brenda_reaction_ids"]:
        # IDs like BR1125-a (a part of a BRENDA reaction) are in BKMS-react without the suffix
        for _ec, equation in bkms.equations(br) or bkms.equations(re.sub(r"-[a-z]$", "", br)):
            m = match_equation(registry, [p[2] for p in parts], equation, br, main_idx)
            if m is not None and (best is None or m.score > best.score):
                best = m
    if best is None or best.score == 0:
        return Orientation()
    votes = collections.Counter()
    for i, (role, _key, _c) in enumerate(parts):
        if role == "main" or i not in best.sides:
            continue
        forward = (role == "consumed") == (best.sides[i] == "L")
        votes["+" if forward else "-"] += 1
    side = {parts[i][1]: best.sides[i] for i in main_idx if i in best.sides}
    if len(mains) == 2 and len(side) == 1:
        known = next(iter(side))
        other = next(p[1] for p in mains if p[1] != known)
        side[other] = "R" if side[known] == "L" else "L"
    o = Orientation(reaction_id=best.reaction_id, equation=best.equation)
    used = {best.methods[i] for i in main_idx if i in best.methods}
    o.match = next((m for m in ("similar", "single", "kegg") if m in used), None)
    for i, (role, key, _c) in enumerate(parts):
        if i in best.coef:
            o.coef[(role, key)] = best.coef[i]
        if i in best.location:
            o.location[(role, key)] = best.location[i]
    if len(side) != len(mains) or not mains:
        return o  # sides not known
    if votes["+"] and votes["-"]:
        o.evidence, flip = "conflict", False
    elif votes["+"]:
        o.evidence, flip, o.oriented = "cof+eq", False, True
    elif votes["-"]:
        o.evidence, flip, o.oriented = "cof-eq", True, True
    else:
        o.evidence, flip, o.oriented = "eq", False, True
    first = "R" if flip else "L"
    o.side = {node: ("S" if s == first else "P") for node, s in side.items()}
    if o.oriented and bkms is not None:
        sig = equation_signature(registry, best.equation)
        if sig is not None:
            for ec in r.get("ec_numbers", []):
                for other in bkms.ec_equations(ec):
                    if equation_signature(registry, other) == (sig[1], sig[0]):
                        o.reversible = "eq"
                        break
    return o


# ---------------------------------------------------------------------------------- per-subsystem model


class Sub:
    """A subsystem: reactions, drawn main-metabolite nodes with tokens, orientation of the reactions."""

    def __init__(self, sid: str, ctx: Context):
        self.sid = sid
        self.ctx = ctx
        self.locs: set[str] = set()
        self.meta = ctx.subs[sid]
        self.rs = ctx.by_sub[sid]
        self.rx = {r["node_id"]: r for r in self.rs}
        self.nodes = {}
        for r in self.rs:
            for m in r["main_metabolites"]:
                nd = ctx.mnode.get(m["node_id"])
                pos = (m["position"]["x"], m["position"]["y"]) if m.get("position") else (nd["x"], nd["y"])
                if m["node_id"] not in self.nodes:
                    info = ctx.registry.compound(m)
                    self.nodes[m["node_id"]] = dict(node_id=m["node_id"], compound=info, pos=pos,
                                                    loc=split_location(m.get("display_name") or "")[1])
        self.assign_tokens()
        self.compute_flow()
        self.orient = {rid: ctx.orientation(r) for rid, r in self.rx.items()}
        self.mark_reversible_pairs()

    def assign_tokens(self):
        order = sorted(self.nodes.values(), key=lambda n: (n["pos"][1], n["pos"][0], n["node_id"]))
        copies = collections.defaultdict(list)
        for n in order:
            copies[n["compound"].token].append(n)
        for token, lst in copies.items():
            for i, n in enumerate(lst):
                n["tok"] = token if i == 0 else f"{token}#{i + 1}"
        self.tok = {n["node_id"]: n["tok"] for n in self.nodes.values()}

    def compute_flow(self):
        import numpy as np

        pts = np.array([[r["layout"]["x"], r["layout"]["y"]] for r in self.rs])
        self.flow = "down"
        if len(pts) >= 3:
            w, v = np.linalg.eigh(np.cov(pts.T))
            if (w[0] > 1e-9 and w[1] / w[0] >= 1.35) or (w[0] <= 1e-9 and w[1] > 0):
                ax = v[:, 1]
                self.flow = "down" if abs(ax[1]) >= abs(ax[0]) else "right"

    def mark_reversible_pairs(self):
        """Two drawn reactions with the same EC label and main compounds in opposite directions are reversible."""
        groups = collections.defaultdict(list)
        for rid, r in self.rx.items():
            o = self.orient[rid]
            if not o.oriented:
                continue
            subs = frozenset(self.nodes[n]["compound"].cid for n, s in o.side.items() if s == "S")
            prods = frozenset(self.nodes[n]["compound"].cid for n, s in o.side.items() if s == "P")
            groups[(r["ec_label"], frozenset((subs, prods)))].append((rid, subs))
        for items in groups.values():
            if len({s for _rid, s in items}) > 1:
                for rid, _s in items:
                    self.orient[rid].reversible = self.orient[rid].reversible or "map"

    def npos(self, nid):
        return self.nodes[nid]["pos"]

    def participant(self, nid, role_key=None):
        """Text of a drawn node in a reaction: token, coefficient and compartment."""
        text = self.tok[nid]
        o = self.orient.get(role_key[0]) if role_key else None
        loc = self.nodes[nid]["loc"]
        coef = None
        if o is not None:
            loc = loc or o.location.get(("main", nid))
            coef = o.coef.get(("main", nid))
        if loc:
            self.locs.add(loc)
            text += "@" + loc
        return f"{coef}*{text}" if coef else text

    # --- reaction body
    def cofactor_items(self, rid):
        r = self.rx[rid]
        o = self.orient[rid]
        reg = self.ctx.registry

        def text(role, k, entry):
            tok = reg.token(entry)
            loc = split_location(entry.get("display_name") or "")[1] or o.location.get((role, k))
            if loc:
                self.locs.add(loc)
                tok += "@" + loc
            coef = o.coef.get((role, k))
            return f"{coef}*{tok}" if coef else tok

        cons = list(enumerate(r["cofactors"]["consumed"]))
        prod = list(enumerate(r["cofactors"]["produced"]))
        used = [False] * len(prod)
        items, single = [], []
        for k, c in cons:
            name = c["display_name"]
            hit = None
            for cand in partners(name):
                for i, (_kp, p) in enumerate(prod):
                    if not used[i] and p["display_name"] == cand:
                        hit = i
                        break
                if hit is not None:
                    break
            if hit is None:
                single.append((k, c))
            else:
                used[hit] = True
                items.append(text("consumed", k, c) + ">" + text("produced", prod[hit][0], prod[hit][1]))
        items += ["-" + text("consumed", k, c) for k, c in single]
        items += ["+" + text("produced", k, p) for i, (k, p) in enumerate(prod) if not used[i]]
        items += ["-" + reg.token(u) for u in r["cofactors"].get("undirected", [])]
        return items

    def rbody(self, rid, direction=None, side=None):
        r = self.rx[rid]
        o = self.orient[rid]
        ec = r["ec_label"].strip()
        spont = ec.lower() == "spontaneous"
        enz = "?" if (spont or ec.startswith("?")) else ec
        txt = f"{rid}: " + (direction + " " if direction else "") + ("~" if spont else "") + enz
        cofs = self.cofactor_items(rid)
        if cofs:
            txt += " | " + ", ".join(cofs)
        if side and cofs:
            txt += " ^" + side
        attrs = []
        if o.reaction_id and (o.sided or o.oriented):
            attrs.append("br=" + o.reaction_id)
        if o.evidence:
            attrs.append("ev=" + o.evidence)
        if o.reversible:
            attrs.append("rev=" + o.reversible)
        if o.match and o.sided:
            attrs.append("match=" + o.match)
        if attrs:
            txt += " {" + ", ".join(attrs) + "}"
        return txt

    def cof_side(self, r, h):
        ex, ey = r["layout"]["x"], r["layout"]["y"]
        sides = set()
        for k in ("consumed", "produced"):
            for c in r["cofactors"][k]:
                nd = self.ctx.mnode.get(c["node_id"])
                if not nd:
                    continue
                v = (nd["x"] - ex, nd["y"] - ey)
                dot = v[0] * cw(h)[0] + v[1] * cw(h)[1]
                sides.add("R" if dot > 0 else "L")
        return sides.pop() if len(sides) == 1 else None

    def comment(self, rids):
        return "  // " + "; ".join(self.rx[rid]["enzyme_name"] or "" for rid in rids)

    # --- directions of reactions between two drawn nodes
    def arrow(self, rid, a, b):
        """Direction of reaction `rid` written from node a to node b."""
        o = self.orient[rid]
        if not o.oriented:
            return None
        if o.reversible:
            return "<->"
        if o.side.get(a) == "S" and o.side.get(b) == "P":
            return "->"
        if o.side.get(a) == "P" and o.side.get(b) == "S":
            return "<-"
        return None

    def edge_dir(self, rids, a, b) -> int:
        arrows = {self.arrow(rid, a, b) for rid in rids}
        if arrows == {"->"}:
            return 1
        if arrows == {"<-"}:
            return -1
        return 0


_OPEN = {None: "-", "->": "-", "<-": "<-", "<->": "<-", "=>": "-", "<=": "<="}
_CLOSE = {None: "-", "->": "->", "<-": "-", "<->": "->", "=>": "=>", "<=": "-"}


def build(sid: str, ctx: Context) -> tuple[str, dict, dict]:
    """PWL text of a subsystem, its statistics and its node index {token: drawn node ID}."""
    import networkx as nx
    import numpy as np

    S = Sub(sid, ctx)
    out = []
    emitted = []
    stats = collections.Counter()
    stats["reactions"] = len(S.rs)
    stats["nodes"] = len(S.nodes)

    # --- reactions between two drawn nodes on opposite sides (pairs) and the others (hyper)
    pair = collections.defaultdict(list)
    hyper = []
    for r in S.rs:
        ids = list(dict.fromkeys(m["node_id"] for m in r["main_metabolites"]))
        o = S.orient[r["node_id"]]
        if len(ids) == 2 and not (o.sided and o.side[ids[0]] == o.side[ids[1]]):
            pair[frozenset(ids)].append(r["node_id"])
        else:
            hyper.append(r["node_id"])
        if o.oriented:
            stats["oriented"] += 1
            stats["evidence_" + o.evidence] += 1
        elif o.sided:
            stats["sided_unoriented"] += 1
        if o.reversible:
            stats["reversible"] += 1
        if len(ids) == 2 and o.sided and o.side[ids[0]] == o.side[ids[1]]:
            stats["pair_same_side"] += 1

    def sides_for(rid):
        """Two sides of a reaction with several main metabolites: from the orientation, else geometric."""
        r = S.rx[rid]
        o = S.orient[rid]
        ids = list(dict.fromkeys(m["node_id"] for m in r["main_metabolites"]))
        if o.sided:
            return [i for i in ids if o.side[i] == "S"], [i for i in ids if o.side[i] == "P"], True
        ex, ey = r["layout"]["x"], r["layout"]["y"]
        P = np.array([S.npos(i) for i in ids])
        sp = P.max(axis=0) - P.min(axis=0)
        axis = 0 if sp[0] >= sp[1] else 1
        center = (ex, ey)[axis]
        left = [i for i in ids if S.npos(i)[axis] < center - 1e-6]
        right = [i for i in ids if i not in left]
        if not left or not right:
            srt = sorted(ids, key=lambda i: (S.npos(i)[axis], S.npos(i)[1 - axis]))
            k = max(1, len(srt) // 2)
            left, right = srt[:k], srt[k:]
        return left, right, False

    hyper1 = {h for h in hyper if len(set(m["node_id"] for m in S.rx[h]["main_metabolites"])) == 1}
    hyper3 = [h for h in hyper if h not in hyper1]
    G = nx.Graph()
    for k in pair:
        G.add_edge(*sorted(k))
    remaining = set(pair.keys())
    derived = collections.defaultdict(list)
    for rid in hyper3:
        left, right, _known = sides_for(rid)
        for a in left:
            for b in right:
                e = frozenset((a, b))
                if e not in pair:
                    derived[e].append(rid)
                    G.add_edge(a, b)
    head0 = HEAD[S.flow]

    def edge_rids(a, b):
        return sorted(pair[frozenset((a, b))], key=lambda x: (S.rx[x]["layout"]["y"], S.rx[x]["layout"]["x"]))

    def link_text(a, b, h, with_side=True):
        pa, pb = S.npos(a), S.npos(b)
        s = dir4(pb[0] - pa[0], pb[1] - pa[1]) if (abs(pb[0] - pa[0]) + abs(pb[1] - pa[1])) > 1e-6 else h
        letter = rel_letter(h, s)
        rids = edge_rids(a, b)
        arrows = [S.arrow(rid, a, b) for rid in rids]
        common = arrows[0] if len(set(arrows)) == 1 else None
        bodies = []
        for rid, arr in zip(rids, arrows):
            sd = S.cof_side(S.rx[rid], s) if (with_side and len(rids) == 1) else None
            bodies.append(S.rbody(rid, None if len(set(arrows)) == 1 else arr, sd))
        t = _OPEN[common] + ("" if letter == "F" else letter) + "[" + " || ".join(bodies) + "]" + _CLOSE[common]
        emitted.extend(rids)
        if len(rids) > 1:
            stats["fan"] += 1
        return t, s, rids

    def node_text(nids, rid=None):
        toks = [S.participant(i, (rid,) if rid else None) for i in nids]
        if not toks:
            return "_"
        return toks[0] if len(toks) == 1 else "(" + " + ".join(toks) + ")"

    def hyper_stmt(rid, prefix=""):
        """A reaction with one or more than two main metabolites, or two on the same side."""
        r = S.rx[rid]
        o = S.orient[rid]
        ids = list(dict.fromkeys(m["node_id"] for m in r["main_metabolites"]))
        emitted.append(rid)
        if o.sided:
            subs = sorted([i for i in ids if o.side[i] == "S"], key=lambda i: (S.npos(i)[1], S.npos(i)[0]))
            prods = sorted([i for i in ids if o.side[i] == "P"], key=lambda i: (S.npos(i)[1], S.npos(i)[0]))
            arrow = ("<->" if o.reversible else "->") if o.oriented else None
            stats["group_sided"] += 1
            return "%s%s %s[%s]%s %s;%s" % (prefix, node_text(subs, rid), _OPEN[arrow], S.rbody(rid), _CLOSE[arrow],
                                            node_text(prods, rid), S.comment([rid]))
        if len(ids) == 1:
            stats["open_end"] += 1
            return "%s%s -[%s]- _;%s" % (prefix, node_text(ids, rid), S.rbody(rid), S.comment([rid]))
        stats["unsided"] += 1
        parts = ", ".join(S.participant(i, (rid,)) for i in sorted(ids, key=lambda i: (S.npos(i)[1], S.npos(i)[0])))
        return "%sreact [%s] {%s};%s" % ("", S.rbody(rid), parts, S.comment([rid]))

    # --- ports
    ports = []
    seen_ports = set()
    for link in ctx.data["metabolite_links"]:
        for me, other in (("start", "end"), ("end", "start")):
            if link[me + "_subsystem"] == sid and link[other + "_subsystem"] != sid:
                nid = link[me + "_node_id"]
                if nid not in S.nodes:
                    continue
                pa = S.npos(nid)
                opos = link[other]
                side = {(1, 0): "R", (-1, 0): "L", (0, 1): "D", (0, -1): "U"}[dir4(opos[0] - pa[0], opos[1] - pa[1])]
                key = (nid, side, link[other + "_subsystem"])
                if key in seen_ports:
                    continue
                seen_ports.add(key)
                ports.append((S.tok[nid], side, link[other + "_subsystem"]))
    ports.sort()
    stats["ports"] = len(ports)

    # --- rings (cycles of the drawing), chosen by enclosed area
    ring_blocks = []

    def pick_ring(Gc):
        best = None
        for comp in nx.biconnected_components(Gc):
            if len(comp) < 4:
                continue
            H = Gc.subgraph(comp)
            if H.number_of_nodes() > 70:
                continue
            for cyc in itertools.islice(nx.simple_cycles(H, length_bound=24), 30000):
                if len(cyc) < 4:
                    continue
                pts = [S.npos(n) for n in cyc]
                ar = abs(poly_area(pts))
                ha = hull_area(pts)
                if ha <= 0 or ar < 40000 or ar / ha < 0.55:
                    continue
                sc = ar * (len(cyc) ** 0.5)
                if best is None or sc > best[0]:
                    best = (sc, cyc)
        return best

    Gr = G.copy()
    while True:
        b = pick_ring(Gr)
        if not b:
            break
        cyc = b[1]
        if poly_area([S.npos(n) for n in cyc]) < 0:
            cyc = cyc[::-1]
        i0 = min(range(len(cyc)), key=lambda i: (S.npos(cyc[i])[1], S.npos(cyc[i])[0]))
        cyc = cyc[i0:] + cyc[:i0]
        ring_blocks.append(cyc)
        used_h = set()
        for i in range(len(cyc)):
            e = frozenset((cyc[i], cyc[(i + 1) % len(cyc)]))
            if e not in pair:
                used_h.add(derived[e][0])
            Gr.remove_edge(*tuple(e))
            remaining.discard(e)
        for e2, rl in list(derived.items()):
            if any(x in used_h for x in rl):
                derived[e2] = [x for x in rl if x not in used_h] or rl
                if Gr.has_edge(*tuple(e2)) and e2 not in pair and all(x in used_h for x in rl):
                    Gr.remove_edge(*tuple(e2))
        Gr.remove_nodes_from([n for n in list(Gr.nodes) if Gr.degree(n) == 0])
        if len(ring_blocks) >= 4:
            break

    ring_infos = []
    consumed_hyper = set()
    for cyc in ring_blocks:
        cycset = set(cyc)
        n_ = len(cyc)
        # traverse in the direction of the reactions if more of them point against the drawing order
        fwd = bwd = 0
        for i in range(n_):
            e = frozenset((cyc[i], cyc[(i + 1) % n_]))
            if e in pair:
                d = S.edge_dir(pair[e], cyc[i], cyc[(i + 1) % n_])
                fwd += d > 0
                bwd += d < 0
        sense = "cw"
        if bwd > fwd:
            cyc = [cyc[0]] + cyc[1:][::-1]
            sense = "ccw"
        etype = []
        for i in range(n_):
            e = frozenset((cyc[i], cyc[(i + 1) % n_]))
            etype.append(("N", None) if e in pair else ("H", derived[e][0]))
        chords = [e for e in sorted(remaining, key=sorted) if e <= cycset]
        for e in chords:
            remaining.discard(e)
        closure = []
        for t_, x_ in etype:
            if t_ == "H" and x_ not in closure:
                closure.append(x_)
        hch = []
        for rid in hyper3:
            if rid in closure or rid in consumed_hyper:
                continue
            left, right, _k = sides_for(rid)
            if any(a in cycset for a in left) and any(b in cycset for b in right):
                hch.append(rid)
        consumed_hyper |= set(closure) | set(hch)
        ring_infos.append((cyc, sense, chords, etype, closure, hch))

    # --- stars (hubs with 5 or more remaining neighbours)
    def rem_adj():
        # sorted, so that the output does not depend on the order of the set (hash randomization)
        adj = collections.defaultdict(list)
        for e in sorted(remaining, key=sorted):
            a, b = sorted(e)
            adj[a].append(b)
            adj[b].append(a)
        return adj

    star_blocks = []
    while True:
        adj = rem_adj()
        hubs = [(len(v), n) for n, v in adj.items() if len(v) >= 5]
        if not hubs:
            break
        hubs.sort(key=lambda t: (-t[0], t[1]))
        hub = hubs[0][1]
        hp = S.npos(hub)
        spokes = []
        for nb in sorted(adj[hub], key=lambda n: math.atan2(S.npos(n)[1] - hp[1], S.npos(n)[0] - hp[0])):
            path = [hub, nb]
            remaining.discard(frozenset((hub, nb)))
            cur = nb
            while True:
                adj2 = rem_adj()
                nxt = [x for x in adj2.get(cur, []) if x != path[-2] and len(adj2.get(cur, [])) == 1]
                if not nxt or nxt[0] == hub or len(adj2.get(nxt[0], [])) >= 5:
                    break
                remaining.discard(frozenset((cur, nxt[0])))
                path.append(nxt[0])
                cur = nxt[0]
            spokes.append(path)
        star_blocks.append((hub, spokes))

    # --- path cover of the rest, following the direction of the reactions where it is known
    placed = set()
    for info in ring_infos:
        placed |= set(info[0])
    for _hub, sp in star_blocks:
        for p in sp:
            placed |= set(p)
    axis = 1 if S.flow == "down" else 0

    def incoming(adj, n):
        return any(S.edge_dir(pair[frozenset((m, n))], m, n) > 0 for m in adj[n])

    def pick_start(adj):
        cand = [n for n in adj if adj[n]]
        return min(cand, key=lambda n: (n not in placed, incoming(adj, n), len(adj[n]) % 2 == 0,
                                        S.npos(n)[axis], S.npos(n)[1 - axis], n))

    chains = []
    while remaining:
        adj = rem_adj()
        start = pick_start(adj)
        path = [start]
        h = head0
        cur = start
        while True:
            adj = rem_adj()
            cand = adj.get(cur, [])
            if not cand:
                break
            pc = S.npos(cur)

            def score(n, pc=pc, h=h, cur=cur):
                pn = S.npos(n)
                v = (pn[0] - pc[0], pn[1] - pc[1])
                ln = math.hypot(*v) or 1
                d = S.edge_dir(pair[frozenset((cur, n))], cur, n)
                return (d >= 0, d > 0, (v[0] * h[0] + v[1] * h[1]) / ln, -pn[1] if S.flow == "down" else -pn[0],
                        n)

            nxt = max(cand, key=score)
            remaining.discard(frozenset((cur, nxt)))
            pn = S.npos(nxt)
            v = (pn[0] - pc[0], pn[1] - pc[1])
            if abs(v[0]) + abs(v[1]) > 1e-6:
                h = dir4(*v)
            path.append(nxt)
            cur = nxt
        placed |= set(path)
        # write a chain in the direction of most of its reactions
        dirs = [S.edge_dir(pair[frozenset((a, b))], a, b) for a, b in zip(path, path[1:])]
        if sum(d < 0 for d in dirs) > sum(d > 0 for d in dirs):
            path = path[::-1]
        chains.append(path)

    # --- header
    reg = ctx.registry
    out.append(f"pwl {VERSION};")
    out.append(f"subsystem {sid} {quote(S.meta['label_text'])} {{")
    n_or = stats["oriented"]
    out.append(f"  // BRENDA subsystem: {len(S.rs)} reactions ({n_or} with direction), "
               f"{len(S.nodes)} drawn main-metabolite nodes.")
    out.append("  // flow, turn letters, ring, star and wrap describe the drawing; "
               "arrows and groups describe the reactions.")
    out.append(f"  flow {S.flow};")
    body = []

    # rings
    for cyc, sense, chords, etype, closure, hch in ring_infos:
        stats["ring"] += 1
        stats["chord"] += len(chords) + len(hch)
        n = len(cyc)
        body.append(f"  // ring: cycle of the drawing ({n} metabolite nodes); see the motifs for biochemical cycles")
        body.append(f"  ring {sense} {{")
        hidx = [i for i, (t_, _x) in enumerate(etype) if t_ == "H"]
        start = (hidx[0] + 1) % n if hidx else 0
        order = [(start + i) % n for i in range(n)]
        runs, cur = [], []
        for i in order:
            if etype[i][0] == "N":
                cur.append(i)
            elif cur:
                runs.append(cur)
                cur = []
        if cur:
            runs.append(cur)
        for run in runs:
            body.append("    " + S.tok[cyc[run[0]]])
            h = head0
            for j, i in enumerate(run):
                a_, b_ = cyc[i], cyc[(i + 1) % n]
                t, h, rids = link_text(a_, b_, h)
                body.append("      %s %s%s%s" % (t, S.tok[b_], ";" if j == len(run) - 1 else "", S.comment(rids)))
        for rid in closure:
            body.append("    " + hyper_stmt(rid))
        for e in sorted(chords, key=lambda e: tuple(sorted(S.tok[x] for x in e))):
            a_, b_ = sorted(e, key=lambda x: (S.npos(x)[1], S.npos(x)[0]))
            if S.edge_dir(pair[e], a_, b_) < 0:
                a_, b_ = b_, a_
            t, _hh, rids = link_text(a_, b_, head0)
            body.append("    chord %s %s %s;%s" % (S.tok[a_], t, S.tok[b_], S.comment(rids)))
        for rid in hch:
            stmt = hyper_stmt(rid, "chord ")
            body.append("    " + stmt)
        body.append("  }")
        body.append("")

    # stars
    for hub, spokes in star_blocks:
        stats["star"] += 1
        hp = S.npos(hub)
        angs = sorted(math.degrees(math.atan2(S.npos(p[1])[1] - hp[1], S.npos(p[1])[0] - hp[0])) % 360 for p in spokes)
        gaps = [(angs[(i + 1) % len(angs)] - angs[i]) % 360 for i in range(len(angs))]
        arc = int(round((360 - max(gaps)) / 30.0) * 30) or 30
        body.append(f"  star {S.tok[hub]} arc {min(arc, 360)} order cw {{")
        for p in spokes:
            h = head0
            parts = []
            for i in range(1, len(p)):
                t, h, rids = link_text(p[i - 1], p[i], h)
                parts.append((t, S.tok[p[i]], rids))
            line = "    " + parts[0][0] + " " + parts[0][1]
            if len(parts) == 1:
                body.append(line + ";" + S.comment(parts[0][2]))
            else:
                body.append(line + S.comment(parts[0][2]))
                for j, (t, nd, rids) in enumerate(parts[1:], 1):
                    body.append("      %s %s%s%s" % (t, nd, ";" if j == len(parts) - 1 else "", S.comment(rids)))
        body.append("  }")
        body.append("")

    # chains
    for path in chains:
        stats["chain"] += 1
        h = head0
        body.append("  " + S.tok[path[0]])
        rows = []
        for i in range(1, len(path)):
            t, h, rids = link_text(path[i - 1], path[i], h)
            rows.append((t, S.tok[path[i]], rids))
        seq = ["|".join(sorted(S.rx[x]["ec_label"] for x in rids)) for _t, _n, rids in rows]
        wrap = None
        L = len(seq)
        if L >= 9:
            best = None
            for p in range(2, 9):
                run = bl = 0
                for i in range(L - p):
                    if seq[i] == seq[i + p]:
                        run += 1
                        bl = max(bl, run)
                    else:
                        run = 0
                tot = bl + p
                if tot >= 3 * p and tot >= 0.8 * L and (best is None or tot > best[1]):
                    best = (p, tot)
            if best:
                wrap = best[0]
                stats["wrap"] += 1
        for j, (t, nd, rids) in enumerate(rows):
            end = ((f" wrap serpentine period {wrap}" if wrap else "") + ";") if j == len(rows) - 1 else ""
            body.append("    %s %s%s%s" % (t, nd, end, S.comment(rids)))
        body.append("")

    # reactions with one, more than two, or two main metabolites on the same side
    for rid in sorted(hyper, key=lambda x: (S.rx[x]["layout"]["y"], S.rx[x]["layout"]["x"])):
        if rid in consumed_hyper:
            continue
        body.append("  " + hyper_stmt(rid))

    assert sorted(emitted) == sorted(S.rx), (sid, len(emitted), len(S.rx))

    # declarations (after the body: compartments are collected while writing it)
    used = sorted({n["compound"].token for n in S.nodes.values()}
                  | {reg.token(c) for r in S.rs for k in ("consumed", "produced", "undirected")
                     for c in r["cofactors"].get(k, [])})
    for loc in sorted(S.locs):
        out.append(f"  compartment {loc} {quote(compartment_label(loc))};")
    for token in used:
        info = reg.by_token(token)
        xrefs = " ".join(f"brenda:{lig}" for lig in sorted(info.ligand_ids))
        xrefs += "".join(f" kegg:{x}" for x in sorted(reg.xref_ids(info)))
        label = "" if info.label == token else " " + quote(info.label)
        out.append(f"  compound {token}{label}{(' ' + xrefs) if xrefs.strip() else ''};")
    for t, sd, tgt in ports:
        out.append(f"  port {t} side {sd} to {tgt};")
    out.append("")
    out.extend(body)
    text = "\n".join(out).rstrip() + "\n}\n"

    # motifs, detected on the parsed description
    doc = parse(text)
    motifs = detect_motifs(doc)
    if motifs:
        lines = text.rstrip().split("\n")[:-1]
        lines.append("  // motifs (detected from the reactions above, see PWL.md)")
        lines += ["  " + format_motif(m) for m in motifs]
        text = "\n".join(lines) + "\n}\n"
        for m in motifs:
            stats["motif_" + m.kind] += 1
    stats["compounds_declared"] = len(used)
    index = {S.tok[nid]: nid for nid in S.nodes}
    return text, dict(stats), index


# ---------------------------------------------------------------------------------- whole map


class Context:
    """The export, the compound registry and the orientation of all reactions."""

    def __init__(self, data: dict, bkms=None, kegg_names: dict | None = None):
        self.data = data
        self.subs = {s["id"]: s for s in data["subsystems"]}
        self.by_sub = collections.defaultdict(list)
        for r in data["reactions"]:
            self.by_sub[r["subsystem"]].append(r)
        self.mnode = {m["node_id"]: m for m in data["metabolite_nodes"]}
        self.registry = CompoundRegistry(data, kegg_names)
        self.bkms = bkms
        self._orient: dict[str, Orientation] = {}

    def orientation(self, r: dict) -> Orientation:
        if r["node_id"] not in self._orient:
            self._orient[r["node_id"]] = orient(self.registry, r, self.bkms)
        return self._orient[r["node_id"]]


def generate(data: dict, outdir, bkms=None, kegg_names: dict | None = None, subsystems=None) -> dict:
    """Write ``subsystems/<id>.pwl``, ``index.json``, ``node_index.json`` and ``compounds.json`` to `outdir`.

    `bkms` is a :class:`~ecmapper.sources.BKMSSource` for the BRENDA equations;
    without it no reaction gets sides or a direction. `kegg_names` (the data
    of :class:`~ecmapper.sources.KEGGCompoundNames`) adds synonyms for the
    matching and KEGG compound IDs. Returns the summary.
    """
    outdir = Path(outdir)
    (outdir / "subsystems").mkdir(parents=True, exist_ok=True)
    ctx = Context(data, bkms, kegg_names)
    index, node_index = [], {}
    total = collections.Counter()
    for sid in sorted(subsystems or [s for s in ctx.subs if ctx.by_sub[s]]):
        text, st, idx = build(sid, ctx)
        path = outdir / "subsystems" / f"{sid}.pwl"
        path.write_text(text, encoding="utf-8")
        node_index[sid] = idx
        index.append(dict(id=sid, title=ctx.subs[sid]["label_text"], file=f"subsystems/{sid}.pwl",
                          bytes=len(text.encode()), **st))
        total.update(st)
    summary = dict(pwl_version=VERSION, source=data.get("source"), subsystem_count=len(index),
                   compounds=len(ctx.registry), totals=dict(sorted(total.items())), subsystems=index)
    _dump(outdir / "index.json", summary)
    _dump(outdir / "node_index.json", node_index)
    _dump(outdir / "compounds.json", ctx.registry.to_json())
    return summary


def _dump(path: Path, obj) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False)
        f.write("\n")


__all__ = ["CompoundInfo", "Context", "Orientation", "build", "generate", "orient"]
