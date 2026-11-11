"""Detection of metabolic motifs in a parsed PWL document.

The detection works on the reactions of a :class:`~ecmapper.pwl.syntax.Document`
(their sides, directions and cofactors), not on the drawing, so it applies to
hand-written descriptions as well. Main metabolites are compared by compound
token (drawn copies ``X#2`` are the same compound); cofactors are not part
of the graph of main metabolites.

Kinds (members and roles as in :data:`~ecmapper.pwl.syntax.MOTIF_KINDS`):

``cycle``
    A directed cycle of at least three main metabolites; every step is a
    reaction that runs in this direction (oriented or reversible) and at most
    half of the steps are only reversible reactions. Members: the compounds in
    the order of the cycle.
``substrate_cycle``
    Different enzymes convert the same main metabolites in opposite
    directions (phosphofructokinase and fructose-bisphosphatase); none of
    them is known to be reversible. Members: the compounds; roles ``forward``
    and ``backward``: the reactions.
``shunt``
    A route of at most three steps from u to v that bypasses another route
    from u to v which is at least two steps longer and shares no metabolite
    with it in between (glyoxylate shunt). Short routes whose steps are all
    converted back by other reactions are left out. Members: the short route;
    role ``bypasses``: the long route.
``repeat``
    A chain whose EC numbers repeat with a period of 2 to 8 steps, at least
    three times (fatty acid elongation); the reactions of the chain with a
    direction all point the same way. Members: one reaction per step, in
    the direction of the reactions; parameters ``period`` and ``count``.
``carrier``
    A cofactor couple that is turned over in both directions within the
    subsystem (NAD+ > NADH in one reaction, NADH > NAD+ in another).
    Members: the two forms; roles ``forward`` and ``backward``: the reactions.
``transport``
    A reaction with participants in different compartments. Members: the
    compounds with compartments; role ``via``: the reaction.
"""

from __future__ import annotations

import collections

from .syntax import Chain, Document, Motif, Reaction, Ring, Star

MAX_CYCLES = 3
MAX_SHUNTS = 3


def cofactor_tokens(doc: Document) -> set[str]:
    """Compounds that are cofactors of some reaction of the document (e.g. acetyl-CoA, NH3)."""
    return {p.token for r in doc.reactions() for p in r.consumed + r.produced}


def _steps(doc: Document, exclude: set[str] = frozenset()) -> dict[tuple[str, str], list[tuple[Reaction, bool]]]:
    """{(a, b): [(reaction, one_way)]}: reactions that convert main metabolite a into b.

    Compounds in `exclude` (cofactors of the document) are left out, so that
    cycles and shunts do not run through acetyl-CoA, NH3 or ATP.
    """
    steps = collections.defaultdict(list)
    for r in doc.reactions():
        subs, prods = r.substrates, r.products
        if subs is None:
            continue
        for a in {p.token for p in subs} - exclude:
            for b in {p.token for p in prods} - exclude:
                if a == b:
                    continue
                steps[(a, b)].append((r, not r.reversible))
                if r.reversible:
                    steps[(b, a)].append((r, False))
    return steps


def _graph(steps):
    succ = collections.defaultdict(set)
    for a, b in steps:
        succ[a].add(b)
    return succ


def _simple_cycles(succ, max_len: int, limit: int = 20000):
    """Elementary directed cycles of at most `max_len` nodes, each once, starting at its smallest node."""
    nodes = sorted(set(succ) | {b for v in succ.values() for b in v})
    count = 0
    for start in nodes:
        stack = [(start, [start], {start})]
        while stack:
            node, path, seen = stack.pop()
            for nxt in sorted(succ.get(node, ()), reverse=True):
                if nxt == start and len(path) >= 3:
                    yield list(path)
                    count += 1
                    if count >= limit:
                        return
                elif nxt not in seen and nxt > start and len(path) < max_len:
                    stack.append((nxt, path + [nxt], seen | {nxt}))


def find_cycles(doc: Document, max_len: int = 16) -> list[list[str]]:
    steps = _steps(doc, cofactor_tokens(doc))
    one_way = {k for k, v in steps.items() if any(ow for _r, ow in v)}
    candidates = []
    for cyc in _simple_cycles(_graph(steps), max_len):
        edges = [(cyc[i], cyc[(i + 1) % len(cyc)]) for i in range(len(cyc))]
        if sum(e in one_way for e in edges) * 2 < len(edges):
            continue
        candidates.append(cyc)
    candidates.sort(key=lambda c: (-len(c), c))
    chosen, covered = [], set()
    for cyc in candidates:
        edges = {(cyc[i], cyc[(i + 1) % len(cyc)]) for i in range(len(cyc))}
        if len(edges - covered) >= 2:
            chosen.append(cyc)
            covered |= edges
        if len(chosen) >= MAX_CYCLES:
            break
    return chosen


def find_substrate_cycles(doc: Document) -> list[Motif]:
    groups = collections.defaultdict(lambda: ([], []))
    for r in doc.reactions():
        if r.substrates is None or r.reversible or not r.substrates or not r.products:
            continue
        a = frozenset(p.token for p in r.substrates)
        b = frozenset(p.token for p in r.products)
        if a == b:
            continue
        key = tuple(sorted((tuple(sorted(a)), tuple(sorted(b)))))
        groups[key][0 if tuple(sorted(a)) == key[0] else 1].append(r)
    out = []
    for key, (fwd, bwd) in sorted(groups.items()):
        if not fwd or not bwd:
            continue
        # the same enzyme in both directions is a reversible reaction, not a substrate cycle
        ecs_f = {r.enzyme_label for r in fwd}
        ecs_b = {r.enzyme_label for r in bwd}
        fwd = [r for r in fwd if r.enzyme_label not in ecs_b and r.enzyme_label != "?"]
        bwd = [r for r in bwd if r.enzyme_label not in ecs_f and r.enzyme_label != "?"]
        if fwd and bwd:
            out.append(Motif("substrate_cycle", list(key[0]) + list(key[1]),
                             roles={"forward": [r.rid for r in fwd], "backward": [r.rid for r in bwd]}))
    return out


def _shortest(succ, src, dst, banned_nodes, banned_edges, max_len):
    prev = {src: None}
    frontier = [src]
    depth = 0
    while frontier and depth < max_len:
        depth += 1
        nxt = []
        for u in frontier:
            for v in sorted(succ.get(u, ())):
                if (u, v) in banned_edges or v in prev or (v in banned_nodes and v != dst):
                    continue
                prev[v] = u
                if v == dst:
                    path = [v]
                    while prev[path[-1]] is not None:
                        path.append(prev[path[-1]])
                    return path[::-1]
                nxt.append(v)
        frontier = nxt
    return None


def find_shunts(doc: Document, max_short: int = 3, max_long: int = 10) -> list[tuple[list[str], list[str]]]:
    steps = _steps(doc, cofactor_tokens(doc))
    succ = _graph(steps)
    found = []
    for u in sorted(succ):
        if len(succ[u]) < 2:
            continue
        reach = {}
        frontier, depth = [u], 0
        while frontier and depth < max_short:
            depth += 1
            frontier = [v for x in frontier for v in succ.get(x, ()) if v not in reach and v != u]
            for v in frontier:
                reach.setdefault(v, depth)
        for v in sorted(reach):
            short = _shortest(succ, u, v, set(), set(), max_short)
            if short is None:
                continue
            inner = set(short[1:-1])
            edges = set(zip(short, short[1:]))
            long = _shortest(succ, u, v, inner, edges, max_long)
            if long is None or len(long) < len(short) + 2:
                continue
            # the long route must not just walk back through the short one
            if set(long[1:-1]) & set(short):
                continue
            # a short route whose steps are converted back as well is a substrate cycle or reversible
            if all((b, a) in steps for a, b in zip(short, short[1:])):
                continue
            found.append((short, long))
    # keep the candidates with the longest bypassed routes, not two for the same pair of routes' nodes
    found.sort(key=lambda sl: (-len(sl[1]), len(sl[0]), sl))
    chosen, seen = [], []
    for short, long in found:
        nodes = set(short) | set(long)
        if any(nodes <= s for s in seen):
            continue
        chosen.append((short, long))
        seen.append(nodes)
        if len(chosen) >= MAX_SHUNTS:
            break
    return chosen


def find_repeats(doc: Document, min_count: int = 3) -> list[Motif]:
    out = []
    chains = []
    for s in doc.statements:
        if isinstance(s, Chain):
            chains.append(s)
        elif isinstance(s, Ring):
            chains += [c for c in s.statements if isinstance(c, Chain)]
        elif isinstance(s, Star):
            chains += s.spokes
    for chain in chains:
        if len(chain.links) < 2 * min_count:
            continue
        # the chain must not contradict itself: reactions without direction are allowed
        dirs = {r.direction for link in chain.links for r in link.reactions} - {None, "<->"}
        if dirs and dirs <= {"<-", "<="}:
            links = chain.links[::-1]
        elif not dirs or dirs <= {"->", "=>"}:
            links = chain.links
        else:
            continue
        seq = ["|".join(sorted(r.enzyme_label for r in link.reactions)) for link in links]
        best = None
        for p in range(2, 9):
            for start in range(len(seq)):
                n = 0
                unit = seq[start: start + p]
                while start + (n + 1) * p <= len(seq) and seq[start + n * p: start + (n + 1) * p] == unit:
                    n += 1
                if n >= min_count and (best is None or n * p > best[1] * best[0]):
                    best = (p, n, start)
        if best:
            p, n, start = best
            rids = [links[i].reactions[0].rid for i in range(start, start + n * p)]
            out.append(Motif("repeat", rids, params={"period": str(p), "count": str(n)}))
    return out


def find_carriers(doc: Document) -> list[Motif]:
    couples = collections.defaultdict(list)
    for r in doc.reactions():
        for c in r.cofactors:
            if c.kind == "couple" and c.a.token != c.b.token:
                couples[(c.a.token, c.b.token)].append(r.rid)
    out = []
    for (a, b), rids in sorted(couples.items()):
        if a < b and (b, a) in couples:
            out.append(Motif("carrier", [a, b], roles={"forward": sorted(set(rids)),
                                                        "backward": sorted(set(couples[(b, a)]))}))
    return out


def find_transport(doc: Document) -> list[Motif]:
    out = []
    for r in doc.reactions():
        located = [p for p in r.participants + r.consumed + r.produced if p.location]
        if len({p.location for p in located}) >= 2:
            members = sorted({p.species for p in located})
            out.append(Motif("transport", members, roles={"via": [r.rid]}))
    return out


def detect_motifs(doc: Document) -> list[Motif]:
    """All motifs of a document, see the module docstring."""
    motifs = [Motif("cycle", c) for c in find_cycles(doc)]
    motifs += find_substrate_cycles(doc)
    motifs += [Motif("shunt", s, roles={"bypasses": lng}) for s, lng in find_shunts(doc)]
    motifs += find_repeats(doc)
    motifs += find_carriers(doc)
    motifs += find_transport(doc)
    return motifs


def format_motif(m: Motif) -> str:
    params = "".join(f" {k}={v}" for k, v in m.params.items())
    roles = "".join(f" {role} {{ {' '.join(members)} }}" for role, members in m.roles.items())
    return f"motif {m.kind}{params} {{ {' '.join(m.members)} }}{roles};"




# ---------------------------------------------------------------------------------- checks


def check_motif(doc: Document, motif: Motif) -> list[str]:
    """Problems of a motif of a known kind: does the document support it? (Works for hand-written motifs.)"""
    where = f"line {motif.line}: motif {motif.kind}"
    reactions = doc.reaction_index()
    steps = _steps(doc)
    problems = []

    def route(nodes, closed=False):
        pairs = list(zip(nodes, nodes[1:])) + ([(nodes[-1], nodes[0])] if closed else [])
        for a, b in pairs:
            if (a, b) not in steps:
                problems.append(f"{where}: no reaction converts {a} into {b}")

    if motif.kind == "cycle":
        if len(motif.members) < 3:
            problems.append(f"{where}: a cycle has at least three members")
        route(motif.members, closed=True)
    elif motif.kind == "substrate_cycle":
        fwd = [reactions[r] for r in motif.roles.get("forward", []) if r in reactions]
        bwd = [reactions[r] for r in motif.roles.get("backward", []) if r in reactions]
        if not fwd or not bwd:
            problems.append(f"{where}: needs forward and backward reactions")
        sig = {(frozenset(p.token for p in r.substrates), frozenset(p.token for p in r.products))
               for r in fwd if r.substrates is not None}
        back = {(frozenset(p.token for p in r.products), frozenset(p.token for p in r.substrates))
                for r in bwd if r.substrates is not None}
        if len(sig) != 1 or sig != back:
            problems.append(f"{where}: forward and backward reactions do not convert the same compounds "
                            "in opposite directions")
        if any(r.reversible for r in fwd + bwd):
            problems.append(f"{where}: contains a reversible reaction")
    elif motif.kind == "shunt":
        long = motif.roles.get("bypasses", [])
        short = motif.members
        if len(short) < 2 or len(long) < len(short) + 2:
            problems.append(f"{where}: the bypassed route must be at least two steps longer")
        elif short[0] != long[0] or short[-1] != long[-1]:
            problems.append(f"{where}: routes do not have the same ends")
        elif set(short[1:-1]) & set(long[1:-1]):
            problems.append(f"{where}: routes share metabolites between their ends")
        route(short)
        route(long)
    elif motif.kind == "repeat":
        try:
            period, count = int(motif.params["period"]), int(motif.params["count"])
        except (KeyError, ValueError):
            return [f"{where}: needs integer parameters period and count"]
        rs = [reactions[r] for r in motif.members if r in reactions]
        if len(rs) != period * count:
            problems.append(f"{where}: {len(rs)} reactions, expected period*count = {period * count}")
        elif any(rs[i].enzyme_label != rs[i + period].enzyme_label for i in range(len(rs) - period)):
            problems.append(f"{where}: EC numbers do not repeat with period {period}")
    elif motif.kind == "carrier":
        if len(motif.members) != 2:
            return [f"{where}: a carrier has two forms"]
        a, b = motif.members
        for role, (x, y) in (("forward", (a, b)), ("backward", (b, a))):
            for rid in motif.roles.get(role, []):
                r = reactions.get(rid)
                if r and not any(c.kind == "couple" and c.a.token == x and c.b.token == y for c in r.cofactors):
                    problems.append(f"{where}: {rid} does not turn {x} into {y}")
    elif motif.kind == "transport":
        for rid in motif.roles.get("via", []):
            r = reactions.get(rid)
            if r is None:
                continue
            species = {p.species for p in r.participants + r.consumed + r.produced}
            missing = [m for m in motif.members if m not in species]
            if missing or len({m.split("@")[1] for m in motif.members if "@" in m}) < 2:
                problems.append(f"{where}: {rid} does not move {', '.join(motif.members)} between compartments")
    return problems


def check_motifs(doc: Document) -> list[str]:
    return [p for m in doc.motifs for p in check_motif(doc, m)]


__all__ = ["check_motif", "check_motifs", "detect_motifs", "format_motif", "find_cycles", "find_shunts",
           "find_substrate_cycles", "find_repeats", "find_carriers", "find_transport"]
