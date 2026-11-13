"""Warnings about a spec that is formally right but breaks a rule of PROMPT.md (check-spec prints them).

These are the rules that are easy to overlook in the instructions: no cofactor merges, a copy for reactions with very
many metabolites, no long lines, no nodes on top of each other.
"""
import math
from .inspect_map import base, COFACTORS

#: a reaction with at least this many primary metabolites is "biomass-like": it gets copies, not the hub nodes
MANY = 6
#: a reaction whose primary nodes are further apart than this draws a long line (map units)
LONG = 3000
CLOSE = 70


def lint(spec, L):
    """Warnings (strings) for `spec`, applied to Layout `L` by declspec.compile_spec (partial or complete)."""
    W = []
    N, R = L.N, L.R
    prim = {rid: [n for _, n in I['subs'] + I['prods'] if n in N and N[n]['node_is_primary']] for rid, I in L.info.items()}
    alias = getattr(L, 'spec_alias', {})
    for m in spec.get('merges', []):
        keep = alias.get(m['keep'], m['keep'])
        if keep not in N: continue
        name = N[keep]['bigg_id']
        if base(name) in COFACTORS:
            W.append(f"merge {m['keep']}: {name} is a cofactor; cofactors are never merged (remove this merge)")
        users = [rid for rid, ns in prim.items() if keep in ns]
        big = [rid for rid in users if len(prim[rid]) >= MANY]
        if big and len(users) > len(big):
            W.append(f"merge {m['keep']} ({name}): it pulls {', '.join(R[r]['bigg_id'] for r in big)} ({len(prim[big[0]])} primary metabolites) onto the "
                     f"hub node; a reaction with that many metabolites keeps its own copy (leave its node out of \"dups\")")
    for rid, ns in prim.items():
        pts = [L.pos[n] for n in ns if n in L.pos]
        if len(pts) < 2 or len(ns) >= MANY: continue
        span = max(math.hypot(p[0] - q[0], p[1] - q[1]) for p in pts for q in pts)
        if span > LONG:
            W.append(f"reaction {rid} ({R[rid]['bigg_id']}): its metabolites are {span:.0f} units apart (a long line); move the blocks together, "
                     f"or keep a copy of the far metabolite (undo its merge, or a `splits` entry)")
    placed = [(n, L.pos[n]) for n in getattr(L, 'spec_placed', {}) if n in L.pos]
    close = [(a, b) for i, (a, p) in enumerate(placed) for b, q in placed[i + 1:] if math.hypot(p[0] - q[0], p[1] - q[1]) < CLOSE]
    if close:
        W.append(f'{len(close)} pairs of primary nodes are closer than {CLOSE} units: '
                 + ', '.join(f"{N[a]['bigg_id']}#{a} / {N[b]['bigg_id']}#{b}" for a, b in close[:8]))
    return W
