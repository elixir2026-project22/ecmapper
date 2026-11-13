"""A connected starting spec: the subsystems of the draft, arranged and joined through shared metabolites.

`start.start_spec` gives the draft as it is: every subsystem a separate block. Arranging the blocks and choosing the
merges is the part of the layout that needs spatial planning, and the part a weaker language model fails at. So this
module does a first version of it, deterministically:

1. the first block (glycolysis if the map has it, else the block that shares most) stays where it is;
2. the other blocks are attached one after the other, the one that shares most metabolites with the placed blocks
   first. A block is moved as a whole (its drawing is kept). For every metabolite it shares, its own node of the
   metabolite is merged into the node of the placed blocks, and the block is put where that costs least: close to
   the kept node (short lines), not on top of what is placed, where the placement catalogue wants the subsystem,
   and with few crossings of the new lines;
3. further copies of a metabolite are merged where the two nodes ended up near each other and the lines stay free;
   the other copies remain copies (the build shows them, a later patch can merge them).

The geometry used here is a sketch (straight lines from each metabolite to the middle of its reaction); the build
measures the real map. Blocks that share nothing but cofactors are left alone: the packer places them.
"""
import math
import numpy as np
from .declspec import SpecError
from .engine import Layout
from .inspect_map import base, COFACTORS
from .start import start_spec
from . import directions as D

CELL = 50                 # raster of the occupancy grid (map units)
MARGIN = 5                # free cells kept around everything that is placed
STEP = 100                # raster of the positions tried for a block
REACH = 4500              # a block is tried up to this far from the node it is merged into
REACH_FAR = 9000          # ... and this far if there is no free place nearer
COFACTOR_REACH = 220      # room of the cofactors around the middle of a reaction
W_DIRECTION = 500         # cost of one penalty point of a violated direction, in map units of line length
W_CROSSING = 400          # cost of one crossing of the sketch
W_COMPACT = 0.15          # cost of the distance to the centre of what is placed
NEAR = 1500               # further copies are merged up to this distance
MANY = 6                  # a reaction with this many primary metabolites keeps its own copies (lint.MANY)


def _crossings(P, Q):
    """Number of proper crossings between the segments P (n, 4) and Q (m, 4): rows x1, y1, x2, y2."""
    if len(P) == 0 or len(Q) == 0: return 0
    P = np.asarray(P, float)[:, None, :]; Q = np.asarray(Q, float)[None, :, :]
    def orient(ax, ay, bx, by, cx, cy): return np.sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
    o1 = orient(P[..., 0], P[..., 1], P[..., 2], P[..., 3], Q[..., 0], Q[..., 1])
    o2 = orient(P[..., 0], P[..., 1], P[..., 2], P[..., 3], Q[..., 2], Q[..., 3])
    o3 = orient(Q[..., 0], Q[..., 1], Q[..., 2], Q[..., 3], P[..., 0], P[..., 1])
    o4 = orient(Q[..., 0], Q[..., 1], Q[..., 2], Q[..., 3], P[..., 2], P[..., 3])
    return int(np.sum((o1 * o2 < 0) & (o3 * o4 < 0)))


class _Sketch:
    """The blocks of a start spec as rigid bodies: global node positions, reaction middles, straight-line segments."""

    def __init__(self, spec, L):
        self.spec, self.L = spec, L
        N = L.N
        self.blocks = {b['block']: b for b in spec['layouts']}
        self.block_of = {}                              # reaction / node -> block
        self.pos = {}                                   # primary node -> global position
        self.prim = {}                                  # reaction -> its primary nodes
        self.entry = {}                                 # reaction -> its entry in the spec
        for b in spec['layouts']:
            o = b['origin']
            for n, (x, y) in b['nodes'].items(): self.pos[n] = (o[0] + x, o[1] + y); self.block_of['n' + n] = b['block']
            for r in b['reactions']:
                I = L.info[r['id']]
                self.prim[r['id']] = list(dict.fromkeys(n for _, n in I['subs'] + I['prods'] if N[n]['node_is_primary']))
                self.block_of[r['id']] = b['block']; self.entry[r['id']] = r
        self.users = {}                                 # primary node -> reactions
        for rid, ns in self.prim.items():
            for n in ns: self.users.setdefault(n, []).append(rid)
        self.keep = {}                                  # merged node -> the node that stands for it

    def rep(self, n):
        while n in self.keep: n = self.keep[n]
        return n

    def name(self, n): return self.L.N[n]['bigg_id']

    def carbon(self, n):
        """A node that may be merged: no cofactor, and not only used by reactions with very many metabolites."""
        return base(self.name(n)) not in COFACTORS and any(len(self.prim[r]) < MANY for r in self.users.get(n, []))

    def middle(self, rid):
        """Where the sketch puts the middle of a reaction."""
        r = self.entry[rid]; o = self.blocks[self.block_of[rid]]['origin']
        if 'mid' in r: return (o[0] + r['mid'][0], o[1] + r['mid'][1])
        pts = [self.pos[self.rep(n)] for n in self.prim[rid]]
        s = r.get('shift', (0, 0))
        return (sum(p[0] for p in pts) / len(pts) + s[0], sum(p[1] for p in pts) / len(pts) + s[1])

    def segments(self, rids, skip=None):
        """The straight lines node - middle of the reactions `rids` (without those that end at node `skip`)."""
        out = []
        for rid in rids:
            c = self.middle(rid)
            for n in self.prim[rid]:
                n = self.rep(n)
                if n != skip: out.append((self.pos[n][0], self.pos[n][1], c[0], c[1]))
        return out

    def points(self, block, skip=None):
        """What a block occupies: its nodes, samples of its lines, the room of its cofactors (global coordinates)."""
        b = self.blocks[block]; pts = [self.pos[n] for n in b['nodes'] if n != skip]
        for r in b['reactions']:
            c = self.middle(r['id'])
            pts += [(c[0] + dx, c[1] + dy) for dx in (-COFACTOR_REACH, 0, COFACTOR_REACH) for dy in (-COFACTOR_REACH, 0, COFACTOR_REACH)]
        for x1, y1, x2, y2 in self.segments([r['id'] for r in b['reactions']], skip):
            k = max(2, int(math.hypot(x2 - x1, y2 - y1) / CELL))
            pts += [(x1 + (x2 - x1) * i / k, y1 + (y2 - y1) * i / k) for i in range(k + 1)]
        return np.array(pts, float).reshape(-1, 2)

    def move(self, block, dx, dy):
        b = self.blocks[block]; b['origin'] = [b['origin'][0] + dx, b['origin'][1] + dy]
        for n in b['nodes']: self.pos[n] = (self.pos[n][0] + dx, self.pos[n][1] + dy)

    def merge(self, k, d):
        """Node d is drawn as node k from now on. A reaction that hung on d alone moves with it."""
        b = self.blocks[self.block_of['n' + d]]; o = b['origin']
        shift = (self.pos[k][0] - self.pos[d][0], self.pos[k][1] - self.pos[d][1])
        for rid in self.users.get(d, []):
            r = self.entry[rid]
            if 'mid' in r and all(self.rep(n) == d for n in self.prim[rid]) and self.block_of[rid] == b['block']:
                r['mid'] = [round(r['mid'][0] + shift[0]), round(r['mid'][1] + shift[1])]
        self.keep[d] = k
        b['nodes'].pop(d, None); b.get('node_names', {}).pop(d, None); self.pos.pop(d, None)
        self.users.setdefault(k, []).extend(self.users.pop(d, []))


class _Grid:
    """Occupancy raster with a margin around every point that is added."""

    def __init__(self, extent=40000):
        self.n = int(2 * extent / CELL); self.off = extent; self.g = np.zeros((self.n, self.n), bool)

    def cells(self, pts):
        return np.clip(((np.asarray(pts, float) + self.off) // CELL).astype(int), 0, self.n - 1)

    def add(self, pts):
        for j, i in {tuple(c) for c in self.cells(pts)}:
            self.g[max(i - MARGIN, 0):i + MARGIN + 1, max(j - MARGIN, 0):j + MARGIN + 1] = True

    def free(self, pts, shifts):
        """For every shift (m, 2): whether all points (p, 2), moved by it, are on free cells."""
        c = self.cells(pts[None, :, :] + shifts[:, None, :])
        return ~self.g[c[..., 1], c[..., 0]].any(axis=1)


def connected_start(input_path, taxa='all'):
    """(spec, report): the start spec of `input_path` with the subsystems arranged and connected; `report` is a list of
    lines saying what was merged, where the blocks went and which copies remain."""
    spec = start_spec(input_path); L = Layout(input_path)
    S = _Sketch(spec, L); report = []
    blocks = list(S.blocks)
    if len(blocks) < 2: return spec, ['only one subsystem: nothing to connect']

    def shared(block, others):
        """[(metabolite, node of `block`, node of `others`)] for the metabolites both draw."""
        mine = {}
        for n in S.blocks[block]['nodes']:
            if S.carbon(n): mine.setdefault(S.name(n), []).append(n)
        out = []
        for o in others:
            for n in S.blocks[o]['nodes']:
                if S.carbon(n) and S.name(n) in mine: out += [(S.name(n), d, n) for d in mine[S.name(n)]]
        return out

    # ---- the catalogue's directions between the subsystems
    arr = None; cat = D.Catalogue.load()
    if cat is not None:
        labels = {z['name']: z.get('label') or z.get('title') or z['name'] for z in spec['zones']}
        subs = {b: dict(title=labels.get(S.blocks[b]['zone'], b), metabolites=sorted({S.name(n) for n in S.blocks[b]['nodes']})) for b in blocks}
        arr = D.Arrangement(subs, cat, taxa)
    wanted = [d for d in (arr.subsystem_directions if arr else []) if d.strength in D.PENALTY]

    def centre(block, skip=None):
        p = [S.pos[n] for n in S.blocks[block]['nodes'] if n != skip]
        return (sum(q[0] for q in p) / len(p), sum(q[1] for q in p) / len(p)) if p else None

    def direction_cost(block, c, shifts, placed):
        """Penalty points of the directions between `block` (centre c, moved by each shift) and the placed blocks."""
        cost = np.zeros(len(shifts)); tol = 50.0
        for d in wanted:
            if d.anchor == block and d.target in placed: other, me_target = centre(d.target), False
            elif d.target == block and d.anchor in placed: other, me_target = centre(d.anchor), True
            else: continue
            if other is None: continue
            mx, my = c[0] + shifts[:, 0], c[1] + shifts[:, 1]
            # the target seen from the anchor
            dx = (mx - other[0]) if me_target else (other[0] - mx); north = (other[1] - my) if me_target else (my - other[1])
            parts = d.direction.split('_'); bad = np.zeros(len(shifts), bool)
            if 'above' in parts: bad |= north <= tol
            if 'below' in parts: bad |= north >= -tol
            if 'right' in parts: bad |= dx <= tol
            if 'left' in parts: bad |= dx >= -tol
            cost += bad * D.PENALTY[d.strength]
        return cost

    # ---- 1. the first block
    share_count = {b: len({m for m, _, _ in shared(b, [o for o in blocks if o != b])}) for b in blocks}
    glyc = [b for b in blocks if arr and 'glycolysis' in (arr.entities.get(b) or []) and share_count[b]]
    first = glyc[0] if glyc else max(blocks, key=lambda b: (share_count[b], len(S.blocks[b]['reactions'])))
    if not share_count[first]: return spec, ['the subsystems share no metabolite (cofactors apart): nothing to connect']
    placed = [first]; grid = _Grid(); grid.add(S.points(first))
    report.append(f"'{first}' stays where the draft draws it; the other subsystems were attached to it one after the other:")
    lattice = lambda reach: np.array([(x, y) for x in range(-reach, reach + 1, STEP) for y in range(-reach, reach + 1, STEP)], float)

    # ---- 2. attach the others
    while True:
        cands = {b: shared(b, placed) for b in blocks if b not in placed}
        cands = {b: v for b, v in cands.items() if v}
        if not cands: break
        block = max(cands, key=lambda b: (len({m for m, _, _ in cands[b]}), len(S.blocks[b]['reactions']), -blocks.index(b)))
        rids = [r['id'] for r in S.blocks[block]['reactions']]
        placed_rids = [r['id'] for p in placed for r in S.blocks[p]['reactions']]
        placed_segs = S.segments(placed_rids)
        pc = np.mean([centre(p) for p in placed], axis=0)
        best = None
        for m, d, k in sorted(set((m, d, S.rep(k)) for m, d, k in cands[block])):
            pts = S.points(block, skip=d); c = centre(block, skip=d)
            if c is None or len(pts) == 0: continue
            own = np.array(S.segments(rids, skip=d), float).reshape(-1, 4)
            mids = [S.middle(r) for r in S.users[d] if S.block_of[r] == block]           # the lines that will end at k
            for reach in (REACH, REACH_FAR):
                # shift of the block such that d comes to lie at k + offset
                off = lattice(reach); shifts = off + (np.array(S.pos[k]) - np.array(S.pos[d]))
                ok = grid.free(pts, shifts)
                if ok.any(): break
            if not ok.any(): continue
            off, shifts = off[ok], shifts[ok]
            cost = np.hypot(off[:, 0], off[:, 1]) + W_DIRECTION * direction_cost(block, c, shifts, placed) \
                + W_COMPACT * np.hypot(c[0] + shifts[:, 0] - pc[0], c[1] + shifts[:, 1] - pc[1])
            for i in np.argsort(cost, kind='stable')[:40]:                               # the sketch's crossings for the best ones
                s = shifts[i]
                new = [(S.pos[k][0], S.pos[k][1], x + s[0], y + s[1]) for x, y in mids]
                moved = own + np.array([s[0], s[1], s[0], s[1]]) if len(own) else own
                total = cost[i] + W_CROSSING * (_crossings(new, placed_segs) + _crossings(new, moved) + _crossings(moved, placed_segs))
                if best is None or total < best[0] - 1e-9: best = (total, m, d, k, (float(s[0]), float(s[1])))
        if best is None:                                                                # no room anywhere: it stays a separate block
            report.append(f"'{block}' found no free place next to the others: it stays a separate block")
            placed.append(block); grid.add(S.points(block)); continue
        _, m, d, k, s = best
        S.move(block, round(s[0]), round(s[1])); S.merge(k, d)
        placed.append(block); grid.add(S.points(block))
        report.append(f"'{block}' through {m} (node {d} merged into {k}, in '{S.block_of['n' + k]}')")

    # ---- 3. further copies that ended up near each other
    all_rids = list(S.entry)
    extra = []
    while True:
        groups = {}
        for b in placed:
            for n in S.blocks[b]['nodes']:
                if S.carbon(n): groups.setdefault(S.name(n), []).append(n)
        best = None
        for m, ns in sorted(groups.items()):
            for d in ns:
                for k in ns:
                    if d == k or (m, d, k) in extra or set(S.users[d]) & set(S.users[k]): continue
                    if len(S.users[k]) < len(S.users[d]) or (len(S.users[k]) == len(S.users[d]) and int(k) > int(d)): continue      # the busier node is kept
                    dist = math.hypot(S.pos[d][0] - S.pos[k][0], S.pos[d][1] - S.pos[k][1])
                    if dist <= NEAR and (best is None or dist < best[0]): best = (dist, m, d, k)
        if best is None: break
        dist, m, d, k = best; extra.append((m, d, k))
        mine = [r for r in S.users[d]]; others = [r for r in all_rids if r not in mine]
        before = _crossings(S.segments(mine), S.segments(others))
        new = []
        for r in mine:
            c = S.middle(r)
            new += [(S.pos[k][0], S.pos[k][1], c[0], c[1])] + [(S.pos[S.rep(n)][0], S.pos[S.rep(n)][1], c[0], c[1]) for n in S.prim[r] if S.rep(n) != d]
        if _crossings(new, S.segments(others)) > before: continue
        S.merge(k, d); report.append(f'also merged: {m} (node {d} into {k}, {dist:.0f} units apart)')

    # ---- the spec
    members = {}
    for d in S.keep: members.setdefault(S.rep(d), []).append(d)
    spec['merges'] = [dict(keep=k, dups=sorted(v, key=int)) for k, v in sorted(members.items(), key=lambda kv: int(kv[0]))]
    for b in placed:
        xs = [S.pos[n] for n in S.blocks[b]['nodes']]
        if xs: spec['titles'][S.blocks[b]['zone']] = [round(min(p[0] for p in xs)), round(min(p[1] for p in xs) - 130)]
    copies = {}
    for b in spec['layouts']:
        for n in b['nodes']:
            if S.carbon(n): copies.setdefault(S.name(n), []).append(n)
    left = {m: ns for m, ns in copies.items() if len(ns) > 1}
    if left: report.append('still drawn more than once (copies; merge them with a patch where the line stays short): '
                           + '; '.join(f"{m}: {', '.join(sorted(ns, key=int))}" for m, ns in sorted(left.items())))
    spec['notes'] = ('Starting spec written by `run.py start`: the subsystems of the draft, arranged by the placement catalogue and '
                     'connected through shared metabolites (a first arrangement: improve it in small steps).')
    return spec, report
