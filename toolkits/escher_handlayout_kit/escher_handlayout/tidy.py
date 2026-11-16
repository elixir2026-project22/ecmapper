"""Tidy the drawing of a spec without changing its design: `run.py tidy`.

Two mechanical steps that are tedious to do by hand, coordinate by coordinate (from the BRENDA drawings of ecmapper,
LAYOUT_LEARNINGS.md 6 and 8):

- **rows and columns** (:func:`align`): primary metabolites that one reaction joins and that are almost above or beside
  each other are put exactly into one column or row, along chains of such reactions;
- **cofactors** (:func:`place_cofactors`): `side` and `dist` of a reaction are chosen so that its cofactors and their
  branches lie where the fewest nodes and lines of the map are;
- **untangling** and **copies** (:func:`tidy`): local changes with the rest of the map pinned (LAYOUT_LEARNINGS.md 2
  and 4): a loose end is turned around its reaction, a reaction with a long line to a hub gets its own copy of it.

All write plain values into the spec (node coordinates, `side`, `dist`, `splits`), which stay editable; nothing moves a
block or changes a merge. :func:`tidy` keeps a change only if the metrics of the drawing get better.
"""
import copy, math
import numpy as np
from .engine import Layout, bez
from .declspec import compile_spec

TOLERANCE, SPREAD = 60.0, 130.0        # align: a link so far off the upright; a row at most so wide before
DISTS = (110, 150, 190)                # cofactors: the distances tried besides the reaction's own
CELL = 60.0


class _Grid:
    """Weighted points of the map, each with the reaction it belongs to (None: a node)."""
    def __init__(self): self.cells = {}; self.of = {}
    def add(self, x, y, w, tag=None):
        key = (int(x // CELL), int(y // CELL)); self.cells.setdefault(key, []).append((x, y, w, tag))
        if tag is not None: self.of.setdefault(tag, set()).add(key)
    def remove(self, tag, weight=None):
        """Drop the points of `tag` (only those of `weight` if given)."""
        for key in self.of.get(tag, ()):
            self.cells[key] = [p for p in self.cells[key] if p[3] != tag or (weight is not None and p[2] != weight)]
    def near(self, x, y, r, skip=None):
        cx, cy = int(x // CELL), int(y // CELL); k = int(r // CELL) + 1; total = 0.0
        for i in range(cx - k, cx + k + 1):
            for j in range(cy - k, cy + k + 1):
                for px, py, w, tag in self.cells.get((i, j), ()):
                    if tag != skip and (px - x) ** 2 + (py - y) ** 2 < r * r: total += w
        return total


def _compiled(input_path, spec):
    """(Layout with `spec` applied, {reaction id: (zone, the keywords of its drawing)})."""
    L = Layout(input_path); calls = {}; draw = L.reaction
    def reaction(rid, zone=None, **kw):
        calls[str(rid)] = (zone, kw); return draw(rid, zone, **kw)
    L.reaction = reaction
    compile_spec(spec, L, partial=True)
    L.reaction = draw
    return L, calls


def _blocks(spec, blocks):
    names = [b['block'] for b in spec.get('layouts', [])]
    unknown = [b for b in blocks or [] if b not in names]
    if unknown: raise ValueError(f"no such block: {', '.join(unknown)} (blocks: {', '.join(names)})")
    return [b for b in spec.get('layouts', []) if not blocks or b['block'] in blocks]


# ------------------------------------------------------------------------------------------------ rows and columns
def align(input_path, spec, blocks=None, tolerance=TOLERANCE, spread=SPREAD):
    """Put the nodes of each block that are almost in a column or a row into it. Changes `spec`; returns notes."""
    L, _ = _compiled(input_path, spec); notes = []
    ident = lambda k: L.spec_alias.get(k, k) if k.startswith('@') else k
    for b in _blocks(spec, blocks):
        nodes = b.get('nodes', {}); key = {ident(k): k for k in nodes}
        if len(nodes) < 2: continue
        links = []
        for r in b.get('reactions', []):
            I = L.info.get(r['id'])
            if I is None or 'arc' in r or r['id'] not in L.spec_done: continue
            s = [n for _, n in I['subs'] if L.N[n]['node_is_primary'] and n in key]
            p = [n for _, n in I['prods'] if L.N[n]['node_is_primary'] and n in key]
            links += [(key[x], key[y]) for x in s for y in p if x != y]
        before = dict(nodes); moved = set()
        for here, across in ((0, 1), (1, 0)):            # columns (the same x), then rows
            parent = {}
            def find(i):
                while parent.setdefault(i, i) != i: parent[i] = parent[parent[i]]; i = parent[i]
                return i
            for x, y in links:
                off = abs(nodes[x][here] - nodes[y][here]); length = abs(nodes[x][across] - nodes[y][across])
                if off < tolerance and length > 4 * off: parent[find(x)] = find(y)
            groups = {}
            for i in list(parent): groups.setdefault(find(i), []).append(i)
            for members in groups.values():
                values = [nodes[i][here] for i in members]
                if len(members) < 2 or max(values) - min(values) > spread or max(values) - min(values) < 0.5: continue
                value = round(float(np.median(values)))  # the place of most of them: a chain keeps its line
                for i in members:
                    if abs(nodes[i][here] - value) >= 0.5:
                        xy = list(nodes[i]); xy[here] = value; nodes[i] = xy; moved.add(i)
        close = lambda P: {frozenset((i, j)) for i in P for j in P if i < j and math.dist(P[i], P[j]) < 70}
        if close(nodes) - close(before):                 # two nodes pushed onto each other: leave the block as it was
            nodes.clear(); nodes.update(before); continue
        if moved: notes.append(f"block '{b['block']}': {len(moved)} nodes put into rows and columns")
    return notes


# ------------------------------------------------------------------------------------------------------ cofactors
def _samples(a, b, step=25.0):
    n = max(1, int(math.dist(a, b) // step))
    return [(a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n) for i in range(1, n)]


def _cofactors(L, rid):
    I = L.info[rid]
    return [(n, I['sm']) for _, n in I['subs'] if not L.N[n]['node_is_primary']] + \
           [(n, I['pm']) for _, n in I['prods'] if not L.N[n]['node_is_primary']]


def _sides(L, rid, kw):
    """The sides the cofactors of a reaction can hang on (the two normals of its line)."""
    if 'arc' in kw or kw.get('side') in ('out', 'in'): return ['out', 'in']
    u = L.rxn_geom[rid]['u']
    return ['L', 'R'] if abs(u[1]) >= abs(u[0]) else ['U', 'D']


def _side_now(L, rid, kw):
    side = kw.get('side', 'R')
    if side in ('out', 'in'): return side
    if 'arc' in kw: return 'in'                           # on an arc every other side is the inside
    n = L.rxn_geom[rid]['n']                              # the engine turns the side into the normal that points there
    return ('R' if n[0] > 0 else 'L') if abs(n[0]) >= abs(n[1]) else ('D' if n[1] > 0 else 'U')


def place_cofactors(input_path, spec, blocks=None, rounds=2):
    """Choose `side` and `dist` of the reactions so that their cofactors lie in free space. Changes `spec`; returns notes."""
    L, calls = _compiled(input_path, spec)
    if len(calls) < len(L.R): return ['cofactors: not placed, the spec is not complete yet']
    chosen = {r['id'] for b in _blocks(spec, blocks) for r in b.get('reactions', [])}
    grid = _Grid()
    for k, v in L.N.items():
        if k not in L.pos: continue
        x, y = L.pos[k]
        if v['node_type'] == 'metabolite' and v['node_is_primary']:
            grid.add(x, y, 10.0)
            for d in (-45.0, 45.0): grid.add(x + d, y - 20, 2.0)        # room for its label
    for rid in calls:
        I = L.info[rid]; mid = L.pos[I['mid']]
        grid.add(*mid, 5.0, rid)
        for nodes, marker in ((I['subs'], I['sm']), (I['prods'], I['pm'])):
            for p in _samples(L.pos[marker], mid) + [L.pos[marker]]: grid.add(*p, 1.5, rid)
            for _, n in nodes:
                if L.N[n]['node_is_primary']:
                    for p in _samples(L.pos[n], L.pos[marker]): grid.add(*p, 1.5, rid)
    # the lines of the map as they are drawn: the main lines fixed, the branches of the cofactors as straight lines
    L.finish(); mainA, mainB, mainR = [], [], []; index = {rid: i for i, rid in enumerate(calls)}
    for rid, r in L.R.items():
        cof = {n for n, _ in _cofactors(L, rid)}
        for sg in r['segments'].values():
            a, b = sg['from_node_id'], sg['to_node_id']
            if a in cof or b in cof: continue
            pa, pb = np.array(L.pos[a]), np.array(L.pos[b])
            pts = bez(pa, [sg['b1']['x'], sg['b1']['y']], [sg['b2']['x'], sg['b2']['y']], pb, max(6, int(np.linalg.norm(pb - pa) / 30))) if sg['b1'] else np.array([pa, pb])
            if len(pts) > 2: pts = pts[(np.linalg.norm(pts - pa, axis=1) > 40) & (np.linalg.norm(pts - pb, axis=1) > 40)]
            for i in range(len(pts) - 1): mainA.append(pts[i]); mainB.append(pts[i + 1]); mainR.append(index[rid])
    mainA, mainB, mainR = np.array(mainA).reshape(-1, 2), np.array(mainB).reshape(-1, 2), np.array(mainR, int)
    def branches(rid):
        out = []
        for n, marker in _cofactors(L, rid):
            a, b = np.array(L.pos[marker], float), np.array(L.pos[n], float); d = np.linalg.norm(b - a)
            if d > 85: out.append((a + (b - a) * 40 / d, b - (b - a) * 40 / d))        # qa ignores the 40 units next to nodes
        return out
    drawn = {rid: branches(rid) for rid in calls}
    def crossings(rid, A, B, Rr):
        """How many lines of other reactions the branches of `rid` cross."""
        total = 0; keep = Rr != index[rid]; A, B = A[keep], B[keep]
        cr = lambda a, b, c: (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
        for p, q in branches(rid):
            total += int((((cr(p, q, A) * cr(p, q, B)) < 0) & ((cr(A, B, p) * cr(A, B, q)) < 0)).sum())
        return total
    def others(rid):
        extra = [(a, b, index[r]) for r, bs in drawn.items() if r != rid for a, b in bs]
        if not extra: return mainA, mainB, mainR
        return (np.vstack([mainA, [e[0] for e in extra]]), np.vstack([mainB, [e[1] for e in extra]]), np.concatenate([mainR, [e[2] for e in extra]]))
    COF = 6.0
    def put(rid):
        for n, marker in _cofactors(L, rid):
            grid.add(*L.pos[n], COF, rid)
            for p in _samples(L.pos[marker], L.pos[n], 30.0): grid.add(*p, 1.0001, rid)
    def take(rid): grid.remove(rid, COF); grid.remove(rid, 1.0001)
    def draw(rid, **change):
        zone, kw = calls[rid]; kw = dict(kw, **change)
        for n, _ in _cofactors(L, rid): L.pos.pop(n, None)
        L.reaction(rid, zone, **kw)
    def cost(rid):
        n = L.rxn_geom[rid]['n']; total = 0.0
        for node, marker in _cofactors(L, rid):
            x, y = L.pos[node]
            total += grid.near(x, y, 48.0, rid) + 0.5 * grid.near(x + 30 * n[0], y + 30 * n[1], 30.0, rid)
            total += sum(grid.near(*p, 22.0, rid) for p in _samples(L.pos[marker], (x, y), 30.0))
        return total
    for rid in calls: put(rid)
    todo = sorted((r for r in calls if r in chosen and _cofactors(L, r)), key=lambda r: (-len(_cofactors(L, r)), int(r)))
    first = {r: (_side_now(L, r, calls[r][1]), calls[r][1].get('dist', 150)) for r in todo}
    for _ in range(rounds):
        for rid in todo:
            zone, kw = calls[rid]; side0, dist0 = _side_now(L, rid, kw), kw.get('dist', 150)
            take(rid); best = None; lines = others(rid)
            for side in _sides(L, rid, kw):
                for dist in sorted({dist0, *DISTS}):
                    draw(rid, side=side, dist=dist)
                    c = 4.0 * crossings(rid, *lines) + 0.2 * cost(rid) + (0.0 if side == side0 else 0.4) + abs(dist - dist0) / 200.0 + dist / 600.0
                    if (side, dist) == (side0, dist0): c -= 0.3                      # a change must pay
                    if best is None or c < best[0] - 1e-9: best = (c, side, dist)
            kw['side'], kw['dist'] = best[1], best[2]
            draw(rid); put(rid); drawn[rid] = branches(rid)
    entry = {r['id']: r for b in spec.get('layouts', []) for r in b.get('reactions', [])}
    sides = dists = 0; names = []
    for rid in todo:
        kw = calls[rid][1]; now = (kw['side'], kw['dist'])
        if now == first[rid]: continue
        sides += now[0] != first[rid][0]; dists += now[0] == first[rid][0]
        entry[rid]['side'], entry[rid]['dist'] = now; names.append(L.R[rid]['bigg_id'])
    if not names: return []
    return [f"cofactors: {sides} reactions got another side, {dists} another distance (" + ', '.join(names[:30]) + (' ...' if len(names) > 30 else '') + ')']


# ------------------------------------------------------------------------------- only what the metrics confirm
WEIGHT = {'crossings': 1, 'labels on lines': 3, 'labels on nodes': 6, 'label overlaps': 4, 'close nodes': 10, 'long lines': 3,
          'unconnected subsystems': 40, 'detached pieces': 20, 'placement penalty': 1}
TURNS = (30, -30, 60, -60, 90, -90, 180)   # untangle: a loose end is tried at these angles around the rest of its reaction
COPY_COST = 3                              # copies: a copy of a metabolite must save more than so many points
COPY_REACH, COPY_AT = 500.0, 400.0         # ... is tried for a line longer than this, and put so far from the reaction's other end


def measure(input_path, spec, curv=0.6, taxa='all'):
    """(value, parts, the reactions named) of the drawing of `spec`: the terms of the score that tidying can change,
    measured as `build` does."""
    from .engine import place_labels
    from . import pack, qa, directions, pipeline
    L = Layout(input_path); compile_spec(spec, L, partial=False); L.curv = curv
    L.finish(); place_labels(L); pack.pack(L)
    m = qa.metrics(L); pl = directions.check_layout(L, taxa) or {}
    parts = {'crossings': m['line_crossings'], 'labels on lines': m['labels_on_lines'], 'labels on nodes': m['labels_on_nodes'],
             'label overlaps': m['label_label_overlaps'], 'close nodes': m['primary_nodes_closer_than_70'], 'long lines': m.get('long_lines', 0),
             'unconnected subsystems': len(pipeline.unconnected_subsystems(L, spec.get('merges', []))),
             'detached pieces': len(pipeline.detached_pieces(L)), 'placement penalty': pl.get('penalty', 0)}
    named = {x for pair, _ in m['worst_crossing_pairs'] for x in pair} | {x.split(' ', 1)[1] for x in m['labels_on_lines_list'] if x.startswith('reaction ')}
    return sum(WEIGHT[k] * v for k, v in parts.items()), parts, named


def _primary(L, rid):
    I = L.info[rid]; prim = lambda xs: [n for _, n in xs if L.N[n]['node_is_primary']]
    return prim(I['subs']), prim(I['prods'])


def tidy(input_path, spec, blocks=None, do_align=True, do_untangle=True, do_cofactors=True, do_copies=False, seconds=90.0, taxa='all'):
    """A tidied copy of `spec` and the notes of what changed. Every change is a proposal that is kept only if the
    metrics of the drawing get better (:func:`measure`):

    - rows and columns (:func:`align`), block by block;
    - untangling: a metabolite that only one reaction uses (a loose end) is turned around the rest of its reaction, for
      the reactions that the metrics name, and a lane of parallel reactions is put on the other side of its line;
      everything else stays pinned;
    - cofactors (:func:`place_cofactors`), block by block, then the other side and other distances for the reactions
      that the metrics still name;
    - copies (only with `do_copies`): a named reaction with a long line to a metabolite that three or more reactions
      use gets its own copy of it beside its other metabolites (a `splits` entry), if that saves more than COPY_COST points.

    `seconds`: when to stop trying (shared between the steps)."""
    import time
    spec = copy.deepcopy(spec); names = [b['block'] for b in _blocks(spec, blocks)]
    value, first, named = measure(input_path, spec, taxa=taxa); notes = []
    block = lambda s, name: next(b for b in s['layouts'] if b['block'] == name)
    steps = [do_align, do_untangle, do_cofactors, do_copies]; start = time.time(); end = start
    def next_step(i):                        # each step gets its share of the time that is left
        nonlocal end
        left = start + seconds - time.time(); end = time.time() + max(0.0, left) / max(1, sum(steps[i:]))
    def accept(trial, margin=0):
        nonlocal spec, value, named
        try: v, _, who = measure(input_path, trial, taxa=taxa)
        except Exception: return False       # a proposal that does not compile is no improvement
        if v + margin < value: spec, value, named = trial, v, who; return True
        return False
    if do_align:
        next_step(0); proposal = copy.deepcopy(spec); align(input_path, proposal, blocks); kept = []
        for name in names:
            if block(proposal, name).get('nodes') == block(spec, name).get('nodes') or time.time() > end: continue
            trial = copy.deepcopy(spec); block(trial, name)['nodes'] = block(proposal, name)['nodes']
            if accept(trial): kept.append(name)
        if kept: notes.append('rows and columns: nodes aligned in ' + ', '.join(kept))
    if do_untangle:
        next_step(1); moved = set(); tried = set(); improved = True
        while improved and time.time() < end:
            improved = False; L, _ = _compiled(input_path, spec)
            uses = {}
            for rid in L.info:
                for n in set(sum(_primary(L, rid), [])): uses.setdefault(n, []).append(rid)
            ident = lambda k: L.spec_alias.get(k, k) if k.startswith('@') else k
            for name in names:
                b = block(spec, name); entry = {r['id']: r for r in b.get('reactions', [])}
                for k in list(b.get('nodes', {})):
                    n = ident(k); rs = uses.get(n, [])
                    if len(rs) != 1 or rs[0] not in entry or L.R[rs[0]]['bigg_id'] not in named or {'arc', 'mid', 'axis'} & set(entry[rs[0]]): continue
                    s_, p_ = _primary(L, rs[0]); rest = p_ if n in s_ else s_
                    if not rest or n in rest: continue
                    a = np.mean([L.pos[x] for x in rest], 0); v = np.array(L.pos[n]) - a
                    for turn in TURNS:
                        if (k, turn, tuple(b['nodes'][k])) in tried or time.time() > end: continue
                        tried.add((k, turn, tuple(b['nodes'][k])))
                        c, s2 = math.cos(math.radians(turn)), math.sin(math.radians(turn))
                        d = np.array([c * v[0] - s2 * v[1], s2 * v[0] + c * v[1]]) - v
                        trial = copy.deepcopy(spec); xy = block(trial, name)['nodes'][k]
                        block(trial, name)['nodes'][k] = [round(xy[0] + d[0]), round(xy[1] + d[1])]
                        if accept(trial): moved.add(L.N[n]['bigg_id']); improved = True; break
                    if improved: break
                if improved: break                       # the positions changed: look at the map again
                # a lane of parallel reactions: on the other side of the line, if no other lane is there
                size = lambda sh: sh['offset'] if isinstance(sh, dict) else None
                for i, r in enumerate(b.get('reactions', [])):
                    sh = r.get('shift')
                    if sh is None or r['id'] not in L.info or L.R[r['id']]['bigg_id'] not in named or 'arc' in r: continue
                    mirror = dict(sh, offset=-sh['offset']) if isinstance(sh, dict) else [-sh[0], -sh[1]]
                    ends = frozenset(map(frozenset, _primary(L, r['id'])))
                    lanes = [q.get('shift', [0, 0]) for q in b['reactions'] if q is not r and q['id'] in L.info and frozenset(map(frozenset, _primary(L, q['id']))) == ends]
                    near = lambda q: (isinstance(q, dict) and isinstance(mirror, dict) and abs(q['offset'] - mirror['offset']) < 100) or                                      (isinstance(q, list) and isinstance(mirror, list) and math.dist(q, mirror) < 100)
                    key = (r['id'], 'lane', str(mirror))
                    if key in tried or any(near(q) for q in lanes) or time.time() > end: continue
                    tried.add(key); tried.add((r['id'], 'lane', str(sh)))
                    trial = copy.deepcopy(spec); block(trial, name)['reactions'][i]['shift'] = mirror
                    if accept(trial): moved.add(L.R[r['id']]['bigg_id'] + ' (lane)'); improved = True; break
                if improved: break
        if moved: notes.append(f'untangled: {len(moved)} loose ends turned around their reaction or lanes put on the other side of their line (' + ', '.join(sorted(moved)[:40]) + ')')
    if do_cofactors:
        next_step(2); changed = set()
        proposal = copy.deepcopy(spec); place_cofactors(input_path, proposal, blocks)
        for name in names:
            if block(proposal, name).get('reactions') == block(spec, name).get('reactions') or time.time() > end: continue
            trial = copy.deepcopy(spec); block(trial, name)['reactions'] = block(proposal, name)['reactions']
            old = {r['id']: r for r in block(spec, name)['reactions']}
            if accept(trial): changed |= {r.get('name', r['id']) for r in block(spec, name)['reactions'] if r != old.get(r['id'])}
        # the reactions that the metrics still name: the other side, other distances
        L, calls = _compiled(input_path, spec); tried = set(); improved = True
        while improved and time.time() < end:
            improved = False
            for name in names:
                for i, r in enumerate(block(spec, name).get('reactions', [])):
                    rid = r['id']; bigg = L.R[rid]['bigg_id'] if rid in L.R else None
                    if bigg not in named or rid not in calls or not _cofactors(L, rid): continue
                    kw = calls[rid][1]; side0, dist0 = _side_now(L, rid, kw), r.get('dist', 150)
                    other = next(s for s in _sides(L, rid, kw) if s != side0)
                    for side, dist in [(other, dist0)] + [(side0, d) for d in DISTS if d != dist0] + [(other, d) for d in DISTS if d != dist0]:
                        if (rid, side, dist) in tried or time.time() > end: continue
                        tried.add((rid, side, dist))
                        trial = copy.deepcopy(spec); block(trial, name)['reactions'][i].update(side=side, dist=dist)
                        if accept(trial):
                            changed.add(bigg); improved = True; kw['side'], kw['dist'] = side, dist
                            r = block(spec, name)['reactions'][i]; side0, dist0 = side, dist; break
        if changed: notes.append(f'cofactors: side or distance changed for {len(changed)} reactions (' + ', '.join(sorted(changed)[:40]) + ')')
    if do_copies:
        next_step(3); copies = []; tried = set(); improved = True; raw = Layout(input_path)
        while improved and time.time() < end:
            improved = False; L, _ = _compiled(input_path, spec)
            uses = {}
            for rid in L.info:
                for n in set(sum(_primary(L, rid), [])): uses.setdefault(n, []).append(rid)
            now = {d: m['keep'] for m in spec.get('merges', []) for d in m.get('dups', [])}
            split = {(x['reaction'], x['node']) for x in spec.get('splits', [])}
            for name in names:
                b = block(spec, name)
                for r in b.get('reactions', []):
                    rid = r['id']
                    if rid not in L.spec_done or L.R[rid]['bigg_id'] not in named or {'arc', 'mid', 'axis'} & set(r): continue
                    s_, p_ = _primary(L, rid); mid = np.array(L.pos[L.info[rid]['mid']])
                    for n in s_ + p_:
                        rest = p_ if n in s_ else s_
                        if len(uses.get(n, [])) < 3 or not rest or n in rest or np.linalg.norm(np.array(L.pos[n]) - mid) < COPY_REACH: continue
                        own = next((x for _, x in raw.info[rid]['subs'] + raw.info[rid]['prods'] if now.get(x, x) == n), None)
                        if own is None or (rid, own) in split or (rid, n) in tried or time.time() > end: continue
                        tried.add((rid, n))
                        a = np.mean([L.pos[x] for x in rest], 0); v = np.array(L.pos[n]) - a; at = a + v / np.linalg.norm(v) * COPY_AT
                        alias = ('@c_' + ''.join(ch if ch.isalnum() else '_' for ch in L.N[n]['bigg_id']))[:24] + '_' + rid
                        trial = copy.deepcopy(spec); o = b.get('origin', [0, 0])
                        trial.setdefault('splits', []).append({'alias': alias, 'node': own, 'reaction': rid})
                        tb = block(trial, name); tb['nodes'][alias] = [round(at[0] - o[0]), round(at[1] - o[1])]
                        if 'node_names' in tb: tb['node_names'][alias] = L.N[n]['bigg_id']
                        if accept(trial, COPY_COST): copies.append(f"{L.N[n]['bigg_id']} for {L.R[rid]['bigg_id']}"); improved = True; break
                    if improved: break
                if improved: break
        if copies: notes.append(f'copies: {len(copies)} reactions got their own copy of a metabolite, beside their other metabolites (' + ', '.join(copies) + '); note them in DECISIONS.md, or undo one by deleting its `splits` entry and the @alias coordinates')
    _, last, _ = measure(input_path, spec, taxa=taxa) if notes else (0, first, 0)
    diff = [f'{k} {first[k]} -> {last[k]}' for k in first if first[k] != last[k]]
    notes.append('measured: ' + ', '.join(diff) if diff else 'nothing to tidy: no change made the drawing better by the metrics')
    return spec, notes
