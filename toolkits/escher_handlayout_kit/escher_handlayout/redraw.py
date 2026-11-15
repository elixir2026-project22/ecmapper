"""Regular shapes for the blocks of a start spec: a cycle becomes a ring, a wobbling chain a straight column or row.

The draft draws every subsystem as its template does, which is often irregular. Before the blocks are connected
(`connect.connected_start`) this module redraws the two shapes a reader looks for:

- a block with a cycle of at least MIN_RING primary metabolites (the longest one) gets this cycle as a circle: the
  nodes of the cycle at equal angles, in their order and sense of rotation, around the old centre; everything else
  of the block keeps its place relative to the cycle (inside stays inside, outside keeps its distance). The
  reactions along the circle become arcs (`arc` with the centre), parallel ones in lanes;
- else a chain of at least MIN_CHAIN metabolites that runs down (or across) the block and wobbles is put on one
  straight line; what hangs on it moves with it.

Both are warps of the block's own coordinates: a function old point -> new point, applied to the nodes; the
entries of the reactions (`shift`, `axis`, `mid`, `side`) are derived again from the warped geometry.
"""
import math
from statistics import median

MIN_RING = 6              # primary metabolites of a cycle that is drawn as a ring
MIN_CHAIN = 4             # primary metabolites of a chain that is straightened
RING_SPACING = (300, 600) # allowed distance between neighbours on a ring (the radius follows from it)
LANE = 150                # radial distance between parallel reactions on a ring
STRAIGHT = 0.08           # a chain is straightened if its nodes leave the line by at most this share of its length
MANY = 6                  # reactions with this many primary metabolites do not define the shape
SIDES = {'L': (-1, 0), 'R': (1, 0), 'U': (0, -1), 'D': (0, 1)}


def _unit(dx, dy):
    n = math.hypot(dx, dy)
    return (dx / n, dy / n) if n > 1e-9 else (1.0, 0.0)


def _sides(L, rid, nodes):
    N = L.N; I = L.info[rid]
    return ([n for _, n in I['subs'] if N[n]['node_is_primary'] and n in nodes], [n for _, n in I['prods'] if N[n]['node_is_primary'] and n in nodes])


def _graph(b, L):
    """Neighbours of the primary nodes of a block, and the reactions between each pair."""
    adj = {n: set() for n in b['nodes']}; pairs = {}
    for r in b['reactions']:
        s, p = _sides(L, r['id'], adj)
        if len(set(s + p)) >= MANY: continue
        for x in s:
            for y in p:
                if x != y: adj[x].add(y); adj[y].add(x); pairs.setdefault(frozenset((x, y)), []).append(r['id'])
    return adj, pairs


def _longest_cycle(adj, P, budget=200000):
    """The longest cycle; of several, the one that encloses most in the draft (P: the coordinates)."""
    def area(c): return abs(sum(P[c[i]][0] * P[c[(i + 1) % len(c)]][1] - P[c[(i + 1) % len(c)]][0] * P[c[i]][1] for i in range(len(c)))) / 2
    order = sorted(adj, key=int); best = []; left = [budget]
    for i, start in enumerate(order):
        allowed = set(order[i:])                      # the start is the smallest node of its cycle
        def dfs(v, path, seen):
            left[0] -= 1
            if left[0] < 0: return
            for w in sorted(adj[v], key=int):
                if w == start:
                    if len(path) >= 3 and (len(path) > len(best) or (len(path) == len(best) and area(path) > area(best) + 1)): best[:] = path
                elif w in allowed and w not in seen:
                    seen.add(w); dfs(w, path + [w], seen); seen.discard(w)
        dfs(start, [start], {start})
    return best


def _geometry(r, s, p, P):
    """(middle, unit axis) of a reaction entry in the block's coordinates `P`."""
    mean = lambda ns: (sum(P[n][0] for n in ns) / len(ns), sum(P[n][1] for n in ns) / len(ns))
    if 'mid' in r and 'axis' in r: return tuple(r['mid']), _unit(*r['axis'])
    if not s or not p: return (mean(s or p) if (s or p) else (0.0, 0.0)), (1.0, 0.0)
    a, c = mean(s), mean(p); sh = r.get('shift', (0, 0))
    return ((a[0] + c[0]) / 2 + sh[0], (a[1] + c[1]) / 2 + sh[1]), _unit(c[0] - a[0], c[1] - a[1])


def _rederive(r, s, p, P, f):
    """The entry of a reaction after the warp f (P: the old coordinates of the block's nodes)."""
    m, u = _geometry(r, s, p, P)
    if 'mid' in r and 'axis' in r:
        m2 = f(m); a, c = f((m[0] - 60 * u[0], m[1] - 60 * u[1])), f((m[0] + 60 * u[0], m[1] + 60 * u[1]))
        u2 = _unit(c[0] - a[0], c[1] - a[1])
        r['mid'] = [round(m2[0]), round(m2[1])]; r['axis'] = [round(u2[0], 3), round(u2[1], 3)]
    else:
        Q = {n: f(P[n]) for n in s + p}
        _, u2 = _geometry({}, s, p, Q)
        if 'shift' in r:
            off = r['shift'][0] * -u[1] + r['shift'][1] * u[0]
            r['shift'] = [round(-u2[1] * off), round(u2[0] * off)]
    if r.get('side') in SIDES:
        v = SIDES[r['side']]; a, c = f(m), f((m[0] + 80 * v[0], m[1] + 80 * v[1]))
        if abs(u2[0]) < abs(u2[1]): r['side'] = 'R' if c[0] - a[0] >= 0 else 'L'
        else: r['side'] = 'D' if c[1] - a[1] >= 0 else 'U'


def _apply(b, L, f, skip=()):
    P = {n: tuple(xy) for n, xy in b['nodes'].items()}
    for r in b['reactions']:
        if r['id'] in skip: continue
        s, p = _sides(L, r['id'], P)
        _rederive(r, s, p, P, f)
    for n, xy in P.items():
        q = f(xy); b['nodes'][n] = [round(q[0]), round(q[1])]


def _ring(b, L, cyc, pairs):
    """Draw the cycle `cyc` of block b as a circle. Returns the note, or None if the cycle does not go round its centre."""
    P = {n: tuple(xy) for n, xy in b['nodes'].items()}; m = len(cyc); two_pi = 2 * math.pi
    C = (sum(P[n][0] for n in cyc) / m, sum(P[n][1] for n in cyc) / m)
    ang = [math.atan2(P[n][1] - C[1], P[n][0] - C[0]) for n in cyc]
    rad = [math.hypot(P[n][0] - C[0], P[n][1] - C[1]) for n in cyc]
    d = [((ang[(i + 1) % m] - ang[i] + math.pi) % two_pi) - math.pi for i in range(m)]
    if all(x > 0.03 for x in d): sgn = 1
    elif all(x < -0.03 for x in d): sgn = -1
    else: return None
    if abs(abs(sum(d)) - two_pi) > 0.01 or min(rad) < 1: return None
    R = min(max(sum(rad) / m, m * RING_SPACING[0] / two_pi), m * RING_SPACING[1] / two_pi)
    step = sgn * two_pi / m
    t0 = math.atan2(sum(math.sin(ang[i] - i * step) for i in range(m)), sum(math.cos(ang[i] - i * step) for i in range(m)))

    def f(q):
        dx, dy = q[0] - C[0], q[1] - C[1]; r = math.hypot(dx, dy)
        if r < 1e-6: return C
        phi = math.atan2(dy, dx); i, t = m - 1, 1.0
        for j in range(m):
            a = ((phi - ang[j]) * sgn) % two_pi; w = d[j] * sgn
            if a <= w + 1e-12: i, t = j, a / w; break
        rb = rad[i] + (rad[(i + 1) % m] - rad[i]) * t
        th = t0 + (i + t) * step
        rn = R * r / rb if r <= rb else R + (r - rb)
        return (C[0] + rn * math.cos(th), C[1] + rn * math.sin(th))

    # the reactions along the circle: arcs, parallel ones in lanes (in the order the draft has them, from inside to outside)
    arcs = {}
    entry = {r['id']: r for r in b['reactions']}
    for i in range(m):
        x, y = cyc[i], cyc[(i + 1) % m]
        lanes = []
        for rid in pairs.get(frozenset((x, y)), []):
            s, p = _sides(L, rid, P)
            mid, _ = _geometry(entry[rid], s, p, P)
            lanes.append((math.hypot(mid[0] - C[0], mid[1] - C[1]), int(rid), rid, (x, y) if x in s else (y, x)))
        lanes.sort()
        for j, (_, _, rid, nodes) in enumerate(lanes):
            off = round((j - (len(lanes) - 1) / 2) * LANE)
            r = entry[rid]; side = r.get('side')
            if off: side = 'out' if off > 0 else 'in'
            elif side in SIDES:
                s, p = _sides(L, rid, P); mid, _ = _geometry(r, s, p, P); v = SIDES[side]
                side = 'out' if v[0] * (mid[0] - C[0]) + v[1] * (mid[1] - C[1]) >= 0 else 'in'
            arcs[rid] = dict(arc={'center': [round(C[0]), round(C[1])], 'offset': off}, arc_nodes=list(nodes), **({'side': side} if side else {}))
    _apply(b, L, f, skip=arcs)
    for rid, a in arcs.items():
        r = entry[rid]
        for k in ('shift', 'axis', 'mid', 'side'): r.pop(k, None)
        r.update(a)
    return f"'{b['block']}': the cycle {' - '.join(L.N[n]['bigg_id'] for n in cyc)} is drawn as a ring (radius {R:.0f})", set(cyc)


def _column(b, L, adj):
    """Put the longest chain that runs down or across block b on a straight line. Returns the note or None."""
    P = {n: tuple(xy) for n, xy in b['nodes'].items()}; best = None
    for ax in (1, 0):
        order = sorted(P, key=lambda n: (P[n][ax], int(n))); length = {}; prev = {}
        for v in order:
            length[v] = 1; prev[v] = None
            for w in adj[v]:
                if w in length and P[v][ax] - P[w][ax] >= 40 and length[w] + 1 > length[v]: length[v] = length[w] + 1; prev[v] = w
        end = max(order, key=lambda n: (length[n], -order.index(n)))
        path = [end]
        while prev[path[-1]] is not None: path.append(prev[path[-1]])
        path.reverse()
        # a node between two that are neighbours anyway is a detour (dhap between fdp and g3p), not a part of the chain
        path = [n for i, n in enumerate(path) if not (0 < i < len(path) - 1 and path[i + 1] in adj[path[i - 1]])]
        while len(path) >= MIN_CHAIN:                # ends that leave the line are no part of it either
            extent = P[path[-1]][ax] - P[path[0]][ax]; x0 = median(P[n][1 - ax] for n in path)
            dev = [abs(P[n][1 - ax] - x0) for n in path]
            if max(dev) <= STRAIGHT * extent: break
            path = path[1:] if dev[0] >= dev[-1] else path[:-1]
            if max(dev[1:-1], default=0) > STRAIGHT * extent and max(dev[0], dev[-1]) <= STRAIGHT * extent: path = []
        if len(path) >= MIN_CHAIN and (best is None or len(path) > len(best[1])): best = (ax, path)
    if best is None: return None
    ax, path = best; ts = [P[n][ax] for n in path]; x0 = median(P[n][1 - ax] for n in path)
    devs = [P[n][1 - ax] - x0 for n in path]
    if max(abs(v) for v in devs) < 8: return None

    def f(q):
        t = q[ax]
        if t <= ts[0]: dev = devs[0]
        elif t >= ts[-1]: dev = devs[-1]
        else:
            i = max(j for j in range(len(ts) - 1) if ts[j] <= t)
            dev = devs[i] + (devs[i + 1] - devs[i]) * (t - ts[i]) / (ts[i + 1] - ts[i])
        out = [q[0], q[1]]; out[1 - ax] -= dev
        return tuple(out)

    _apply(b, L, f)
    return f"'{b['block']}': the chain {' - '.join(L.N[n]['bigg_id'] for n in path)} is drawn as a straight {'column' if ax == 1 else 'row'}"


def redraw(spec, L):
    """Redraw the blocks of the start spec `spec` (changed in place). Returns (notes, ring nodes)."""
    notes = []; ring_nodes = set()
    for b in spec['layouts']:
        if b.get('rings') or len(b['nodes']) < MIN_CHAIN: continue
        adj, pairs = _graph(b, L)
        cyc = _longest_cycle(adj, b['nodes'])
        done = _ring(b, L, cyc, pairs) if len(cyc) >= MIN_RING else None
        if done: notes.append(done[0]); ring_nodes |= done[1]; continue
        note = _column(b, L, adj)
        if note: notes.append(note)
    return notes, ring_nodes
