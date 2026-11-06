"""Layout engine: hand-placed coordinates -> Escher geometry (segments, curves, labels).

Nothing in here decides *where* things go; the spec file does that. The engine only converts the
chosen primary-node coordinates and per-reaction options into midmarkers, multimarkers, cofactor
positions, Bezier control points and label positions.
"""
import json, math, copy
import numpy as np
from ._io import read_json

def unit(v):
    v = np.array(v, float); n = np.linalg.norm(v)
    return v / n if n > 0 else v

DIRS = {'L': (-1, 0), 'R': (1, 0), 'U': (0, -1), 'D': (0, 1)}

def parse_reactions(m):
    N, R = m['nodes'], m['reactions']
    info = {}
    for rid, r in R.items():
        mid = sm = pm = None
        subs, prods = [], []
        for sid, s in r['segments'].items():
            a, b = s['from_node_id'], s['to_node_id']
            ta, tb = N[a]['node_type'], N[b]['node_type']
            if ta == 'midmarker': mid = a
            if tb == 'midmarker': mid = b
            if ta == 'metabolite' and tb == 'multimarker':
                subs.append((sid, a)); sm = b
            if ta == 'multimarker' and tb == 'metabolite':
                prods.append((sid, b)); pm = a
        info[rid] = dict(mid=mid, sm=sm, pm=pm, subs=subs, prods=prods)
    return info

PALETTE = ['#2b5c9e', '#5b8c2a', '#8e3b9e', '#d9822b', '#6a4fb3', '#8a7a2a', '#2f9a57', '#c9a100',
           '#c2415a', '#1a8f9a', '#a0522d', '#4a6fa5', '#7a7a7a', '#b5651d']

class Layout:
    def __init__(self, path):
        self.d = read_json(path)
        self.meta, self.m = self.d
        self.N, self.R = self.m['nodes'], self.m['reactions']
        self._dup_shared()
        self.info = parse_reactions(self.m)
        self.pos = {}          # node id -> (x,y)
        self.lab = {}          # node id -> (x,y,anchor)
        self.rlab = {}         # rid -> (x,y,anchor)
        self.rxn_zone = {}
        self.orig = {k: (v['x'], v['y']) for k, v in self.N.items()}
        self.zone_style = {}   # zone -> dict(color, title, label)
        self.titles = {}       # zone -> (x, y) of the zone title
        self.extra_titles = [] # filled by pack()
        self.small_titles = {} # bigg_id of a reaction in a small component -> (zone, title, reuse_label_text)
        self.bend = 0          # default shift of cofactor dots along the reaction axis (room for the curve)
        self.curv = 0.0        # curvature, 0 = gentle ... 1 = strong

    def zone(self, name, color=None, title=None, label=None):
        """Register a subsystem: line colour, title text, and (optionally) the text of the existing
        title label of the input map that should be moved/reused."""
        color = color or PALETTE[len(self.zone_style) % len(PALETTE)]
        self.zone_style[name] = dict(color=color, title=title or name, label=label)

    def _dup_shared(self):
        N, R = self.N, self.R
        use = {}
        for rid, r in R.items():
            nodes = set()
            for sg in r['segments'].values(): nodes |= {sg['from_node_id'], sg['to_node_id']}
            for n in sorted(nodes, key=int):      # sorted -> node ids of the copies are deterministic
                if N[n]['node_type'] == 'metabolite' and not N[n]['node_is_primary']:
                    use.setdefault(n, []).append(rid)
        nxt = max(int(k) for k in N) + 1
        self.dups = 0
        for n, rids in use.items():
            if len(rids) < 2: continue
            for rid in rids[1:]:
                new = str(nxt); nxt += 1; self.dups += 1
                N[new] = copy.deepcopy(N[n])
                for sg in R[rid]['segments'].values():
                    if sg['from_node_id'] == n: sg['from_node_id'] = new
                    if sg['to_node_id'] == n: sg['to_node_id'] = new

    def split(self, nid, rid):
        new = str(max(int(k) for k in self.N) + 1)
        self.N[new] = copy.deepcopy(self.N[nid]); self.orig[new] = self.orig[nid]
        for sg in self.R[str(rid)]['segments'].values():
            if sg['from_node_id'] == nid: sg['from_node_id'] = new
            if sg['to_node_id'] == nid: sg['to_node_id'] = new
        self.info = parse_reactions(self.m)
        return new

    def merge(self, keep, *others):
        keep = str(keep)
        for o in others:
            o = str(o)
            for r in self.R.values():
                for sg in r['segments'].values():
                    if sg['from_node_id'] == o: sg['from_node_id'] = keep
                    if sg['to_node_id'] == o: sg['to_node_id'] = keep
            del self.N[o]
        self.info = parse_reactions(self.m)

    def P(self, **kw):
        for k, v in kw.items():
            self.pos[k.lstrip('n')] = tuple(v)

    def setp(self, d):
        for k, v in d.items(): self.pos[str(k)] = tuple(v)

    def reaction(self, rid, zone=None, shift=(0, 0), side='R', dist=150, axis=None, mid=None,
                 lab=None, curve=None, step=75, gap=35, sec=None, center=None,
                 labdx=0, labdy=0, secdist=None, arc=None, arc_nodes=None, bend=None, bow=0):
        rid = str(rid); I = self.info[rid]; N = self.N
        if I['mid'] is None or I['sm'] is None or I['pm'] is None:
            raise ValueError(f"reaction {rid} ({self.R[rid]['bigg_id']}) has no multimarker topology (e.g. an exchange reaction "
                             f"drawn straight to its midmarker). This engine expects the Escher draft-map topology "
                             f"substrate -> multimarker -> midmarker -> multimarker -> product.")
        sec = {str(k): v for k, v in (sec or {}).items()}
        self.rxn_zone[rid] = zone
        if zone not in self.zone_style: self.zone(zone)
        subP = [n for _, n in I['subs'] if N[n]['node_is_primary']]
        prodP = [n for _, n in I['prods'] if N[n]['node_is_primary']]
        ps = [np.array(self.pos[n]) for n in subP]
        pp = [np.array(self.pos[n]) for n in prodP]
        bend = getattr(self, 'bend', 0) if bend is None else bend
        arcinfo = None
        if arc is not None:
            C = np.array(arc[0], float); off = arc[1]
            sa, sb = arc_nodes if arc_nodes else (subP[0], prodP[0])
            A = np.array(self.pos[sa]); B = np.array(self.pos[sb]); R0 = np.linalg.norm(A - C)
            ta = math.atan2(*(A - C)[::-1]); tb = math.atan2(*(B - C)[::-1])
            dth = (tb - ta + math.pi) % (2 * math.pi) - math.pi; sg = 1 if dth > 0 else -1
            tm = ta + dth / 2; Rm = R0 + off
            pol = lambda ang, rad: C + rad * np.array([math.cos(ang), math.sin(ang)])
            mid = pol(tm, Rm); dl = gap / Rm
            ms, mp = pol(tm - sg * dl, Rm), pol(tm + sg * dl, Rm)
            u = sg * np.array([-math.sin(tm), math.cos(tm)])
            rad = np.array([math.cos(tm), math.sin(tm)])
            n = rad if side == 'out' else -rad
            arcinfo = dict(C=C, s=sg, nodes=(sa, sb), Rm=Rm)
            self.pos[I['mid']] = tuple(mid); self.pos[I['sm']] = tuple(ms); self.pos[I['pm']] = tuple(mp)
            self.rxn_geom = getattr(self, 'rxn_geom', {})
            self.rxn_geom[rid] = dict(u=u, n=n, mid=mid, ms=ms, mp=mp, curve=False, arc=arcinfo)
            sd = secdist or dist
            subS = [x for x in I['subs'] if not N[x[1]]['node_is_primary']]
            prodS = [x for x in I['prods'] if not N[x[1]]['node_is_primary']]
            rsec = Rm + (sd if side == 'out' else -sd)
            for i, (_, nd) in enumerate(subS):
                self.pos[nd] = tuple(pol(tm + sg * (-(gap + bend + step * i)) / Rm, rsec))
                self.sec_n = getattr(self, 'sec_n', {}); self.sec_n.setdefault(nd, tuple(n))
            for i, (_, nd) in enumerate(prodS):
                self.pos[nd] = tuple(pol(tm + sg * (gap + bend + step * i) / Rm, rsec))
                self.sec_n = getattr(self, 'sec_n', {}); self.sec_n.setdefault(nd, tuple(n))
        else:
            if axis is None:
                u = unit(np.mean(pp, 0) - np.mean(ps, 0))
            else:
                u = unit(axis)
            if mid is None:
                mid = (np.mean(ps, 0) + np.mean(pp, 0)) / 2
            mid = np.array(mid, float) + np.array(shift, float)
            ms, mp = mid - u * gap, mid + u * gap
            if curve is None: curve = (shift != (0, 0))
            perp = np.array([-u[1], u[0]])
            if side in ('out', 'in'):
                c = np.array(center, float)
                want = unit(mid - c) * (1 if side == 'out' else -1)
            else:
                want = np.array(DIRS[side], float)
            n = perp if perp @ want >= 0 else -perp
            self.pos[I['mid']] = tuple(mid); self.pos[I['sm']] = tuple(ms); self.pos[I['pm']] = tuple(mp)
            self.rxn_geom = getattr(self, 'rxn_geom', {})
            self.rxn_geom[rid] = dict(u=u, n=n, mid=mid, ms=ms, mp=mp, curve=curve, bow=bow)
            sd = secdist or dist
            subS = [x for x in I['subs'] if not N[x[1]]['node_is_primary']]
            prodS = [x for x in I['prods'] if not N[x[1]]['node_is_primary']]
            for i, (_, nd) in enumerate(subS):
                if nd in sec: self.pos[nd] = tuple(sec[nd])
                elif nd not in self.pos:
                    self.pos[nd] = tuple(mid + u * (-gap - bend - step * i) + n * sd)
                self.sec_n = getattr(self, 'sec_n', {}); self.sec_n.setdefault(nd, tuple(n))
            for i, (_, nd) in enumerate(prodS):
                if nd in sec: self.pos[nd] = tuple(sec[nd])
                elif nd not in self.pos:
                    self.pos[nd] = tuple(mid + u * (gap + bend + step * i) + n * sd)
                self.sec_n = getattr(self, 'sec_n', {}); self.sec_n.setdefault(nd, tuple(n))
        # reaction label
        w = 0.62 * 22 * len(self.R[rid]['bigg_id'])
        if lab is None: lab = side if side in DIRS else ('R' if n[0] >= 0 else 'L')
        if lab in ('out', 'in'):
            lv = n if lab == side else -n
            lab = 'R' if lv[0] >= 0 else 'L'
        if lab == 'R': lp = (mid[0] + 24, mid[1] + 9, 'start')
        elif lab == 'L': lp = (mid[0] - 24, mid[1] + 9, 'end')
        elif lab == 'U': lp = (mid[0], mid[1] - 26, 'middle')
        else: lp = (mid[0], mid[1] + 44, 'middle')
        self.rlab[rid] = (lp[0] + labdx, lp[1] + labdy, lp[2])

    def finish(self):
        N, R = self.N, self.R
        # nodes
        for k, v in self.pos.items():
            N[k]['x'], N[k]['y'] = float(v[0]), float(v[1])
        # segments
        tcv = getattr(self, 'curv', 0.0)       # 0 = gentle, 1 = strongly curved (earlier version)
        lerp = lambda a, b: a * (1 - tcv) + b * tcv
        def qrule(M, K, u, n, sub):
            dv = M - K; d = dv @ n; dl = dv @ u
            if abs(d) < 25 or abs(dl) < 20: return None
            sg = 1 if dl > 0 else -1
            near_g = K + u * sg * min(0.30 * abs(dl), 0.35 * abs(d), 60)
            far_g = M + (K - M) * 0.33
            near_s = K + u * sg * min(0.45 * abs(dl), 0.4 * abs(d))
            far_s = K + u * (dl + sg * 0.2 * abs(d))
            c_near, c_far = lerp(near_g, near_s), lerp(far_g, far_s)
            return (c_far, c_near) if sub else (c_near, c_far)
        for rid, r in R.items():
            I = self.info[rid]; g = self.rxn_geom[rid]
            u0, n0, curve, arc = g['u'], g['n'], g['curve'], g.get('arc')
            for sid, s in r['segments'].items():
                a, b = s['from_node_id'], s['to_node_id']
                ta, tb = N[a]['node_type'], N[b]['node_type']
                s['b1'] = s['b2'] = None
                if ta == 'metabolite' and tb == 'multimarker':
                    met, mk, sub = a, b, True
                elif ta == 'multimarker' and tb == 'metabolite':
                    met, mk, sub = b, a, False
                else:
                    continue
                pm_ = np.array(self.pos[met]); pk = np.array(self.pos[mk])
                dist = np.linalg.norm(pm_ - pk)
                if dist < 1: continue
                res = None
                if arc is not None:
                    C = arc['C']; sgn = arc['s']
                    tau = lambda p: sgn * np.array([-(p - C)[1], (p - C)[0]]) / np.linalg.norm(p - C)
                    if met in arc['nodes']:
                        h = lerp(0.33, 0.36) * dist
                        if sub: res = (pm_ + tau(pm_) * h, pk - tau(pk) * h)
                        else:   res = (pk + tau(pk) * h, pm_ - tau(pm_) * h)
                    else:
                        uk = tau(pk); nk = (pk - C) / np.linalg.norm(pk - C)
                        res = qrule(pm_, pk, uk, nk, sub)
                elif N[met]['node_is_primary']:
                    if curve:
                        k = lerp(0.32, 0.42) * dist
                        res = (pm_ + u0 * k, pk - u0 * k) if sub else (pk + u0 * k, pm_ - u0 * k)
                    else:
                        res = qrule(pm_, pk, u0, n0, sub)
                        if res is None and g.get('bow') and dist > 300:
                            dd = (pk - pm_); pn = np.array([-dd[1], dd[0]]) / np.linalg.norm(dd) * g['bow'] * dist
                            res = (pm_ + dd / 3 + pn, pm_ + 2 * dd / 3 + pn)
                            if not sub: res = (pk + (pm_ - pk) / 3 + pn, pk + 2 * (pm_ - pk) / 3 + pn)
                else:
                    res = qrule(pm_, pk, u0, n0, sub)
                if res is None: continue
                s['b1'] = {'x': float(res[0][0]), 'y': float(res[0][1])}
                s['b2'] = {'x': float(res[1][0]), 'y': float(res[1][1])}
            lx, ly, anc = self.rlab[rid]
            w = 0.62 * 22 * len(r['bigg_id'])
            r['label_x'] = float(lx - (w if anc == 'end' else w / 2 if anc == 'middle' else 0))
            r['label_y'] = float(ly)
            r['_anchor'] = anc
        # metabolite labels
        for k, v in N.items():
            if v['node_type'] != 'metabolite': continue
            x, y = self.pos[k]
            if k in self.lab:
                lx, ly, anc = self.lab[k]
            elif v['node_is_primary']:
                lx, ly, anc = x + 26, y - 20, 'start'
            else:
                nn = self.sec_n.get(k, (1, 0))
                if abs(nn[0]) >= abs(nn[1]):
                    lx, ly, anc = (x + 18, y + 7, 'start') if nn[0] > 0 else (x - 18, y + 7, 'end')
                else:
                    lx, ly, anc = (x, y + 32, 'middle') if nn[1] > 0 else (x, y - 18, 'middle')
            v['label_x'], v['label_y'], v['_anchor'] = float(lx), float(ly), anc


# ---------------------------------------------------------------- label placement + geometry helpers
def short(b):
    for suf in ('_c', '_e', '_p'):
        if b.endswith(suf): return b[:-2]
    return b

def bez(p0, b1, b2, p1, n=24):
    t = np.linspace(0, 1, n)[:, None]
    p0, b1, b2, p1 = map(np.array, (p0, b1, b2, p1))
    return (1-t)**3*p0 + 3*(1-t)**2*t*b1 + 3*(1-t)*t**2*b2 + t**3*p1

def seg_points(L):
    pts = []
    N = L.N
    for rid, r in L.R.items():
        for s in r['segments'].values():
            a = np.array([N[s['from_node_id']]['x'], N[s['from_node_id']]['y']])
            b = np.array([N[s['to_node_id']]['x'], N[s['to_node_id']]['y']])
            if s['b1'] is not None:
                c1 = np.array([s['b1']['x'], s['b1']['y']]); c2 = np.array([s['b2']['x'], s['b2']['y']])
                pts.append(bez(a, c1, c2, b, 30))
            else:
                k = max(2, int(np.linalg.norm(b - a) / 15))
                pts.append(np.linspace(a, b, k))
    return np.vstack(pts)

def tbox(x, y, anc, text, size):
    w = 0.58 * size * len(text); h = size
    if anc == 'start': x0 = x
    elif anc == 'end': x0 = x - w
    else: x0 = x - w / 2
    return (x0 - 3, y - h * 0.85, x0 + w + 3, y + h * 0.28)

def overlap(a, b):
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])

def place_labels(L, sizes=(26, 20, 22)):
    N, R = L.N, L.R
    pts = seg_points(L)
    boxes = []
    # (reaction labels are placed after the scoring helpers are defined)
    # node circles as small boxes
    circ = []
    for k, v in N.items():
        if v['node_type'] == 'metabolite':
            rr = 20 if v['node_is_primary'] else 10
            circ.append((v['x'] - rr, v['y'] - rr, v['x'] + rr, v['y'] + rr))
    def score(b):
        sc = 0
        inside = ((pts[:, 0] > b[0]) & (pts[:, 0] < b[2]) & (pts[:, 1] > b[1]) & (pts[:, 1] < b[3])).sum()
        sc += 3 * inside
        for o in boxes + circ:
            if overlap(b, o): sc += 40
        return sc
    # reaction labels: free spot close to the midmarker, never on top of a line
    import math as _m
    def rbox(x, y, anc, text, size=22):
        w = 0.64 * size * len(text); h = size
        x0 = x if anc == 'start' else x - w if anc == 'end' else x - w / 2
        return (x0 - 5, y - h * 0.85 - 3, x0 + w + 5, y + h * 0.28 + 3)
    order = sorted(R.keys(), key=lambda k: -len(R[k]['bigg_id']))
    for rid in order:
        r = R[rid]; mx, my = L.rxn_geom[rid]['mid']; t = r['bigg_id']
        best = None
        for rad in (30, 46, 66, 92, 125):
            for ang in (0, 180, 90, 270, 45, 135, 225, 315):
                dx, dy = rad * _m.cos(_m.radians(ang)), rad * _m.sin(_m.radians(ang))
                if abs(dx) > 0.5 * rad:
                    cand = (mx + dx, my + dy + 8, 'start' if dx > 0 else 'end')
                else:
                    cand = (mx + dx, my + dy + (24 if dy > 0 else -6), 'middle')
                bx = rbox(cand[0], cand[1], cand[2], t)
                sc = score(bx) + 0.10 * rad
                if best is None or sc < best[0]: best = (sc, cand)
        L.rlab[rid] = best[1]
        boxes.append(rbox(best[1][0], best[1][1], best[1][2], t))
    # secondary labels (rule-based first)
    for k, v in N.items():
        if v['node_type'] != 'metabolite' or v['node_is_primary']: continue
        x, y = v['x'], v['y']; t = short(v['bigg_id'])
        nn = L.sec_n.get(k, (1, 0))
        cands = []
        if abs(nn[0]) >= abs(nn[1]):
            cands = [(x + 17, y + 7, 'start'), (x - 17, y + 7, 'end')] if nn[0] > 0 else [(x - 17, y + 7, 'end'), (x + 17, y + 7, 'start')]
            cands += [(x, y - 17, 'middle'), (x, y + 32, 'middle')]
        else:
            cands = [(x, y + 32, 'middle'), (x, y - 17, 'middle')] if nn[1] > 0 else [(x, y - 17, 'middle'), (x, y + 32, 'middle')]
            cands += [(x + 17, y + 7, 'start'), (x - 17, y + 7, 'end')]
        best = min(range(len(cands)), key=lambda i: score(tbox(*cands[i][:2], cands[i][2], t, sizes[1])) + 0.5 * i)
        c = cands[best]
        v['_lab'] = c; boxes.append(tbox(c[0], c[1], c[2], t, sizes[1]))
    # primary labels
    for k, v in N.items():
        if v['node_type'] != 'metabolite' or not v['node_is_primary']: continue
        x, y = v['x'], v['y']; t = short(v['bigg_id']); s = sizes[0]
        pref = L.lab.get(k)
        cands = [(x + 28, y + 9, 'start'), (x - 28, y + 9, 'end'), (x, y - 30, 'middle'), (x, y + 52, 'middle'),
                 (x + 22, y - 24, 'start'), (x - 22, y - 24, 'end'), (x + 22, y + 48, 'start'), (x - 22, y + 48, 'end')]
        if pref: cands = [pref] + cands
        best = min(range(len(cands)), key=lambda i: score(tbox(*cands[i][:2], cands[i][2], t, s)) + 0.4 * i)
        c = cands[best]
        v['_lab'] = c; boxes.append(tbox(c[0], c[1], c[2], t, s))
    # exit labels
    prim = [(k, np.array(L.orig[k])) for k, v in N.items() if v['node_type'] == 'metabolite' and v['node_is_primary']]
    L.exit_assign = {}
    for tid, tl in L.m['text_labels'].items():
        if not tl['text'].startswith('→'): continue
        p = np.array([tl['x'], tl['y']])
        kbest = min(prim, key=lambda kv: np.linalg.norm(kv[1] - p))[0]
        L.exit_assign[tid] = kbest
    groups = {}
    for tid, k in L.exit_assign.items(): groups.setdefault(k, []).append(tid)
    for k, tids in groups.items():
        v = N[k]; x, y = v['x'], v['y']
        for j, tid in enumerate(sorted(tids)):
            t = L.m['text_labels'][tid]['text']; s = sizes[2]
            cands = [(x + 28, y + 40 + 26 * j, 'start'), (x - 28, y + 40 + 26 * j, 'end'),
                     (x + 28, y - 34 - 26 * j, 'start'), (x - 28, y - 34 - 26 * j, 'end'),
                     (x, y + 80 + 26 * j, 'middle'), (x, y - 50 - 26 * j, 'middle'),
                     (x + 28, y + 9, 'start'), (x - 28, y + 9, 'end')]
            best = min(range(len(cands)), key=lambda i: score(tbox(*cands[i][:2], cands[i][2], t, s)) + 0.6 * i)
            c = cands[best]
            L.m['text_labels'][tid]['_lab'] = c
            boxes.append(tbox(c[0], c[1], c[2], t, s))

def zone_boxes(L):
    N = L.N; zb = {}
    for rid, z in L.rxn_zone.items():
        I = L.info[rid]; ids = {I['mid'], I['sm'], I['pm']} | {n for _, n in I['subs']} | {n for _, n in I['prods']}
        zb.setdefault(z, set()).update(ids)
    out = {}
    for z, ids in zb.items():
        xs = [N[i]['x'] for i in ids]; ys = [N[i]['y'] for i in ids]
        out[z] = [min(xs), min(ys), max(xs), max(ys)]
    return out
