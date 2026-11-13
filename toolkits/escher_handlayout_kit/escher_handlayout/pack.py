"""Pack disconnected sub-graphs into free space so the overall rectangle is minimal."""
import numpy as np, math
from .engine import bez, short
from .inspect_map import base, COFACTORS

CELL = 40
def tb(x, y, anc, text, size, f=0.6):
    w = f * size * len(text); h = size
    x0 = x if anc == 'start' else x - w if anc == 'end' else x - w / 2
    return (x0, y - h * 0.85, x0 + w, y + h * 0.28)

def components(L):
    par = {}
    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for rid, I in L.info.items():
        for n in {I['mid'], I['sm'], I['pm']} | {n for _, n in I['subs'] + I['prods']}:
            par[f('n' + n)] = f('r' + rid)
    comps = {}
    for rid in L.R: comps.setdefault(f('r' + rid), []).append(rid)
    out = []
    for rids in comps.values():
        nodes = set()
        for rid in rids:
            I = L.info[rid]; nodes |= {I['mid'], I['sm'], I['pm']} | {n for _, n in I['subs'] + I['prods']}
        out.append((rids, nodes))
    return sorted(out, key=lambda c: -len(c[0]))

def elements(L, rids, nodes, extra_boxes=()):
    """points (lines + node circles) and boxes (labels) of a set of reactions/nodes in global coords"""
    N, R = L.N, L.R; pts = []; boxes = list(extra_boxes)
    for rid in rids:
        for s in R[rid]['segments'].values():
            a = np.array([N[s['from_node_id']]['x'], N[s['from_node_id']]['y']]); b = np.array([N[s['to_node_id']]['x'], N[s['to_node_id']]['y']])
            if s['b1'] is not None:
                pts.append(bez(a, [s['b1']['x'], s['b1']['y']], [s['b2']['x'], s['b2']['y']], b, 40))
            else:
                pts.append(np.linspace(a, b, max(2, int(np.linalg.norm(b - a) / 10))))
        x, y, anc = L.rlab[rid]; boxes.append(tb(x, y, anc, R[rid]['bigg_id'], 22, 0.64))
    for k in nodes:
        v = N[k]
        if v['node_type'] == 'metabolite':
            r = 26 if v['node_is_primary'] else 12
            th = np.linspace(0, 2 * math.pi, 12, endpoint=False)
            pts.append(np.c_[v['x'] + r * np.cos(th), v['y'] + r * np.sin(th)]); pts.append(np.array([[v['x'], v['y']]]))
            if '_lab' in v:
                x, y, anc = v['_lab']; boxes.append(tb(x, y, anc, short(v['bigg_id']), 26 if v['node_is_primary'] else 20))
    for tid, k in L.exit_assign.items():
        if k in nodes:
            tl = L.m['text_labels'][tid]; x, y, anc = tl['_lab']; boxes.append(tb(x, y, anc, tl['text'], 22))
    return (np.vstack(pts) if pts else np.zeros((0, 2))), boxes

def bbox_of(pts, boxes):
    xs = list(pts[:, 0]) + [b[0] for b in boxes] + [b[2] for b in boxes]
    ys = list(pts[:, 1]) + [b[1] for b in boxes] + [b[3] for b in boxes]
    return min(xs), min(ys), max(xs), max(ys)

def raster(pts, boxes, x0, y0, H, W, dil=0):
    g = np.zeros((H, W), bool)
    def mark(i0, i1, j0, j1):
        i0, j0, i1, j1 = max(i0 - dil, 0), max(j0 - dil, 0), min(i1 + dil + 1, H), min(j1 + dil + 1, W)
        if i1 > i0 and j1 > j0: g[i0:i1, j0:j1] = True
    for x, y in pts:
        i = int((y - y0) // CELL); j = int((x - x0) // CELL); mark(i, i, j, j)
    for b in boxes:
        mark(int((b[1] - y0) // CELL), int((b[3] - y0) // CELL), int((b[0] - x0) // CELL), int((b[2] - x0) // CELL))
    return g

def translate(L, rids, nodes, dx, dy):
    N, R = L.N, L.R
    for k in nodes:
        v = N[k]; v['x'] += dx; v['y'] += dy
        if '_lab' in v: x, y, a = v['_lab']; v['_lab'] = (x + dx, y + dy, a)
        L.pos[k] = (v['x'], v['y'])
    for rid in rids:
        for s in R[rid]['segments'].values():
            if s['b1'] is not None:
                for b in ('b1', 'b2'): s[b]['x'] += dx; s[b]['y'] += dy
        x, y, a = L.rlab[rid]; L.rlab[rid] = (x + dx, y + dy, a)
        g = L.rxn_geom[rid]
    for tid, k in L.exit_assign.items():
        if k in nodes:
            tl = L.m['text_labels'][tid]; x, y, a = tl['_lab']; tl['_lab'] = (x + dx, y + dy, a)

def pack(L, margin_cells=3):
    """Moves the small disconnected pieces into free space around the main component. The largest piece of a subsystem
    that has no reaction in the main component is a subsystem block of its own: it stays where the spec puts it (with
    the subsystem's title) unless it overlaps something there, so that the arrangement of the blocks survives while
    the subsystems are not connected yet."""
    comps = components(L); main_r, main_n = comps[0]; small = comps[1:]
    L.extra_titles = []
    zones_main = {L.rxn_zone[r] for r in main_r}
    def zone_of(rids): return L.rxn_zone[sorted(rids, key=int)[0]]
    def title_box(z): return [tb(L.titles[z][0], L.titles[z][1], 'start', L.zone_style[z]['title'], 48, 0.62)] if z in L.titles else []
    def carbon(nodes): return {L.N[n]['bigg_id'] for n in nodes if L.N[n]['node_type'] == 'metabolite' and L.N[n]['node_is_primary']
                               and base(L.N[n]['bigg_id']) not in COFACTORS}
    cands = {}                                         # zone outside the main component -> its largest piece
    for c in small:
        z = zone_of(c[0])
        if z not in zones_main and all(L.rxn_zone[r] == z for r in c[0]) and carbon(c[1]) & set().union(*[carbon(o[1]) for o in comps if o is not c]):
            cands.setdefault(z, c)                     # (a piece that shares only cofactors with the rest is packed)
    tboxes = [b for z in L.titles if z in zones_main for b in title_box(z)]
    mp, mb = elements(L, main_r, main_n, tboxes)
    mbb = bbox_of(mp, mb)
    cel = {z: elements(L, c[0], c[1], title_box(z)) for z, c in cands.items()}
    cbb = {z: bbox_of(*e) for z, e in cel.items()}
    pad = 3200
    far = lambda bb: bb[0] < mbb[0] - 4 * pad or bb[1] < mbb[1] - 4 * pad or bb[2] > mbb[2] + 4 * pad or bb[3] > mbb[3] + 4 * pad
    for z in [z for z in cands if far(cbb[z])]: cands.pop(z)            # far away from everything: pack it
    gbb = [min([mbb[i]] + [cbb[z][i] for z in cands]) for i in (0, 1)] + [max([mbb[i]] + [cbb[z][i] for z in cands]) for i in (2, 3)]
    gx0, gy0 = gbb[0] - pad, gbb[1] - pad
    W = int((gbb[2] - gbb[0] + 2 * pad) // CELL) + 2; H = int((gbb[3] - gbb[1] + 2 * pad) // CELL) + 2
    occ = raster(mp, mb, gx0, gy0, H, W, dil=margin_cells)
    cur = list(mbb)
    fixed = []
    for z, c in list(cands.items()):
        own = raster(*cel[z], gx0, gy0, H, W, dil=0)
        if (own & occ).any(): continue                 # it overlaps what is already there: pack it like a small piece
        occ |= raster(*cel[z], gx0, gy0, H, W, dil=margin_cells); fixed.append(c)
        cur = [min(cur[0], cbb[z][0]), min(cur[1], cbb[z][1]), max(cur[2], cbb[z][2]), max(cur[3], cbb[z][3])]
        zones_main.add(z)
    small = [c for c in small if not any(c is f for f in fixed)]
    L.fixed_pieces = fixed
    # a zone without a reaction in the main component or in a block that stays hands its title over to a small piece
    for z in [z for z in L.titles if z not in zones_main]: L.titles.pop(z)
    info = []
    # prepare small comps with their titles
    items = []
    for rids, nodes in small:
        key = next((L.R[r]['bigg_id'] for r in sorted(rids, key=int) if L.R[r]['bigg_id'] in L.small_titles), None)
        if key is not None:
            zone, text, reuse = L.small_titles[key]
        else:
            first = sorted(rids, key=int)[0]; zone = L.rxn_zone[first]
            names = ', '.join(dict.fromkeys(L.R[r]['bigg_id'] for r in sorted(rids, key=int)[:3]))
            text = f"{L.zone_style[zone]['title']}: {names}"; reuse = None
        pts, boxes = elements(L, rids, nodes)
        bb = bbox_of(pts, boxes)
        tbx = tb(bb[0], bb[1] - 55, 'start', text, 40, 0.62)
        boxes.append(tbx); bb = bbox_of(pts, boxes)
        items.append(dict(rids=rids, nodes=nodes, pts=pts, boxes=boxes, bb=bb, zone=zone, text=text, reuse=reuse,
                          title_xy=(bb[0], bb[1] - 55 + 0), orig_c=((bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2)))
        items[-1]['title_xy'] = (tbx[0], tbx[1] + 40 * 0.85)
    items.sort(key=lambda it: -(it['bb'][2] - it['bb'][0]) * (it['bb'][3] - it['bb'][1]))
    for it in items:
        bb = it['bb']; fx0, fy0 = bb[0], bb[1]
        h = int((bb[3] - bb[1]) // CELL) + 2; w = int((bb[2] - bb[0]) // CELL) + 2
        fp = raster(it['pts'], it['boxes'], fx0, fy0, h, w, dil=0)
        # overlap count via FFT correlation
        A = np.zeros((H + h, W + w)); A[:H, :W] = occ
        B = np.zeros_like(A); B[:h, :w] = fp
        corr = np.fft.irfft2(np.fft.rfft2(A) * np.conj(np.fft.rfft2(B)), s=A.shape)
        if h > H or w > W: ii = jj = np.zeros(0, int)            # piece larger than the whole grid
        else: ii, jj = np.nonzero(corr[:H - h + 1, :W - w + 1] < 0.5)
        if len(ii) == 0:
            # no overlap-free position inside the grid (e.g. the piece is huge): put it right of everything placed so far
            x0 = np.array([cur[2] + 300.0]); y0 = np.array([cur[1]]); best = 0
        else:
            x0 = gx0 + jj * CELL; y0 = gy0 + ii * CELL
        x1 = x0 + (bb[2] - bb[0]); y1 = y0 + (bb[3] - bb[1])
        if len(ii) > 0:
            ux0 = np.minimum(x0, cur[0]); uy0 = np.minimum(y0, cur[1]); ux1 = np.maximum(x1, cur[2]); uy1 = np.maximum(y1, cur[3])
            area = (ux1 - ux0) * (uy1 - uy0)
            amin = area.min()
            cand = np.nonzero(area <= amin * 1.0005)[0]
            cx = (x0 + x1) / 2; cy = (y0 + y1) / 2
            dist = np.hypot(cx[cand] - it['orig_c'][0], cy[cand] - it['orig_c'][1])
            best = cand[np.argmin(dist)]
        dx = x0[best] - bb[0]; dy = y0[best] - bb[1]
        translate(L, it['rids'], it['nodes'], dx, dy)
        tx, ty = it['title_xy']; it['title_xy'] = (tx + dx, ty + dy)
        L.extra_titles.append((it['text'], it['zone'], it['title_xy'][0], it['title_xy'][1], it['reuse']))
        # update occupancy and bbox
        pts2 = it['pts'] + np.array([dx, dy]); boxes2 = [(b[0] + dx, b[1] + dy, b[2] + dx, b[3] + dy) for b in it['boxes']]
        occ |= raster(pts2, boxes2, gx0, gy0, H, W, dil=margin_cells)
        cur = [min(cur[0], x0[best]), min(cur[1], y0[best]), max(cur[2], x1[best]), max(cur[3], y1[best])]
        info.append((it['text'], round(dx), round(dy)))
    return mbb, cur, info
