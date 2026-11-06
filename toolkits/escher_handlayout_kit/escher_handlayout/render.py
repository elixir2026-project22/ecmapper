"""Preview renderer (matplotlib) and Escher JSON export."""
import json, copy
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle, Polygon
import matplotlib.patheffects as pe
from .engine import unit, short
from ._io import write_json

def arrow(ax, tip, d, color, ln=34, wd=24):
    d = unit(d); p = np.array([-d[1], d[0]]); base = np.array(tip) - d * ln
    ax.add_patch(Polygon([tip, base + p * wd / 2, base - p * wd / 2], closed=True, fc=color, ec=color, lw=0.3, zorder=5))

def hub_set(L):
    cnt = {}
    for rid, I in L.info.items():
        for _, n in I['subs'] + I['prods']:
            if L.N[n]['node_is_primary']: cnt.setdefault(n, set()).add(L.rxn_zone[rid])
    return {n for n, z in cnt.items() if len(z) > 1}

def render(L, path, bbox=None, scale=0.3, titles=True):
    N, R = L.N, L.R
    xs = [v['x'] for v in N.values()]; ys = [v['y'] for v in N.values()]
    if bbox is None: bbox = (min(xs) - 300, min(ys) - 450, max(xs) + 500, max(ys) + 300)
    x0, y0, x1, y1 = bbox
    fig = plt.figure(figsize=((x1 - x0) * scale / 100, (y1 - y0) * scale / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    pt = lambda u: u * scale * 0.72
    hubs = hub_set(L)
    for rid, r in R.items():
        col = L.zone_style[L.rxn_zone[rid]]['color']; rev = r['reversibility']
        for sid, s in r['segments'].items():
            a, b = N[s['from_node_id']], N[s['to_node_id']]
            A = np.array([a['x'], a['y']]); B = np.array([b['x'], b['y']])
            if s['b1'] is not None:
                c1 = np.array([s['b1']['x'], s['b1']['y']]); c2 = np.array([s['b2']['x'], s['b2']['y']])
                verts = [A, c1, c2, B]; codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            else:
                c1 = c2 = None; verts = [A, B]; codes = [Path.MOVETO, Path.LINETO]
            ax.add_patch(PathPatch(Path(verts, codes), fc='none', ec=col, lw=pt(7), capstyle='round', zorder=2))
            ta, tb = a['node_type'], b['node_type']
            if ta == 'multimarker' and tb == 'metabolite':
                rr = (20 if b['node_is_primary'] else 10) + 4
                dv = B - (c2 if c2 is not None else A); arrow(ax, B - unit(dv) * rr, dv, col)
            if rev and ta == 'metabolite' and tb == 'multimarker':
                rr = (20 if a['node_is_primary'] else 10) + 4
                dv = A - (c1 if c1 is not None else B); arrow(ax, A - unit(dv) * rr, dv, col)
    for k, v in N.items():
        if v['node_type'] == 'metabolite':
            if v['node_is_primary']:
                hub = k in hubs
                ax.add_patch(Circle((v['x'], v['y']), 26 if hub else 20, fc='#ffd27a' if hub else '#f9c9a0',
                             ec='#8a4b12', lw=pt(5 if hub else 3), zorder=6))
            else:
                ax.add_patch(Circle((v['x'], v['y']), 10, fc='#d7dde6', ec='#7f8a99', lw=pt(3), zorder=6))
        elif v['node_type'] == 'midmarker':
            ax.add_patch(Circle((v['x'], v['y']), 7, fc='#333', ec='#333', zorder=6))
    for k, v in N.items():
        if v['node_type'] == 'metabolite' and '_lab' in v:
            x, y, anc = v['_lab']; s = 26 if v['node_is_primary'] else 20
            ax.text(x, y, short(v['bigg_id']), fontsize=pt(s * (1.15 if k in hubs else 1)), ha={'start': 'left', 'end': 'right', 'middle': 'center'}[anc],
                    va='baseline', color='#222' if v['node_is_primary'] else '#555', fontweight='bold' if v['node_is_primary'] else 'normal', zorder=8,
                    path_effects=[pe.withStroke(linewidth=pt(5), foreground='white')])
    for rid, r in R.items():
        x, y, anc = L.rlab[rid]
        ax.text(x, y, r['bigg_id'], fontsize=pt(22), ha={'start': 'left', 'end': 'right', 'middle': 'center'}[anc], va='baseline',
                color=L.zone_style[L.rxn_zone[rid]]['color'], fontstyle='italic', fontweight='bold', zorder=8,
                path_effects=[pe.withStroke(linewidth=pt(6), foreground='white')])
    for tid, tl in L.m['text_labels'].items():
        if '_lab' in tl:
            x, y, anc = tl['_lab']
            ax.text(x, y, tl['text'], fontsize=pt(22), ha={'start': 'left', 'end': 'right', 'middle': 'center'}[anc], va='baseline',
                    color='#6b7a8f', fontstyle='italic', zorder=8, path_effects=[pe.withStroke(linewidth=pt(5), foreground='white')])
    if titles:
        for z, (x, y) in L.titles.items():
            ax.text(x, y, L.zone_style[z]['title'], fontsize=pt(48), fontweight='bold', color=L.zone_style[z]['color'], ha='left', va='baseline', zorder=9, alpha=0.95)
        for text, z, x, y, _ in getattr(L, 'extra_titles', []):
            ax.text(x, y, text, fontsize=pt(40), fontweight='bold', color=L.zone_style[z]['color'], ha='left', va='baseline', zorder=9, alpha=0.95)
    fig.savefig(path, dpi=100); plt.close(fig)

def export_json(L, path):
    d = copy.deepcopy(L.d); meta, m = d; N, R = m['nodes'], m['reactions']
    for k, v in N.items():
        if v['node_type'] == 'metabolite' and '_lab' in v:
            x, y, anc = v['_lab']; t = short(v['bigg_id']); s = 26 if v['node_is_primary'] else 20; w = 0.58 * s * len(t)
            v['label_x'] = float(x - (w if anc == 'end' else w / 2 if anc == 'middle' else 0)); v['label_y'] = float(y)
        v.pop('_lab', None); v.pop('_anchor', None)
    for rid, r in R.items():
        r.pop('_anchor', None)
        x, y, anc = L.rlab[rid]; w = 0.62 * 22 * len(r['bigg_id'])
        r['label_x'] = float(x - (w if anc == 'end' else w / 2 if anc == 'middle' else 0)); r['label_y'] = float(y)
    TL = m['text_labels']; ids = {tl['text']: tid for tid, tl in TL.items() if not tl['text'].startswith('→')}
    for z, (x, y) in L.titles.items():
        st = L.zone_style[z]; lab = st.get('label')
        if lab and lab in ids:
            TL[ids[lab]].update(text=lab, x=float(x), y=float(y))
        else:
            nid = str(max(int(k) for k in TL) + 1); TL[nid] = dict(text=st['title'], x=float(x), y=float(y))
    touched = {ids[L.zone_style[z]['label']] for z in L.titles if L.zone_style[z].get('label') in ids}
    for text, z, x, y, reuse in getattr(L, 'extra_titles', []):
        if reuse and reuse in ids:
            TL[ids[reuse]].update(text=text, x=float(x), y=float(y)); touched.add(ids[reuse])
        else:
            nid = str(max(int(k) for k in TL) + 1); TL[nid] = dict(text=text, x=float(x), y=float(y)); touched.add(nid)
    # leftover title labels of the input map that this layout does not use would float far outside the canvas: drop them
    L.removed_labels = []
    for tid in list(TL):
        if tid not in touched and '_lab' not in TL[tid] and not TL[tid]['text'].startswith('→'):
            L.removed_labels.append((tid, TL[tid]['text'])); del TL[tid]
    for tid, tl in TL.items():
        if '_lab' in tl:
            x, y, anc = tl['_lab']; w = 0.58 * 22 * len(tl['text'])
            tl['x'] = float(x - (w if anc == 'end' else w / 2 if anc == 'middle' else 0)); tl['y'] = float(y); del tl['_lab']
    cb = getattr(L, 'canvas_bbox', None)
    if cb is None:
        xs = [v['x'] for v in N.values()]; ys = [v['y'] for v in N.values()]
        cb = (min(xs) - 400, min(ys) - 600, max(xs) + 700, max(ys) + 700)
    m['canvas'] = dict(x=cb[0], y=cb[1], width=cb[2] - cb[0], height=cb[3] - cb[1])
    m['reactions'] = R
    write_json(path, d)
    return d

