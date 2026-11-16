"""Pictures to look at: the whole map or the surroundings of a metabolite or reaction, with the defects marked.

`run.py view` renders the layout of a spec without writing a map: the overview, or one picture per `--around` name at
a scale where every label can be read. With marks, the places that the metrics count are drawn into the picture: a red
cross on a crossing of two reactions, an orange circle around a crowded reaction, a magenta frame around a label that
lies on a line. So the numbers of a build and what the eye sees point at the same spots.
"""
import os
import numpy as np
from .engine import Layout, place_labels
from . import pack, qa, render as rendermod
from .pipeline import apply_spec, check_complete, content_bbox


def laid_out(input_path, spec, curv=0.6, do_pack=True):
    """The finished Layout of a spec (as `build` makes it, without writing files)."""
    L = Layout(input_path); apply_spec(L, spec); L.curv = curv
    check_complete(L); L.finish(); place_labels(L)
    if do_pack: pack.pack(L)
    return L


def _cross_point(P, Q):
    """A point where polyline P crosses polyline Q (None if they do not cross)."""
    p1, p2 = P[:-1][:, None, :], P[1:][:, None, :]; q1, q2 = Q[:-1][None, :, :], Q[1:][None, :, :]
    cr = lambda a, b, c: (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
    hit = np.argwhere(((cr(p1, p2, q1) * cr(p1, p2, q2)) < 0) & ((cr(q1, q2, p1) * cr(q1, q2, p2)) < 0))
    if not len(hit): return None
    i, j = hit[0]
    return (P[i] + P[i + 1] + Q[j] + Q[j + 1]) / 4


def defects(L):
    """Where the metrics count something: dict(crossings=[(x, y, reaction, reaction)], crowded=[(x, y, reaction)],
    labels=[(box, name)]) for a finished Layout. The same rules as :func:`qa.metrics`."""
    N, R = L.N, L.R; pls = qa._polylines(L); trimmed = []
    for rid, a, b, pts in pls:
        t = pts[(np.linalg.norm(pts - a, axis=1) > 40) & (np.linalg.norm(pts - b, axis=1) > 40)]
        trimmed.append((rid, t if len(t) >= 2 else None))
    boxes = [(t.min(0), t.max(0)) if t is not None else None for _, t in trimmed]; crossings = []
    for i in range(len(trimmed)):
        if trimmed[i][1] is None: continue
        for j in range(i + 1, len(trimmed)):
            if trimmed[j][1] is None or trimmed[i][0] == trimmed[j][0]: continue
            (a0, a1), (b0, b1) = boxes[i], boxes[j]
            if a1[0] < b0[0] or b1[0] < a0[0] or a1[1] < b0[1] or b1[1] < a0[1]: continue
            p = _cross_point(trimmed[i][1], trimmed[j][1])
            if p is not None: crossings.append((float(p[0]), float(p[1]), R[trimmed[i][0]]['bigg_id'], R[trimmed[j][0]]['bigg_id']))
    m = qa.metrics(L); crowd = set(m.get('crowded_list', []))
    crowded = [(N[I['mid']]['x'], N[I['mid']]['y'], R[rid]['bigg_id']) for rid, I in L.info.items() if R[rid]['bigg_id'] in crowd]
    allpts = np.vstack([p for _, _, _, p in pls]); labels = []
    for rid, r in R.items():
        x, y, anc = L.rlab[rid]; labels.append((qa._box(x, y, anc, r['bigg_id'], 22, 0.64), 'reaction ' + r['bigg_id']))
    for k, v in N.items():
        if v['node_type'] == 'metabolite' and '_lab' in v:
            x, y, anc = v['_lab']; t = qa.short(v['bigg_id'])
            labels.append((qa._box(x, y, anc, t, 26 if v['node_is_primary'] else 20, 0.58), 'metabolite ' + t))
    on_line = [(b, name) for b, name in labels
               if ((allpts[:, 0] > b[0]) & (allpts[:, 0] < b[2]) & (allpts[:, 1] > b[1]) & (allpts[:, 1] < b[3])).sum() >= 2]
    return dict(crossings=crossings, crowded=crowded, labels=on_line)


def _mark(ax, d, scale):
    from matplotlib.patches import Circle, Rectangle
    w = max(1.5, 14 * scale)
    for x, y, _, _ in d['crossings']:
        ax.plot([x - 28, x + 28], [y - 28, y + 28], color='#e00000', lw=w, zorder=20, solid_capstyle='round')
        ax.plot([x - 28, x + 28], [y + 28, y - 28], color='#e00000', lw=w, zorder=20, solid_capstyle='round')
    for x, y, _ in d['crowded']:
        ax.add_patch(Circle((x, y), 60, fc='none', ec='#ff8c00', lw=w, zorder=19))
    for b, _ in d['labels']:
        ax.add_patch(Rectangle((b[0] - 6, b[1] - 6), b[2] - b[0] + 12, b[3] - b[1] + 12, fc='none', ec='#d000d0', lw=w * 0.7, zorder=19))


def _render(L, path, bbox, scale, d):
    """rendermod.render with the marks of `d` drawn in (None: no marks)."""
    import matplotlib.pyplot as plt
    if d is None: return rendermod.render(L, path, bbox=bbox, scale=scale)
    save = plt.Figure.savefig
    def savefig(fig, *a, **k):                # draw the marks into the picture before it is written
        _mark(fig.axes[0], d, scale); return save(fig, *a, **k)
    plt.Figure.savefig = savefig
    try: rendermod.render(L, path, bbox=bbox, scale=scale)
    finally: plt.Figure.savefig = save


def find(L, name):
    """The point a name stands for: a node id, the BiGG id of a metabolite (its primary node with the most
    reactions) or of a reaction (its middle). ValueError if nothing has the name."""
    N, R = L.N, L.R
    if name in N and N[name]['node_type'] == 'metabolite': return N[name]['x'], N[name]['y']
    uses = {}
    for rid, I in L.info.items():
        for _, n in I['subs'] + I['prods']: uses[n] = uses.get(n, 0) + 1
    cand = [k for k, v in N.items() if v['node_type'] == 'metabolite' and v['node_is_primary'] and v.get('bigg_id') in (name, name + '_c')]
    if cand:
        k = max(cand, key=lambda k: (uses.get(k, 0), -int(k))); return N[k]['x'], N[k]['y']
    rs = [rid for rid, r in R.items() if r['bigg_id'] == name]
    if rs: return N[L.info[rs[0]]['mid']]['x'], N[L.info[rs[0]]['mid']]['y']
    raise ValueError(f'nothing is called {name!r}: give the BiGG id of a metabolite (pyr_c) or of a reaction (PYK), or a node id')


def view(input_path, spec, out_dir, around=(), radius=900.0, marks=True, curv=0.6, name='view'):
    """Write the pictures and return [(path, what it shows, the defects in it)]: the overview if `around` is empty,
    else one picture of 2 x `radius` map units around each name."""
    os.makedirs(out_dir, exist_ok=True)
    L = laid_out(input_path, spec, curv); d = defects(L); out = []
    inside = lambda box: dict(crossings=[c for c in d['crossings'] if box[0] <= c[0] <= box[2] and box[1] <= c[1] <= box[3]],
                              crowded=[c for c in d['crowded'] if box[0] <= c[0] <= box[2] and box[1] <= c[1] <= box[3]],
                              labels=[c for c in d['labels'] if box[0] <= c[0][0] <= box[2] and box[1] <= c[0][1] <= box[3]])
    if not around:
        bb = content_bbox(L); box = (bb[0] - 150, bb[1] - 150, bb[2] + 150, bb[3] + 150)
        path = os.path.join(out_dir, f'{name}_overview.png'); _render(L, path, box, min(0.4, 3200.0 / (box[2] - box[0])), d if marks else None)
        out.append((path, 'the whole map', inside(box)))
    for a in around:
        x, y = find(L, a); box = (x - radius, y - radius, x + radius, y + radius)
        safe = ''.join(ch if ch.isalnum() else '_' for ch in a)
        path = os.path.join(out_dir, f'{name}_{safe}.png'); _render(L, path, box, min(1.0, 1600.0 / (2 * radius)), d if marks else None)
        out.append((path, f'{2 * radius:.0f} x {2 * radius:.0f} map units around {a}', inside(box)))
    return out


def describe(found):
    """The text of what a picture's marks mean."""
    c, w, l = found['crossings'], found['crowded'], found['labels']
    pairs = {}
    for _, _, a, b in c: pairs[tuple(sorted((a, b)))] = pairs.get(tuple(sorted((a, b))), 0) + 1
    lines = [f"{len(c)} crossings (red crosses)" + (': ' + ', '.join(f'{a} x {b}' + (f' ({n})' if n > 1 else '') for (a, b), n in sorted(pairs.items(), key=lambda kv: -kv[1])[:10]) if c else ''),
             f"{len(w)} crowded reactions (orange circles)" + (': ' + ', '.join(sorted(x[2] for x in w)[:14]) if w else ''),
             f"{len(l)} labels on lines (magenta frames)" + (': ' + ', '.join(n for _, n in l[:10]) if l else '')]
    return lines
