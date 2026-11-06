"""Structural validation of the output map and numeric layout-quality metrics.

The metrics replace visual inspection: they tell an automated run whether a change of the spec made
the map better or worse (fewer line crossings, fewer labels sitting on lines, ...).
"""
import json, math
import numpy as np
from .engine import bez, short
from ._io import read_json


# ------------------------------------------------------------------------------------------ validation
def validate(input_path, output_path):
    """Return (errors, warnings). errors == [] means the output is structurally valid and keeps the biochemistry."""
    errors, warns = [], []
    a = read_json(input_path); b = read_json(output_path)
    ma, mb = a[1], b[1]
    Na, Ra = ma['nodes'], ma['reactions']; Nb, Rb = mb['nodes'], mb['reactions']
    fin = lambda v: isinstance(v, (int, float)) and math.isfinite(v)

    if set(Ra) != set(Rb): errors.append(f'reaction id set differs: missing {sorted(set(Ra)-set(Rb))}, new {sorted(set(Rb)-set(Ra))}')
    for rid in set(Ra) & set(Rb):
        ra, rb = Ra[rid], Rb[rid]
        for key in ('bigg_id', 'reversibility', 'name', 'gene_reaction_rule'):
            if ra.get(key) != rb.get(key): errors.append(f'reaction {rid}: field {key} changed')
        sa = sorted((m['bigg_id'], m['coefficient']) for m in ra.get('metabolites', []))
        sb = sorted((m['bigg_id'], m['coefficient']) for m in rb.get('metabolites', []))
        if sa != sb: errors.append(f'reaction {rid}: stoichiometry changed')
        mids = [n for s in rb['segments'].values() for n in (s['from_node_id'], s['to_node_id']) if n in Nb and Nb[n]['node_type'] == 'midmarker']
        if len(set(mids)) != 1: errors.append(f'reaction {rid}: expected exactly 1 midmarker, found {len(set(mids))}')
        # metabolites drawn must be the same multiset of bigg_ids as in the input
        def drawn(R, N):
            out = []
            for s in R['segments'].values():
                for n in (s['from_node_id'], s['to_node_id']):
                    if N[n]['node_type'] == 'metabolite': out.append(n)
            return sorted(N[n]['bigg_id'] for n in set(out))
        if drawn(ra, Na) != drawn(rb, Nb): errors.append(f'reaction {rid}: drawn metabolites differ from input')
        for sid, s in rb['segments'].items():
            if s['from_node_id'] not in Nb or s['to_node_id'] not in Nb: errors.append(f'reaction {rid} segment {sid}: dangling node reference')
            for k in ('b1', 'b2'):
                if s[k] is not None and not (fin(s[k]['x']) and fin(s[k]['y'])): errors.append(f'reaction {rid} segment {sid}: bad {k}')
        if not (fin(rb.get('label_x')) and fin(rb.get('label_y'))): errors.append(f'reaction {rid}: bad label position')
        for k in rb:
            if k.startswith('_'): errors.append(f'reaction {rid}: temporary key {k}')
    c = mb['canvas']
    for nid, v in Nb.items():
        for k in ('x', 'y'):
            if not fin(v[k]): errors.append(f'node {nid}: non-finite {k}')
        for k in v:
            if k.startswith('_'): errors.append(f'node {nid}: temporary key {k}')
        if v['node_type'] == 'metabolite' and not (fin(v.get('label_x')) and fin(v.get('label_y'))): errors.append(f'node {nid}: bad label position')
        if not (c['x'] <= v['x'] <= c['x'] + c['width'] and c['y'] <= v['y'] <= c['y'] + c['height']): errors.append(f'node {nid}: outside canvas')
    for tid, t in mb['text_labels'].items():
        if not (fin(t['x']) and fin(t['y'])): errors.append(f'text label {tid}: bad position')
        if any(k.startswith('_') for k in t): errors.append(f'text label {tid}: temporary key')
        if not (c['x'] <= t['x'] <= c['x'] + c['width'] and c['y'] <= t['y'] <= c['y'] + c['height']): warns.append(f'text label {tid} outside canvas')
    used = {n for r in Rb.values() for s in r['segments'].values() for n in (s['from_node_id'], s['to_node_id'])}
    orphan = [n for n in Nb if n not in used]
    if orphan: warns.append(f'{len(orphan)} nodes are not used by any reaction: {orphan[:10]}')
    return errors, warns


# ------------------------------------------------------------------------------------------ metrics
def _polylines(L, step=12):
    out = []
    N = L.N
    for rid, r in L.R.items():
        for sid, s in r['segments'].items():
            a = np.array([N[s['from_node_id']]['x'], N[s['from_node_id']]['y']]); b = np.array([N[s['to_node_id']]['x'], N[s['to_node_id']]['y']])
            if s['b1'] is not None:
                pts = bez(a, [s['b1']['x'], s['b1']['y']], [s['b2']['x'], s['b2']['y']], b, max(10, int(np.linalg.norm(b - a) / step)))
            else:
                pts = np.linspace(a, b, max(2, int(np.linalg.norm(b - a) / step)))
            out.append((rid, a, b, pts))
    return out

def _inter(P, Q):
    """True if polyline P intersects polyline Q (vectorised segment-segment test)."""
    p1, p2 = P[:-1][:, None, :], P[1:][:, None, :]; q1, q2 = Q[:-1][None, :, :], Q[1:][None, :, :]
    def cr(a, b, c): return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
    d1, d2, d3, d4 = cr(p1, p2, q1), cr(p1, p2, q2), cr(q1, q2, p1), cr(q1, q2, p2)
    return bool((((d1 * d2) < 0) & ((d3 * d4) < 0)).any())

def _box(x, y, anc, text, size, f):
    w = f * size * len(text); h = size
    x0 = x if anc == 'start' else x - w if anc == 'end' else x - w / 2
    return (x0, y - h * 0.85, x0 + w, y + h * 0.28)

def metrics(L):
    """Numeric layout-quality metrics of a finished Layout (after place_labels / pack)."""
    N, R = L.N, L.R
    pls = _polylines(L)
    # --- crossings between different reactions (ignoring the 40 units next to nodes)
    trimmed = []
    for rid, a, b, pts in pls:
        keep = (np.linalg.norm(pts - a, axis=1) > 40) & (np.linalg.norm(pts - b, axis=1) > 40)
        t = pts[keep]
        trimmed.append((rid, t if len(t) >= 2 else None))
    cross = {}; ncross = 0
    boxes = [(t.min(0), t.max(0)) if t is not None else None for _, t in trimmed]
    for i in range(len(trimmed)):
        if trimmed[i][1] is None: continue
        for j in range(i + 1, len(trimmed)):
            if trimmed[j][1] is None or trimmed[i][0] == trimmed[j][0]: continue
            (a0, a1), (b0, b1) = boxes[i], boxes[j]
            if a1[0] < b0[0] or b1[0] < a0[0] or a1[1] < b0[1] or b1[1] < a0[1]: continue
            if _inter(trimmed[i][1], trimmed[j][1]):
                ncross += 1; key = tuple(sorted((R[trimmed[i][0]]['bigg_id'], R[trimmed[j][0]]['bigg_id']))); cross[key] = cross.get(key, 0) + 1
    worst = sorted(cross.items(), key=lambda kv: -kv[1])[:12]
    # --- labels
    allpts = np.vstack([p for _, _, _, p in pls])
    circ = [(k, v['x'], v['y'], 20 if v['node_is_primary'] else 10) for k, v in N.items() if v['node_type'] == 'metabolite']
    labels = []
    for rid, r in R.items():
        x, y, anc = L.rlab[rid]; labels.append(('reaction ' + r['bigg_id'], _box(x, y, anc, r['bigg_id'], 22, 0.64), None))
    for k, v in N.items():
        if v['node_type'] == 'metabolite' and '_lab' in v:
            x, y, anc = v['_lab']; s = 26 if v['node_is_primary'] else 20
            labels.append(('metabolite ' + short(v['bigg_id']), _box(x, y, anc, short(v['bigg_id']), s, 0.58), k))
    for tid, tl in L.m['text_labels'].items():
        if '_lab' in tl:
            x, y, anc = tl['_lab']; labels.append(('exit ' + tl['text'], _box(x, y, anc, tl['text'], 22, 0.58), None))
    titles = []
    for z, (x, y) in L.titles.items(): titles.append(('title ' + L.zone_style[z]['title'], _box(x, y, 'start', L.zone_style[z]['title'], 48, 0.62), None))
    for text, z, x, y, _ in getattr(L, 'extra_titles', []): titles.append(('title ' + text, _box(x, y, 'start', text, 40, 0.62), None))
    on_line, hit_node = [], []
    for name, b, own in labels + titles:
        inside = ((allpts[:, 0] > b[0]) & (allpts[:, 0] < b[2]) & (allpts[:, 1] > b[1]) & (allpts[:, 1] < b[3])).sum()
        if inside >= 2: on_line.append(name)
        for k, cx, cy, rr in circ:
            if k == own: continue
            if not (b[2] < cx - rr or cx + rr < b[0] or b[3] < cy - rr or cy + rr < b[1]): hit_node.append(name); break
    allb = labels + titles; lab_over = []
    for i in range(len(allb)):
        for j in range(i + 1, len(allb)):
            a, b = allb[i][1], allb[j][1]
            if not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1]): lab_over.append((allb[i][0], allb[j][0]))
    prim = [(x, y) for k, x, y, r in circ if r == 20]
    close = sum(1 for i in range(len(prim)) for j in range(i + 1, len(prim)) if math.hypot(prim[i][0] - prim[j][0], prim[i][1] - prim[j][1]) < 70)
    xs = [b[0] for _, b, _ in allb] + [b[2] for _, b, _ in allb] + [v['x'] for v in N.values()]
    ys = [b[1] for _, b, _ in allb] + [b[3] for _, b, _ in allb] + [v['y'] for v in N.values()]
    W, H = max(xs) - min(xs), max(ys) - min(ys)
    return dict(line_crossings=ncross, worst_crossing_pairs=[[list(k), v] for k, v in worst],
                labels_on_lines=len(on_line), labels_on_lines_list=on_line[:25],
                labels_on_nodes=len(hit_node), labels_on_nodes_list=hit_node[:25],
                label_label_overlaps=len(lab_over), label_label_overlaps_list=[list(p) for p in lab_over[:25]],
                primary_nodes_closer_than_70=close, width=round(W), height=round(H), bbox_area_Munits2=round(W * H / 1e6, 2))
