"""A complete starting spec: the draft map as it is, one block per subsystem.

A model that cannot write a whole layout from an empty page changes this spec in small steps instead (PROMPT.md):

- one zone and one block per subsystem (the title labels of the draft; a component without title is "Module n"); the
  block's nodes are the primary metabolite nodes at their positions in the draft, relative to the block's origin
  (their top left corner), so moving a subsystem is one change of `origin`;
- every reaction once, in the block of its subsystem: `side` is where the draft draws its cofactors, `shift` its lane
  (parallel reactions between the same nodes), `axis` and `mid` where one side has no primary node;
- `name` / `node_names` say what the ids are (check-spec compares them with the map);
- titles at the draft's title labels; no merges and no splits (decisions of the layout).
"""
import math
from .declspec import SpecError
from .engine import Layout
from .inspect_map import analyse, zone_name

#: lane offsets below this are no lanes (map units)
MIN_SHIFT = 20


def _mean(points):
    points = list(points)
    return (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points)) if points else None


def reaction_entry(L, rid, origin):
    """`side`, `shift` (and `axis`, `mid` if needed) of a reaction as the draft draws it."""
    I = L.info[rid]; N = L.N; pos = L.orig
    entry = {'id': rid, 'name': L.R[rid]['bigg_id']}
    prim = lambda side: [n for _, n in side if N[n]['node_is_primary']]
    a, b = _mean(pos[n] for n in prim(I['subs'])), _mean(pos[n] for n in prim(I['prods']))
    mid = pos[I['mid']]
    ms, mp = pos[I['sm']], pos[I['pm']]
    if a is None or b is None or math.hypot(b[0] - a[0], b[1] - a[1]) < 1:
        # one side without primary node (or both sides at one point): the axis as the draft's markers give it
        dx, dy = mp[0] - ms[0], mp[1] - ms[1]
        if math.hypot(dx, dy) < 1 and (a or b) is not None:
            side = a or b
            dx, dy = (mid[0] - side[0], mid[1] - side[1]) if a is not None else (side[0] - mid[0], side[1] - mid[1])
        if math.hypot(dx, dy) < 1e-9: dx, dy = 1.0, 0.0
        length = math.hypot(dx, dy)
        ux, uy = dx / length, dy / length
        entry['axis'] = [round(ux, 3), round(uy, 3)]
        entry['mid'] = [round(mid[0] - origin[0]), round(mid[1] - origin[1])]
    else:
        length = math.hypot(b[0] - a[0], b[1] - a[1])
        ux, uy = (b[0] - a[0]) / length, (b[1] - a[1]) / length
        # the lane: the midmarker's distance from the line between the two sides
        centre = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        off = (mid[0] - centre[0]) * -uy + (mid[1] - centre[1]) * ux
        if abs(off) >= MIN_SHIFT: entry['shift'] = [round(-uy * off), round(ux * off)]
    # the side of the cofactors, on the screen
    c = _mean(pos[n] for _, n in I['subs'] + I['prods'] if not N[n]['node_is_primary'])
    if c is not None:
        if abs(ux) < abs(uy): entry['side'] = 'R' if c[0] - mid[0] >= 0 else 'L'       # vertical reaction
        else: entry['side'] = 'D' if c[1] - mid[1] >= 0 else 'U'
    return entry


def skeleton_spec(input_path):
    """The spec of the draft without any layout: zones, one block per subsystem with its reactions (`id`, `name`) and
    `node_names` (which metabolite each primary node id is), but no coordinates, no sides, no merges and no titles.
    For someone who designs the layout from scratch: every position is theirs to decide."""
    spec = start_spec(input_path)
    spec['notes'] = ('Skeleton written by `run.py start --skeleton`: what has to be placed. Give every id of node_names a position in '
                     'nodes (or on a ring), every reaction its keys, every zone a title position; add merges and splits.')
    for b in spec['layouts']:
        b['origin'] = [0, 0]; b['nodes'] = {}
        b['reactions'] = [{'id': r['id'], 'name': r['name']} for r in b['reactions']]
    spec['titles'] = {}
    return spec


def start_spec(input_path):
    """The spec of the draft map `input_path` as it is drawn (a dict; write it with declspec.format_spec)."""
    a = analyse(input_path)
    if a['unsupported_reactions']:
        raise SpecError('unsupported reaction topology (no multimarkers): ' + ', '.join(r['bigg_id'] for r in a['unsupported_reactions'][:20]))
    L = Layout(input_path); N = L.N
    labels = {t['text']: t for t in a['title_labels']}
    by_comp = {}
    for r in a['reactions']: by_comp.setdefault(r['component'], []).append(r['id'])
    zones, order = {}, []
    for c in a['components']:
        title = c['suggested_subsystem'] or f"Module {c['index'] + 1}"
        z = zone_name(title)
        if z not in zones: zones[z] = dict(title=title, rids=[]); order.append(z)
        zones[z]['rids'] += by_comp.get(c['index'], [])
    spec = {'spec_version': 1,
            'notes': 'Starting spec written by `run.py start`: the draft map as it is (no merges yet). Change it in small steps.',
            'defaults': {'bend': 60}, 'zones': [], 'splits': [], 'merges': [], 'layouts': [], 'titles': {}}
    for z in order:
        title, rids = zones[z]['title'], sorted(zones[z]['rids'], key=int)
        spec['zones'].append(dict(name=z, title=z, label=title if title in labels else None))
        members = []
        for rid in rids:
            for _, n in L.info[rid]['subs'] + L.info[rid]['prods']:
                if N[n]['node_is_primary'] and n not in members: members.append(n)
        if not members: continue
        ox = math.floor(min(L.orig[n][0] for n in members) / 10) * 10
        oy = math.floor(min(L.orig[n][1] for n in members) / 10) * 10
        spec['layouts'].append(dict(block=z, zone=z, origin=[ox, oy],
                                    nodes={n: [round(L.orig[n][0] - ox), round(L.orig[n][1] - oy)] for n in members},
                                    node_names={n: N[n]['bigg_id'] for n in members},
                                    reactions=[reaction_entry(L, rid, (ox, oy)) for rid in rids]))
        spec['titles'][z] = [round(labels[title]['x']), round(labels[title]['y'])] if title in labels else [ox, oy - 80]
    return spec
