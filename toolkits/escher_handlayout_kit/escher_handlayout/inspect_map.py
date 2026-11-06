"""Phase-0 analysis of a draft Escher map: everything an automated run needs to plan the layout."""
import json, math
from collections import defaultdict, Counter
import numpy as np
from ._io import read_json

HUB_BASE = {'glc__D', 'g6p', 'f6p', 'fdp', 'g3p', 'dhap', 'pep', 'pyr', 'accoa', 'acald', 'oaa', 'akg', 'succ', 'mal__L',
            'cit', 'icit', 'r5p', 'ru5p__D', 'xu5p__D', 's7p', 'e4p', '6pgc', 'ac', 'actp', 'lac__D', 'lac__L', 'etoh', 'fum'}

def base(b):
    for suf in ('_c', '_e', '_p'):
        if b.endswith(suf): return b[:-2]
    return b

def analyse(path):
    d = read_json(path); meta, m = d
    N, R, TL = m['nodes'], m['reactions'], m['text_labels']
    par = {}
    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    rx = {}
    for rid, r in R.items():
        mid = None; subs = []; prods = []
        nodes = set()
        for s in r['segments'].values():
            a, b = s['from_node_id'], s['to_node_id']; nodes |= {a, b}
            ta, tb = N[a]['node_type'], N[b]['node_type']
            if ta == 'midmarker': mid = a
            if tb == 'midmarker': mid = b
            if ta == 'metabolite' and tb == 'multimarker': subs.append(a)
            if ta == 'multimarker' and tb == 'metabolite': prods.append(b)
        has_mult = any(N[n]['node_type'] == 'multimarker' for n in nodes)
        if mid is None: mid = next(iter(nodes))          # fallback, only used for the nearest-title vote
        rx[rid] = dict(mid=mid, subs=subs, prods=prods, nodes=nodes, supported=has_mult and any(N[n]['node_type'] == 'midmarker' for n in nodes))
        for n in nodes:
            if N[n]['node_type'] != 'metabolite' or N[n]['node_is_primary']: par[f('n' + n)] = f('r' + rid)
    # text labels: titles vs exit labels
    exits = {t: v for t, v in TL.items() if v['text'].startswith('→')}
    titles = {t: v for t, v in TL.items() if t not in exits}
    comps = defaultdict(list)
    for rid in R: comps[f('r' + rid)].append(rid)
    comp_list = sorted(comps.values(), key=lambda c: -len(c))
    def nearest_title(rids):
        votes = Counter()
        for rid in rids:
            mx, my = N[rx[rid]['mid']]['x'], N[rx[rid]['mid']]['y']
            if titles:
                votes[min(titles, key=lambda t: math.hypot(titles[t]['x'] - mx, titles[t]['y'] - my))] += 1
        return votes.most_common(1)[0][0] if votes else None
    out = dict(file=path, map_name=meta.get('map_name'), n_reactions=len(R), n_nodes=len(N), n_text_labels=len(TL), canvas=m['canvas'],
               title_labels=[dict(id=t, text=v['text'], x=round(v['x']), y=round(v['y'])) for t, v in titles.items()],
               exit_labels=[dict(id=t, text=v['text'], x=round(v['x']), y=round(v['y'])) for t, v in exits.items()])
    # components with suggested subsystem
    out['components'] = []
    rid_comp = {}
    for ci, rids in enumerate(comp_list):
        xs = [N[n]['x'] for rid in rids for n in rx[rid]['nodes']]; ys = [N[n]['y'] for rid in rids for n in rx[rid]['nodes']]
        t = nearest_title(rids)
        out['components'].append(dict(index=ci, n_reactions=len(rids), suggested_subsystem=titles[t]['text'] if t else None,
                                      reactions=[R[r]['bigg_id'] for r in sorted(rids, key=int)],
                                      bbox=[round(min(xs)), round(min(ys)), round(max(xs)), round(max(ys))]))
        for r in rids: rid_comp[r] = ci
    # reactions
    def lab(n): return N[n]['bigg_id'] + ('*' if N[n]['node_is_primary'] else '')
    out['reactions'] = []
    for rid in sorted(R, key=int):
        r = R[rid]; x = rx[rid]
        out['reactions'].append(dict(id=rid, bigg_id=r['bigg_id'], name=r.get('name'), reversible=r['reversibility'],
                                     component=rid_comp[rid], mid=[round(N[x['mid']]['x']), round(N[x['mid']]['y'])],
                                     substrates=[lab(n) for n in x['subs']], products=[lab(n) for n in x['prods']],
                                     primary_nodes=[n for n in x['subs'] + x['prods'] if N[n]['node_is_primary']]))
    # traps
    use = defaultdict(set)
    for rid, x in rx.items():
        for n in x['nodes']:
            if N[n]['node_type'] == 'metabolite' and not N[n]['node_is_primary']: use[n].add(rid)
    out['unsupported_reactions'] = [dict(id=rid, bigg_id=R[rid]['bigg_id']) for rid, x in rx.items() if not x['supported']]
    out['shared_secondary_nodes'] = [dict(node=n, bigg_id=N[n]['bigg_id'], reactions=sorted(rs, key=int)) for n, rs in use.items() if len(rs) > 1]
    bybig = defaultdict(list)
    for rid, r in R.items(): bybig[r['bigg_id']].append(rid)
    out['duplicate_reaction_bigg_ids'] = {b: sorted(v, key=int) for b, v in bybig.items() if len(v) > 1}
    fam = defaultdict(list)
    for rid, x in rx.items():
        key = (tuple(sorted(n for n in x['subs'] if N[n]['node_is_primary'])), tuple(sorted(n for n in x['prods'] if N[n]['node_is_primary'])))
        if key[0] and key[1]: fam[key].append(rid)
    out['parallel_reaction_families'] = [[R[r]['bigg_id'] for r in sorted(v, key=int)] for v in fam.values() if len(v) > 1]
    # merge candidates: same metabolite drawn as several primary nodes
    prim = defaultdict(list)
    for n, v in N.items():
        if v['node_type'] == 'metabolite' and v['node_is_primary']: prim[v['bigg_id']].append(n)
    cand = []
    for b, ns in prim.items():
        if len(ns) < 2: continue
        comps_of = sorted({rid_comp[r] for rid, x in rx.items() for r in [rid] if any(n in x['nodes'] for n in ns)})
        cand.append(dict(bigg_id=b, nodes=ns, components=comps_of, in_default_hub_list=base(b) in HUB_BASE,
                         reactions=sorted({R[rid]['bigg_id'] for rid, x in rx.items() if any(n in x['nodes'] for n in ns)})))
    out['merge_candidates'] = sorted(cand, key=lambda c: (-int(c['in_default_hub_list']), -len(c['nodes'])))
    # exit labels -> nearest primary node
    pr = [(n, np.array([v['x'], v['y']])) for n, v in N.items() if v['node_type'] == 'metabolite' and v['node_is_primary']]
    out['exit_label_assignment'] = [dict(label=v['text'], node=min(pr, key=lambda kv: np.linalg.norm(kv[1] - np.array([v['x'], v['y']])))[0]) for v in exits.values()]
    return out

def to_markdown(a):
    L = [f"# Analysis of `{a['file']}`", '',
         f"- map: {a['map_name']} | reactions: {a['n_reactions']} | nodes: {a['n_nodes']} | text labels: {a['n_text_labels']} | canvas: {a['canvas']}", '',
         '## Title labels (candidate subsystem names)', '']
    L += [f"- `{t['id']}`: **{t['text']}** at ({t['x']}, {t['y']})" for t in a['title_labels']] or ['- none']
    L += ['', '## Connected components (over primary nodes)', '', '| # | reactions | suggested subsystem | bbox | reaction ids |', '|---|---|---|---|---|']
    for c in a['components']:
        L.append(f"| {c['index']} | {c['n_reactions']} | {c['suggested_subsystem']} | {c['bbox']} | {', '.join(c['reactions'][:14])}{' ...' if len(c['reactions']) > 14 else ''} |")
    L += ['', '## Reactions (`*` = primary metabolite node)', '', '| id | bigg_id | rev | comp | mid | substrates -> products |', '|---|---|---|---|---|---|']
    for r in a['reactions']:
        L.append(f"| {r['id']} | {r['bigg_id']} | {'r' if r['reversible'] else '-'} | {r['component']} | {r['mid']} | {' + '.join(r['substrates'])} -> {' + '.join(r['products'])} |")
    L += ['', '## Traps / hints', '',
          f"- reactions with unsupported topology (no multimarkers; cannot be laid out by this engine): {[r['bigg_id'] for r in a['unsupported_reactions']] or 'none'}",
          f"- shared secondary nodes (the engine duplicates them per reaction automatically): {len(a['shared_secondary_nodes'])}",
          f"- duplicate reaction bigg_ids (draw as parallel lanes): {a['duplicate_reaction_bigg_ids']}",
          f"- parallel reaction families (same primary substrates/products; stack on one lane or draw as lanes): {a['parallel_reaction_families']}", '',
          '## Merge candidates (same metabolite drawn as several primary nodes)', '', '| metabolite | node ids | components | default hub | reactions |', '|---|---|---|---|---|']
    for c in a['merge_candidates']:
        L.append(f"| {c['bigg_id']} | {c['nodes']} | {c['components']} | {'yes' if c['in_default_hub_list'] else ''} | {', '.join(c['reactions'][:10])} |")
    L += ['', '## Exit labels -> nearest primary node', '']
    L += [f"- `{e['label']}` -> node {e['node']}" for e in a['exit_label_assignment']]
    return '\n'.join(L) + '\n'
