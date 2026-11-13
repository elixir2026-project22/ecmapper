"""End-to-end pipeline: spec -> geometry -> labels -> packing -> JSON + PNG + QA report. No interaction."""
import os, json, importlib.util, traceback
import numpy as np
from .engine import Layout, place_labels
from . import pack as packmod
from . import render as rendermod
from . import qa
from .inspect_map import analyse, to_markdown, base, COFACTORS
from . import directions as dirmod
from . import score as scoremod


from .declspec import SpecError, load_spec_file, compile_spec, validate_schema
from ._io import read_json, write_json, write_text


def load_python_spec(path):
    """LEGACY and UNSAFE: executes the file. Only for specs you wrote yourself (needs allow_python=True)."""
    s = importlib.util.spec_from_file_location('user_spec', path)
    mod = importlib.util.module_from_spec(s); s.loader.exec_module(mod)
    if not hasattr(mod, 'build'): raise SpecError(f'{path} must define build(L)')
    return mod


def apply_spec(L, spec, allow_python=False):
    """spec: dict (declarative JSON spec), path to a .json spec, or (unsafe) path to a .py spec."""
    if isinstance(spec, dict):
        errs = validate_schema(spec)
        if errs: raise SpecError('spec violates the schema:\n  ' + '\n  '.join(errs[:60]))
        compile_spec(spec, L); return
    spec = str(spec)
    if spec.endswith('.json'):
        compile_spec(load_spec_file(spec), L); return
    if not allow_python:
        raise SpecError('Python specs execute arbitrary code and are disabled. Use a JSON spec, or pass --unsafe-python-spec '
                        '(allow_python=True) for a Python spec you wrote yourself.')
    load_python_spec(spec).build(L)


def check_complete(L):
    """The spec must place every reaction and every primary metabolite. Report all gaps at once."""
    miss_r = [f"{rid} ({L.R[rid]['bigg_id']})" for rid in L.R if rid not in getattr(L, 'rxn_geom', {})]
    miss = {}
    for rid, I in L.info.items():
        for _, n in I['subs'] + I['prods']:
            if L.N[n]['node_is_primary'] and n not in L.pos: miss.setdefault(n, set()).add(L.R[rid]['bigg_id'])
    miss_n = [f"{n} ({L.N[n]['bigg_id']}; used by {', '.join(sorted(r))})" for n, r in sorted(miss.items(), key=lambda kv: int(kv[0]))]
    if miss_r or miss_n:
        raise SpecError(f'spec incomplete: {len(miss_r)} reactions and {len(miss_n)} primary nodes are not placed.\n'
                        '  reactions without L.reaction(...): ' + ', '.join(miss_r[:80]) +
                        '\n  primary nodes without coordinates (L.setp): ' + '; '.join(miss_n[:80]))


def content_bbox(L):
    comps = packmod.components(L)
    allr = [r for c in comps for r in c[0]]; alln = set().union(*[c[1] for c in comps])
    tbx = [packmod.tb(x, y, 'start', L.zone_style[z]['title'], 48, 0.62) for z, (x, y) in L.titles.items()]
    p, b = packmod.elements(L, allr, alln, tbx)
    return packmod.bbox_of(p, b)


def detached_pieces(L):
    """The pieces beside the main component that draw a primary metabolite of it again: candidates for a merge
    (or deliberate copies). [dict(reactions=[BiGG ids], shares=[BiGG ids of the metabolites])]"""
    comps = packmod.components(L)
    prim = lambda nodes: {L.N[n]['bigg_id'] for n in nodes if L.N[n]['node_type'] == 'metabolite' and L.N[n]['node_is_primary']}
    main = prim(comps[0][1]) if comps else set()
    out = []
    for rids, nodes in comps[1:]:
        shares = sorted(m for m in prim(nodes) & main if base(m) not in COFACTORS)
        if shares: out.append(dict(reactions=[L.R[r]['bigg_id'] for r in sorted(rids, key=int)], shares=shares))
    return out


def build(input_path, spec_path, out_dir, curv=0.6, do_pack=True, png=True, png_scale=None, name=None, allow_python=False, taxa='all'):
    os.makedirs(out_dir, exist_ok=True)
    name = name or os.path.splitext(os.path.basename(input_path))[0] + '_handlayout'
    L = Layout(input_path)
    apply_spec(L, spec_path, allow_python)
    L.curv = curv
    check_complete(L)
    L.finish(); place_labels(L)
    bb0 = content_bbox(L)
    rep = dict(input=input_path, spec=spec_path if isinstance(spec_path, str) else '<in-memory spec>', curvature=curv, packed=do_pack,
               area_before_packing_Munits2=round((bb0[2] - bb0[0]) * (bb0[3] - bb0[1]) / 1e6, 2))
    if do_pack:
        mbb, cur, info = packmod.pack(L)
        rep['packing'] = [dict(piece=t, dx=dx, dy=dy) for t, dx, dy in info]
        rep['main_component_bbox'] = [round(v) for v in mbb]
    bb = content_bbox(L)
    rep['detached'] = detached_pieces(L)
    rep['placement'] = dirmod.check_layout(L, taxa)       # None without the placement catalogue / without matching subsystems
    L.canvas_bbox = (bb[0] - 150, bb[1] - 150, bb[2] + 150, bb[3] + 150)
    out_json = os.path.join(out_dir, name + '.json')
    rendermod.export_json(L, out_json)
    rep['output_json'] = out_json
    rep['removed_text_labels'] = [f'{t}: {x}' for t, x in getattr(L, 'removed_labels', [])]
    rep['area_after_packing_Munits2'] = round((bb[2] - bb[0]) * (bb[3] - bb[1]) / 1e6, 2)
    rep['canvas'] = dict(width=round(L.canvas_bbox[2] - L.canvas_bbox[0]), height=round(L.canvas_bbox[3] - L.canvas_bbox[1]))
    rep['n_nodes_input'] = len(read_json(input_path)[1]['nodes']); rep['n_nodes_output'] = len(L.N)
    rep['n_reactions'] = len(L.R)
    errs, warns = qa.validate(input_path, out_json)
    rep['validation_errors'] = errs; rep['validation_warnings'] = warns
    rep['metrics'] = qa.metrics(L)
    if png:
        w = L.canvas_bbox[2] - L.canvas_bbox[0]
        sc = png_scale or min(0.4, 3200.0 / w)
        out_png = os.path.join(out_dir, name + '_overview.png')
        rendermod.render(L, out_png, bbox=L.canvas_bbox, scale=sc); rep['output_png'] = out_png
    write_json(os.path.join(out_dir, 'qa_report.json'), rep)
    write_text(os.path.join(out_dir, 'qa_report.md'), report_md(rep))
    return rep


def report_md(rep):
    m = rep['metrics']
    L = ['# Build report', '',
         f"- spec: `{rep['spec']}` | curvature: {rep['curvature']} | packed: {rep['packed']}",
         f"- reactions: {rep['n_reactions']} | nodes: {rep['n_nodes_input']} -> {rep['n_nodes_output']}",
         f"- area before/after packing (M units^2): {rep['area_before_packing_Munits2']} -> {rep['area_after_packing_Munits2']} | canvas {rep['canvas']['width']} x {rep['canvas']['height']}",
         f"- **validation errors: {len(rep['validation_errors'])}**  warnings: {len(rep['validation_warnings'])}", '']
    for e in rep['validation_errors'][:40]: L.append(f'  - ERROR: {e}')
    for w in rep['validation_warnings'][:20]: L.append(f'  - warning: {w}')
    L += ['', '## Layout-quality metrics (lower is better)', '',
          f"- line crossings between different reactions: **{m['line_crossings']}**",
          f"- labels sitting on lines: **{m['labels_on_lines']}**  | labels touching other nodes: **{m['labels_on_nodes']}**  | overlapping label pairs: **{m['label_label_overlaps']}**",
          f"- primary nodes closer than 70 units: {m['primary_nodes_closer_than_70']}",
          f"- bounding box: {m['width']} x {m['height']} ({m['bbox_area_Munits2']} M units^2)", '']
    if m['worst_crossing_pairs']: L += ['Worst crossing pairs: ' + '; '.join(f"{a}x{b} ({n})" for (a, b), n in m['worst_crossing_pairs']), '']
    if m['labels_on_lines_list']: L += ['Labels on lines: ' + ', '.join(m['labels_on_lines_list']), '']
    if m['labels_on_nodes_list']: L += ['Labels on nodes: ' + ', '.join(m['labels_on_nodes_list']), '']
    if m['label_label_overlaps_list']: L += ['Overlapping labels: ' + '; '.join(f'{a} / {b}' for a, b in m['label_label_overlaps_list'][:15]), '']
    if rep.get('detached'):
        L += ['Detached pieces that draw a metabolite of the main map again (merge it, or keep the copy on purpose): '
              + '; '.join(f"{', '.join(d['reactions'][:6])} (shares {', '.join(d['shares'])})" for d in rep['detached']), '']
    L += ['## Score and targets', '', '- ' + scoremod.score_line(rep)] + ['- ' + t.strip() for t in scoremod.target_lines(rep)]
    fix = scoremod.first_fix(rep)
    L += ['', f'Fix first: {fix}' if fix else 'All targets are met.', '']
    pl = rep.get('placement')
    if pl:
        L += ['## Placement of the subsystems (layout_directions catalogue)', '',
              f"- violated directions (hard/strong/medium): **{len(pl['violations'])}** of {pl['checked']} checked | placement penalty: **{pl['penalty']}**"]
        L += [f'  - {v}' for v in pl['violations'][:25]] + ['']
    if rep.get('removed_text_labels'): L += ['## Removed unused title labels of the input map', ''] + [f'- {x}' for x in rep['removed_text_labels']] + ['']
    if rep.get('packing'): L += ['## Packing', ''] + [f"- {p['piece']}: moved by ({p['dx']}, {p['dy']})" for p in rep['packing']]
    return '\n'.join(L) + '\n'


def inspect(input_path, out_dir, taxa='all'):
    os.makedirs(out_dir, exist_ok=True)
    a = analyse(input_path)
    cat = dirmod.Catalogue.load()
    if cat is not None:
        arr = dirmod.Arrangement(dirmod.subsystems_of_analysis(a), cat, taxa)
        a['layout_directions'] = arr.to_dict()
        write_text(os.path.join(out_dir, 'arrangement.md'), arr.text())
    write_json(os.path.join(out_dir, 'analysis.json'), a)
    write_text(os.path.join(out_dir, 'analysis.md'), to_markdown(a))
    return a


def validate_files(input_path, output_path):
    return qa.validate(input_path, output_path)
