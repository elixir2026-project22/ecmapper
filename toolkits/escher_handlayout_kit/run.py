#!/usr/bin/env python3
"""Command line interface.

  python run.py inspect  --input MAP.json --out DIR
  python run.py start    --input MAP.json --out SPEC.json --skeleton     # what has to be placed, without any layout (PROMPT_AGENT.md)
  python run.py start    --input MAP.json --out SPEC.json [--as-is]      # a complete spec to be changed: the subsystems of the draft,
                                                                         # arranged and connected (--as-is: the draft as it is)
  python run.py patch    --input MAP.json --spec SPEC.json --json '{...}' | --fragment FRAGMENT.json   # change a part of a spec
  python run.py tidy     --input MAP.json --spec SPEC.json [--block NAME ...] [--only align|untangle|cofactors] [--copies]
                                                                         # rows and columns, loose ends and lanes untangled, cofactors in
                                                                         # free space: only what makes the metrics better
  python run.py view     --input MAP.json --spec SPEC.json --out DIR [--around pyr_c PYK ...] [--radius 900] [--no-marks]
                                                                         # pictures to look at: the overview, or the surroundings of
                                                                         # metabolites / reactions, with the defects marked
  python run.py script   FILE.py [ARGS]                                  # run your own script with the kit importable (escher_handlayout.shapes)
  python run.py build    --input MAP.json --spec SPEC.json --out DIR [--curv 0.6] [--no-pack] [--no-png] [--name NAME]
                                                                         # prints score, targets and what to do next
  python run.py check-spec --input MAP.json --spec SPEC.json [--partial]
  python run.py validate --input MAP.json --output OUT.json
  python run.py directions --input MAP.json [--taxa bacteria] [--out arrangement.md]   # usual arrangement of the subsystems
  python run.py directions --input MAP.json --check SPEC.json                           # the directions a spec violates
"""
import argparse, sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escher_handlayout import build, inspect, validate_files
from escher_handlayout.pipeline import SpecError

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = ap.add_subparsers(dest='cmd', required=True)
a = sub.add_parser('inspect'); a.add_argument('--input', required=True); a.add_argument('--out', required=True); a.add_argument('--taxa', default='all')
b = sub.add_parser('build'); b.add_argument('--input', required=True); b.add_argument('--spec', required=True); b.add_argument('--out', required=True)
b.add_argument('--curv', type=float, default=0.6); b.add_argument('--no-pack', action='store_true'); b.add_argument('--no-png', action='store_true')
b.add_argument('--name'); b.add_argument('--png-scale', type=float); b.add_argument('--taxa', default='all', help='organism group for the placement check (all, bacteria, plants, animals, fungi, ...)')
b.add_argument('--unsafe-python-spec', action='store_true', help='allow a .py spec (executes code!)')
b.add_argument('--history', help='history file of the builds (default: builds.json beside the --out folder)')
b.add_argument('--no-history', action='store_true', help='neither record this build nor say what to do next')
b.add_argument('--max-builds', type=int, default=16, help='changed versions of the spec after which the build says STOP')
k = sub.add_parser('check-spec'); k.add_argument('--input', required=True); k.add_argument('--spec', required=True)
k.add_argument('--partial', action='store_true', help='an unfinished spec: report what is wrong and what is still missing')
s = sub.add_parser('start', help='write a complete spec of the draft map: its subsystems arranged and connected'); s.add_argument('--input', required=True); s.add_argument('--out', required=True)
s.add_argument('--as-is', action='store_true', help='the draft as it is: every subsystem a separate block, no merges')
s.add_argument('--skeleton', action='store_true', help='no layout at all: the blocks with their reactions and node names, for a layout written from scratch')
s.add_argument('--taxa', default='all', help='organism group for the arrangement (all, bacteria, plants, animals, fungi, ...)')
p = sub.add_parser('patch', help='merge a spec fragment into a spec file (blocks by name, reactions and nodes by id)')
p.add_argument('--input', required=True); p.add_argument('--spec', required=True, help='the spec file to change (created if missing)')
p.add_argument('--json', help='the fragment as JSON text'); p.add_argument('--fragment', help='a file with the fragment')
t = sub.add_parser('tidy', help='tidy the drawing of a complete spec (rows and columns, side and distance of the cofactors); the design stays')
t.add_argument('--input', required=True); t.add_argument('--spec', required=True, help='the spec file to tidy (changed in place)')
t.add_argument('--block', action='append', help='only this block (may be given several times)'); t.add_argument('--only', choices=['align', 'untangle', 'cofactors'])
t.add_argument('--copies', action='store_true', help='also give a reaction with a long line to a hub metabolite its own copy of it, where that saves crossings')
t.add_argument('--seconds', type=float, default=90.0, help='stop trying after this time'); t.add_argument('--taxa', default='all')
v = sub.add_parser('view', help='render the layout of a spec to look at: the overview or the surroundings of names, defects marked')
v.add_argument('--input', required=True); v.add_argument('--spec', required=True); v.add_argument('--out', required=True, help='folder for the pictures')
v.add_argument('--around', nargs='*', default=[], help='BiGG ids of metabolites (pyr_c) or reactions (PYK), or node ids: one picture each')
v.add_argument('--radius', type=float, default=900.0, help='half the width of a picture around a name, in map units')
v.add_argument('--no-marks', action='store_true'); v.add_argument('--curv', type=float, default=0.6); v.add_argument('--name', default='view')
x = sub.add_parser('script', help='run a Python script of yours with the kit on the path (for a script that writes the spec)')
x.add_argument('file'); x.add_argument('args', nargs=argparse.REMAINDER)
c = sub.add_parser('validate'); c.add_argument('--input', required=True); c.add_argument('--output', required=True)
d = sub.add_parser('directions', help='where the subsystems usually sit (layout_directions catalogue), or --check SPEC')
d.add_argument('--input', required=True); d.add_argument('--taxa', default='all', help='all, bacteria, plants, animals, fungi, archaea, ...')
d.add_argument('--out', help='write the brief to this file instead of printing it')
d.add_argument('--check', metavar='SPEC', help='list the directions a layout spec violates, and their penalty')
args = ap.parse_args()
if args.cmd == 'inspect':
    r = inspect(args.input, args.out, taxa=args.taxa)
    print(f"analysis written to {args.out}/analysis.md and analysis.json  ({r['n_reactions']} reactions, {len(r['components'])} components)"
          + (f"; arrangement of the subsystems in {args.out}/arrangement.md" if 'layout_directions' in r else ''))
elif args.cmd == 'build':
    try:
        rep = build(args.input, args.spec, args.out, curv=args.curv, do_pack=not args.no_pack, png=not args.no_png, png_scale=args.png_scale, name=args.name, allow_python=args.unsafe_python_spec, taxa=args.taxa)
    except SpecError as e:
        print('SPEC ERROR:', e); sys.exit(2)
    m = rep['metrics']
    print(f"output: {rep['output_json']}")
    print(f"validation errors: {len(rep['validation_errors'])}  warnings: {len(rep['validation_warnings'])}")
    print(f"crossings={m['line_crossings']} labels_on_lines={m['labels_on_lines']} labels_on_nodes={m['labels_on_nodes']} "
          f"label_overlaps={m['label_label_overlaps']} area={rep['area_after_packing_Munits2']}M"
          + (f" placement_penalty={rep['placement']['penalty']} ({len(rep['placement']['violations'])} violated directions)" if rep.get('placement') else ''))
    from escher_handlayout import score as S
    print(S.score_line(rep)); print('\n'.join(S.target_lines(rep)))
    if not args.no_history:
        hist = S.record(args.history or os.path.join(os.path.dirname(os.path.abspath(args.out)), 'builds.json'), rep, args.out, args.spec)
        print('builds: ' + ' | '.join(f"{h['out']} {h['score']}" for h in hist[-10:]))
        print(S.verdict(hist, rep, args.max_builds)[1])
    sys.exit(1 if rep['validation_errors'] else 0)
elif args.cmd == 'directions':
    from escher_handlayout import directions as D
    from escher_handlayout.engine import Layout
    from escher_handlayout.declspec import load_spec_file, compile_spec
    from escher_handlayout.inspect_map import analyse
    cat = D.Catalogue.load()
    if cat is None: sys.exit(f'placement catalogue not found: {D.CATALOGUE}')
    L = Layout(args.input)
    if args.check:
        spec = load_spec_file(args.check); compile_spec(spec, L, partial=True)
        arr = D.Arrangement(D.subsystems_of_layout(L), cat, args.taxa)
        res = D.check(arr, *D.spec_geometry(spec, L))
        for v in res['violations']: print('violated:', v)
        print(f"placement penalty: {res['penalty']} ({len(res['violations'])} of {res['checked']} directions violated)")
    else:
        arr = D.Arrangement(D.subsystems_of_analysis(analyse(args.input)), cat, args.taxa)
        if args.out:
            from escher_handlayout._io import write_text
            write_text(args.out, arr.text())
            print(f'arrangement written to {args.out} ({len(arr.subsystem_directions)} directions between subsystems, {len(arr.metabolite_directions)} between metabolites)')
        else: sys.stdout.write(arr.text())
elif args.cmd == 'check-spec':
    from escher_handlayout.engine import Layout
    from escher_handlayout.declspec import load_spec_file, compile_spec
    from escher_handlayout.lint import lint
    try:
        spec = load_spec_file(args.spec); L = Layout(args.input)
        P, W = compile_spec(spec, L, partial=args.partial)
    except SpecError as e:
        print('SPEC ERROR:', e); sys.exit(2)
    W += lint(spec, L)
    if args.partial:
        todo = [f"{r} ({L.R[r]['bigg_id']})" for r in L.R if r not in L.spec_done]
        if P: print('SPEC ERROR: spec problems:\n  ' + '\n  '.join(P[:100]))
        print(f'{len(L.spec_done)} of {len(L.R)} reactions placed' + ('; still to place: ' + ', '.join(todo[:60]) + (' ...' if len(todo) > 60 else '') if todo else ''))
    print('spec ok;' if not P else 'spec not ok;', len(W), 'warnings'); [print('  warning:', w) for w in W]
    if P: sys.exit(2)
elif args.cmd == 'start':
    from escher_handlayout.start import start_spec
    from escher_handlayout.declspec import format_spec
    from escher_handlayout._io import write_text
    try:
        if args.skeleton:
            from escher_handlayout.start import skeleton_spec
            spec = skeleton_spec(args.input); os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True); write_text(args.out, format_spec(spec))
            print(f"skeleton written to {args.out}: {len(spec['layouts'])} blocks, {sum(len(b['reactions']) for b in spec['layouts'])} reactions and "
                  f"{sum(len(b['node_names']) for b in spec['layouts'])} primary metabolite nodes to place (no coordinates, no merges, no titles)")
            sys.exit(0)
        if args.as_is: spec, notes = start_spec(args.input), None
        else:
            from escher_handlayout.connect import connected_start
            spec, notes = connected_start(args.input, args.taxa)
    except SpecError as e:
        print('SPEC ERROR:', e); sys.exit(2)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True); write_text(args.out, format_spec(spec))
    print(f"start spec written to {args.out}: {'the draft as it is' if notes is None else 'the subsystems of the draft, arranged and connected'}, {len(spec['layouts'])} blocks ("
          + ', '.join(f"{b['block']}: {len(b['reactions'])}" for b in spec['layouts']) + ' reactions)' + ('; no merges yet' if notes is None else f"; {len(spec['merges'])} merges"))
    for i, n in enumerate(notes or []): print(('' if i == 0 else '  ') + n)
elif args.cmd == 'patch':
    from escher_handlayout.engine import Layout
    from escher_handlayout.declspec import load_spec_file, compile_spec, parse_json_strict, validate_schema, merge_fragment, format_spec
    from escher_handlayout.lint import lint
    from escher_handlayout._io import read_text, write_text
    try:
        if (args.json is None) == (args.fragment is None): raise SpecError('give the fragment either as --json TEXT or as --fragment FILE')
        frag = parse_json_strict(args.json if args.json is not None else read_text(args.fragment))
        errs = validate_schema(frag)
        if errs: raise SpecError('the fragment violates the schema:\n  ' + '\n  '.join(errs[:40]))
        old = load_spec_file(args.spec) if os.path.exists(args.spec) else {'spec_version': 1}
        new = merge_fragment(old, frag)
        gone = {d for m in frag.get('merges', []) for d in m.get('dups', []) if not m.get('remove')}      # a merged node has no coordinates of its own any more
        dropped = [d for b in new.get('layouts', []) for d in gone if b.get('nodes', {}).pop(d, None) is not None]
        for b in new.get('layouts', []): [b.get('node_names', {}).pop(d, None) for d in gone]
        P0, _ = compile_spec(old, Layout(args.input), partial=True)
        # a node that is no dup any more is a copy again: it gets coordinates beside the node that stood for it
        dups = lambda spec: {d: m['keep'] for m in spec.get('merges', []) for d in m.get('dups', [])}
        back = {d: k for d, k in dups(old).items() if d not in dups(new) and d != k
                and not any(d in b.get('nodes', {}) or any(d in r['nodes'] for r in b.get('rings', [])) for b in new.get('layouts', []))}
        restored = []
        if back:
            L1 = Layout(args.input); compile_spec(new, L1, partial=True); raw = Layout(args.input)
            for i, (d, k) in enumerate(sorted(back.items(), key=lambda kv: int(kv[0]) if kv[0].isdigit() else 0)):
                users = {rid for rid, I in raw.info.items() if d in [n for _, n in I['subs'] + I['prods']]}
                b = next((b for b in new.get('layouts', []) if users & {r['id'] for r in b.get('reactions', [])}), None)
                if b is None or d not in raw.N: continue
                kp = L1.pos.get(k) or L1.pos.get(dups(new).get(k)) or (raw.N[d]['x'], raw.N[d]['y'])
                o = b.get('origin', [0, 0]); xy = [round(kp[0] + 160 - o[0]), round(kp[1] + 160 * (i + 1) - o[1])]
                b.setdefault('nodes', {})[d] = xy; b.setdefault('node_names', {})[d] = raw.N[d]['bigg_id']
                restored.append(f"{d} ({raw.N[d]['bigg_id']}) at {xy} in block '{b['block']}'")
        L = Layout(args.input); P, W = compile_spec(new, L, partial=True)
        fresh = [x for x in P if x not in P0]
        if fresh: raise SpecError('the fragment makes these problems:\n  ' + '\n  '.join(fresh[:60]))
    except SpecError as e:
        print('SPEC ERROR (the spec file is unchanged):', e); sys.exit(2)
    os.makedirs(os.path.dirname(os.path.abspath(args.spec)), exist_ok=True); write_text(args.spec, format_spec(new))
    W += lint(new, L); todo = [r for r in L.R if r not in L.spec_done]
    print(f"patched {args.spec}: {len(L.spec_done)} of {len(L.R)} reactions placed" + (f' ({len(todo)} still to place)' if todo else ' (complete: build it)'))
    if dropped: print(f"  the coordinates of the merged nodes {', '.join(dropped)} were removed (the kept node stands for them)")
    if restored: print('  no longer merged, so copies again, placed beside the node that stood for them (move them where they belong): ' + '; '.join(restored))
    [print('  problem:', x) for x in P]; [print('  warning:', w) for w in W]
elif args.cmd == 'view':
    from escher_handlayout import view as V
    try: shown = V.view(args.input, args.spec, args.out, args.around, args.radius, not args.no_marks, args.curv, args.name)
    except (SpecError, ValueError) as e:
        print('SPEC ERROR:', e); sys.exit(2)
    for path, what, found in shown:
        print(f'{path}: {what}')
        if not args.no_marks: [print('  ' + line) for line in V.describe(found)]
elif args.cmd == 'script':
    import runpy
    sys.argv = [args.file] + args.args
    runpy.run_path(args.file, run_name='__main__')
elif args.cmd == 'tidy':
    from escher_handlayout.tidy import tidy
    from escher_handlayout.declspec import load_spec_file, format_spec
    from escher_handlayout._io import write_text
    try:
        old = load_spec_file(args.spec)
        new, notes = tidy(args.input, old, args.block, do_align=args.only in (None, 'align'), do_untangle=args.only in (None, 'untangle'), do_cofactors=args.only in (None, 'cofactors'),
                           do_copies=args.copies, seconds=args.seconds, taxa=args.taxa)
    except (SpecError, ValueError) as e:
        print('SPEC ERROR (the spec file is unchanged):', e); sys.exit(2)
    if new != old: write_text(args.spec, format_spec(new))
    print(f'tidied {args.spec} (build it):' if new != old else f'{args.spec} is unchanged:')
    for n in notes: print('  ' + n)
else:
    e, w = validate_files(args.input, args.output)
    print('errors:', len(e), 'warnings:', len(w)); [print(' ', x) for x in e[:50]]; sys.exit(1 if e else 0)
