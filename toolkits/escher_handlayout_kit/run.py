#!/usr/bin/env python3
"""Command line interface.

  python run.py inspect  --input MAP.json --out DIR
  python run.py start    --input MAP.json --out SPEC.json                # a complete spec of the draft as it is, to be changed
  python run.py patch    --input MAP.json --spec SPEC.json --json '{...}' | --fragment FRAGMENT.json   # change a part of a spec
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
b.add_argument('--max-builds', type=int, default=8, help='changed versions of the spec after which the build says STOP')
k = sub.add_parser('check-spec'); k.add_argument('--input', required=True); k.add_argument('--spec', required=True)
k.add_argument('--partial', action='store_true', help='an unfinished spec: report what is wrong and what is still missing')
s = sub.add_parser('start', help='write a complete spec of the draft map as it is'); s.add_argument('--input', required=True); s.add_argument('--out', required=True)
p = sub.add_parser('patch', help='merge a spec fragment into a spec file (blocks by name, reactions and nodes by id)')
p.add_argument('--input', required=True); p.add_argument('--spec', required=True, help='the spec file to change (created if missing)')
p.add_argument('--json', help='the fragment as JSON text'); p.add_argument('--fragment', help='a file with the fragment')
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
    try: spec = start_spec(args.input)
    except SpecError as e:
        print('SPEC ERROR:', e); sys.exit(2)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True); write_text(args.out, format_spec(spec))
    print(f"start spec written to {args.out}: the draft as it is, {len(spec['layouts'])} blocks ("
          + ', '.join(f"{b['block']}: {len(b['reactions'])}" for b in spec['layouts']) + ' reactions); no merges yet')
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
        gone = {d for m in frag.get('merges', []) for d in m['dups']}      # a merged node has no coordinates of its own any more
        dropped = [d for b in new.get('layouts', []) for d in gone if b.get('nodes', {}).pop(d, None) is not None]
        for b in new.get('layouts', []): [b.get('node_names', {}).pop(d, None) for d in gone]
        P0, _ = compile_spec(old, Layout(args.input), partial=True)
        L = Layout(args.input); P, W = compile_spec(new, L, partial=True)
        fresh = [x for x in P if x not in P0]
        if fresh: raise SpecError('the fragment makes these problems:\n  ' + '\n  '.join(fresh[:60]))
    except SpecError as e:
        print('SPEC ERROR (the spec file is unchanged):', e); sys.exit(2)
    os.makedirs(os.path.dirname(os.path.abspath(args.spec)), exist_ok=True); write_text(args.spec, format_spec(new))
    W += lint(new, L); todo = [r for r in L.R if r not in L.spec_done]
    print(f"patched {args.spec}: {len(L.spec_done)} of {len(L.R)} reactions placed" + (f' ({len(todo)} still to place)' if todo else ' (complete: build it)'))
    if dropped: print(f"  the coordinates of the merged nodes {', '.join(dropped)} were removed (the kept node stands for them)")
    [print('  problem:', x) for x in P]; [print('  warning:', w) for w in W]
else:
    e, w = validate_files(args.input, args.output)
    print('errors:', len(e), 'warnings:', len(w)); [print(' ', x) for x in e[:50]]; sys.exit(1 if e else 0)
