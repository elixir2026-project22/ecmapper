#!/usr/bin/env python3
"""Command line interface.

  python run.py inspect  --input MAP.json --out DIR
  python run.py build    --input MAP.json --spec SPEC.py --out DIR [--curv 0.6] [--no-pack] [--no-png] [--name NAME]
  python run.py check-spec --input MAP.json --spec SPEC.json
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
k = sub.add_parser('check-spec'); k.add_argument('--input', required=True); k.add_argument('--spec', required=True)
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
    try:
        P, W = compile_spec(load_spec_file(args.spec), Layout(args.input))
    except SpecError as e:
        print('SPEC ERROR:', e); sys.exit(2)
    print('spec ok;', len(W), 'warnings'); [print('  warning:', w) for w in W]
else:
    e, w = validate_files(args.input, args.output)
    print('errors:', len(e), 'warnings:', len(w)); [print(' ', x) for x in e[:50]]; sys.exit(1 if e else 0)
