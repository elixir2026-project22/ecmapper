#!/usr/bin/env python3
"""Command line interface of escher_rules_kit: lay out and judge Escher maps by a rule set (rules/*.yaml).

The rule set
  python run.py rules   [--type motif] [--strength must should] [--group Duplication] [--id M01 E10] [--bindings]
  python run.py rules   --validate                    # check the rule-set and checks files, say what is wrong
  python run.py render  [--out RULES.md]              # the rule set as Markdown tables
  python run.py probes                                # what the kit can find out about a map, with the parameters

A map against the rule set
  python run.py brief   --input DRAFT.json --out DIR            # before the layout: which rules apply where
  python run.py check   --input MAP.json [--draft DRAFT.json] [--previous OLD.json] [--out DIR] [--strict]
                                                                # which rules the map breaks

The layout itself is done by escher_handlayout_kit (the engine); its commands are run through this kit:
  python run.py inspect | start | patch | tidy | view | script | check-spec | build | validate | directions  ...
      inspect also writes the brief (rules_brief.md) into its --out folder,
      build also checks the built map and ends with the verdict of the rules (rules_report.md in its --out folder).
  python run.py doc     [NAME] [--part N]             # read a document of the engine (format_reference, guide, ...)

Every command takes --rules FILE [FILE ...] (default: all rule-set files of rules/); a file NAME.yaml is read
together with NAME.checks.yaml if that exists.
"""
import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from escher_rules import probes as P, report, ruleset  # noqa: E402

ENGINE_COMMANDS = ('inspect', 'start', 'patch', 'tidy', 'view', 'script', 'check-spec', 'build', 'validate', 'directions')
ENGINE_ENV = 'ESCHER_HANDLAYOUT_KIT'
DOCS = {'format_reference': 'escher_handlayout/agent/prompts/format_reference.md', 'guide': 'docs/GUIDE.md',
        'recipes': 'escher_handlayout/agent/prompts/recipes.md', 'prompt': 'PROMPT.md', 'prompt_agent': 'PROMPT_AGENT.md',
        'readme': 'README.md', 'directions': 'layout_directions/README_v2.md'}
PART = 15000  # characters of a document per part


def engine_dir():
    """The folder of escher_handlayout_kit: $ESCHER_HANDLAYOUT_KIT, a folder `engine` in this kit, or the kit's
    sibling folder. None if there is none."""
    for folder in (os.environ.get(ENGINE_ENV), os.path.join(HERE, 'engine'), os.path.join(os.path.dirname(HERE), 'escher_handlayout_kit')):
        if folder and os.path.isfile(os.path.join(folder, 'run.py')) and os.path.isdir(os.path.join(folder, 'escher_handlayout')):
            return os.path.abspath(folder)
    return None


def need_engine():
    folder = engine_dir()
    if folder is None:
        sys.exit('The layout commands need escher_handlayout_kit: put it beside this kit (toolkits/escher_handlayout_kit), '
                 f'into a folder `engine` in this kit, or set {ENGINE_ENV} to its folder. The rule commands (rules, render, '
                 'probes, brief, check) work without it.')
    return folder


def load_rules(paths):
    try:
        return ruleset.load(paths, probes=P.PROBES)
    except ruleset.RuleError as e:
        sys.exit(f'RULE ERROR: {e}')


def take(argv, name, many=False, flag=False):
    """Remove option `name` from argv and return its value(s) (None if absent)."""
    if name not in argv:
        return None
    i = argv.index(name)
    j = i + 1
    if not flag:
        while j < len(argv) and not argv[j].startswith('--') and (many or j == i + 1):
            j += 1
    values = argv[i + 1:j]
    del argv[i:j]
    return True if flag else (values if many else (values[0] if values else None))


def value(argv, name):
    return argv[argv.index(name) + 1] if name in argv and argv.index(name) + 1 < len(argv) else None


def forward(cmd, argv):
    """Run a command of the engine; after `inspect` add the brief, after `build` the check of the built map."""
    engine = need_engine()
    rule_files = take(argv, '--rules', many=True)
    no_rules = take(argv, '--no-rules', flag=True)
    previous = take(argv, '--previous')
    env = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1', MPLBACKEND=os.environ.get('MPLBACKEND', 'Agg'))
    env['PYTHONPATH'] = os.pathsep.join(p for p in (engine, HERE, env.get('PYTHONPATH', '')) if p)
    done = subprocess.run([sys.executable, os.path.join(engine, 'run.py'), cmd, *argv], env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
    sys.stdout.write(done.stdout)
    if done.returncode not in (0, 1) or no_rules or cmd not in ('inspect', 'build'):
        sys.exit(done.returncode)
    draft, out = value(argv, '--input'), value(argv, '--out')
    rs = load_rules(rule_files)
    try:
        if cmd == 'inspect' and done.returncode == 0:
            report.brief(rs, draft, out)
            print(f'rules that apply to this map, and where: {out}/rules_brief.md (read it before you design)')
        elif cmd == 'build':
            built = re.search(r'^output: (.+)$', done.stdout, re.M)
            if built:
                _, _, summary = report.check(rs, built.group(1).strip(), draft=draft, previous=previous, out_dir=out)
                print('\n'.join(summary))
    except (OSError, ValueError, KeyError) as e:
        print(f'rules: not checked ({type(e).__name__}: {e})')
    sys.exit(done.returncode)


def show_doc(name, part):
    engine = need_engine()
    if not name:
        print('documents of the engine (python run.py doc NAME [--part N]):')
        for key, rel in DOCS.items():
            path = os.path.join(engine, rel)
            if os.path.isfile(path):
                print(f'  {key:17s}{rel}  ({os.path.getsize(path) // 1000 + 1} k characters)')
        return
    if name not in DOCS or not os.path.isfile(os.path.join(engine, DOCS[name])):
        sys.exit(f"no document `{name}`; there are: {', '.join(DOCS)}")
    with open(os.path.join(engine, DOCS[name]), encoding='utf-8') as f:
        text = f.read()
    parts, lines, size = [], [], 0
    for line in text.splitlines(keepends=True):  # parts end at line ends
        if size + len(line) > PART and lines:
            parts.append(''.join(lines))
            lines, size = [], 0
        lines.append(line)
        size += len(line)
    parts.append(''.join(lines))
    part = min(max(part, 1), len(parts))
    sys.stdout.write(parts[part - 1])
    if len(parts) > 1:
        print(f'\n[{name}: part {part} of {len(parts)}' + (f'; next: doc {name} --part {part + 1}]' if part < len(parts) else ']'))


def main(argv):
    if argv and argv[0] in ENGINE_COMMANDS:
        forward(argv[0], list(argv[1:]))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--rules', nargs='+', metavar='FILE', help='rule-set files (default: all of rules/)')
    select = argparse.ArgumentParser(add_help=False)
    select.add_argument('--id', nargs='+', help='only these rules')
    select.add_argument('--skip', nargs='+', help='not these rules')
    select.add_argument('--type', nargs='+', help='only rules of these types')
    select.add_argument('--strength', nargs='+', help='only rules of these strengths')
    select.add_argument('--group', nargs='+', help='only rules whose group contains one of these texts')
    a = sub.add_parser('rules', parents=[common, select], help='list the rules')
    a.add_argument('--validate', action='store_true', help='only check the files')
    a.add_argument('--bindings', action='store_true', help='show the probe of every rule')
    a = sub.add_parser('render', parents=[common, select], help='the rule set as Markdown')
    a.add_argument('--out', help='write to this file instead of printing')
    a.add_argument('--bindings', action='store_true', help='add a column with the probe of every rule')
    sub.add_parser('probes', parents=[common], help='list the probes and their parameters')
    a = sub.add_parser('brief', parents=[common, select], help='which rules apply where on a map')
    a.add_argument('--input', required=True, help='the draft map')
    a.add_argument('--out', help='folder for rules_brief.md and .json (default: print the text)')
    a = sub.add_parser('check', parents=[common, select], help='which rules a map breaks')
    a.add_argument('--input', required=True, help='the map to check')
    a.add_argument('--draft', help='the draft it was laid out from: its subsystems are the units')
    a.add_argument('--previous', help='an earlier layout of the map: no placed node may have moved')
    a.add_argument('--out', help='folder for rules_report.md and .json')
    a.add_argument('--strict', action='store_true', help='exit code 1 if a must rule is broken')
    a = sub.add_parser('doc', help='read a document of the engine')
    a.add_argument('name', nargs='?')
    a.add_argument('--part', type=int, default=1)
    args = ap.parse_args(argv)

    if args.cmd == 'doc':
        return show_doc(args.name, args.part)
    if args.cmd == 'rules' and args.validate:
        try:
            rs = ruleset.load(args.rules, probes=P.PROBES, strict=False)
        except ruleset.RuleError as e:
            sys.exit(f'RULE ERROR: {e}')
        for problem in rs.problems:
            print('problem:', problem)
        print(f"{'rule set ok' if not rs.problems else 'rule set NOT ok'}: {len(rs)} rules "
              f"({sum(1 for r in rs if r.probe)} with a probe) in {', '.join(os.path.basename(p) for p in rs.files)}")
        sys.exit(1 if rs.problems else 0)
    rs = load_rules(args.rules)
    if args.cmd == 'probes':
        used = {}
        for r in rs:
            if r.probe:
                used.setdefault(r.probe, []).append(r.id)
        for name, p in P.PROBES.items():
            print(f"{name}  [{', '.join(used.get(name, [])) or 'no rule'}]")
            print(f'    {p.doc}')
            if p.defaults:
                print('    parameters: ' + ', '.join(f'{k}={v!r}' for k, v in p.defaults.items()))
        return
    chosen = rs.select(only=args.id, skip=args.skip, types=args.type, strengths=args.strength, groups=args.group)
    if args.cmd == 'rules':
        group = None
        for r in chosen:
            if (r.type, r.group) != group:
                group = (r.type, r.group)
                print(f'\n== {r.type}: {r.group}')
            probe = f'  [{r.probe}]' if args.bindings and r.probe else ''
            print(f'{r.id} ({r.strength}, {r.basis}) WHEN {r.when} DO {r.do}{probe}')
        print(f'\n{len(chosen)} of {len(rs)} rules')
    elif args.cmd == 'render':
        text = ruleset.to_markdown(rs, chosen, bindings=args.bindings)
        if args.out:
            with open(args.out, 'w', encoding='utf-8') as f:
                f.write(text)
            print(f'{len(chosen)} rules written to {args.out}')
        else:
            sys.stdout.write(text)
    elif args.cmd == 'brief':
        try:
            text, _ = report.brief(rs, args.input, args.out, chosen)
        except (OSError, ValueError) as e:
            sys.exit(f'ERROR: {e}')
        if args.out:
            print(f'rules that apply to this map, and where: {args.out}/rules_brief.md')
        else:
            sys.stdout.write(text)
    elif args.cmd == 'check':
        try:
            _, data, summary = report.check(rs, args.input, args.draft, args.previous, args.out, chosen)
        except (OSError, ValueError) as e:
            sys.exit(f'ERROR: {e}')
        print('\n'.join(summary))
        sys.exit(1 if args.strict and data['state'] == 'not_met' else 0)


if __name__ == '__main__':
    main(sys.argv[1:])
