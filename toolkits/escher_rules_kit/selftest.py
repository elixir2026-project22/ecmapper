#!/usr/bin/env python3
"""Quick test: the bundled rule set is valid, agrees with its Markdown text, and every probe runs on the example map
of escher_handlayout_kit (if that kit is beside this one). The complete suite is `python tests/run_tests.py`."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from escher_rules import probes, report, ruleset  # noqa: E402

ok = True
try:
    rs = ruleset.load(probes=probes.PROBES)
except ruleset.RuleError as e:
    print('FAIL:', e)
    sys.exit(1)
print(f'{len(rs)} rules, {sum(1 for r in rs if r.probe)} with a probe, in {len(rs.files)} files')
with open(os.path.join(HERE, 'rules', 'metabolic_map_rules.md'), encoding='utf-8') as f:
    text = f.read()
for r in rs:
    if r.source == 'metabolic_map_rules.yaml' and ruleset.rule_row(r) not in text:
        print(f'FAIL: rule {r.id} differs between metabolic_map_rules.yaml and metabolic_map_rules.md')
        ok = False
example = os.path.join(os.path.dirname(HERE), 'escher_handlayout_kit', 'examples', 'central_metabolism_escher.json')
if os.path.isfile(example):
    _, data, summary = report.check(rs, example)
    print(summary[0])
    if not any(r['state'] == 'broken' for r in data['rules']) or not any(r['state'] == 'met' for r in data['rules']):
        print('FAIL: the check of the example draft finds nothing')
        ok = False
else:
    print('escher_handlayout_kit is not beside this kit: the layout commands will not work, the example was not checked')
print('SELFTEST', 'PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
