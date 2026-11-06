#!/usr/bin/env python3
"""Quick regression test: rebuild the bundled example from its JSON spec and compare with examples/expected_output.
(The complete suite - schema, security, agent, HTTP clients, fuzzing - is `python tests/run_tests.py`.)"""
import sys, os, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from escher_handlayout import build
from escher_handlayout._io import read_json

inp = os.path.join(HERE, 'examples', 'central_metabolism_escher.json')
spec = os.path.join(HERE, 'specs', 'iml1515_central_carbon.layout.json')
exp = os.path.join(HERE, 'examples', 'expected_output')
with tempfile.TemporaryDirectory() as tmp:
    rep = build(inp, spec, tmp, curv=0.6, png=False, name='out')
    new = read_json(os.path.join(tmp, 'out.json'))[1]
ref = read_json(os.path.join(exp, 'central_metabolism_escher_final_packed.json'))[1]
ok = True
def sig(m):
    out = {}
    for rid, r in m['reactions'].items():
        pts = sorted({(m['nodes'][n].get('bigg_id', m['nodes'][n]['node_type']), round(m['nodes'][n]['x'], 3), round(m['nodes'][n]['y'], 3))
                      for s in r['segments'].values() for n in (s['from_node_id'], s['to_node_id'])})
        out[rid] = (pts, round(r['label_x'], 3), round(r['label_y'], 3))
    return out
if sig(new) != sig(ref): print('FAIL: geometry differs from expected output'); ok = False
refm = read_json(os.path.join(exp, 'qa_report.json'))['metrics']
for k in ('line_crossings', 'labels_on_lines', 'labels_on_nodes', 'label_label_overlaps', 'width', 'height'):
    if rep['metrics'][k] != refm[k]: print(f'FAIL: metric {k}: {rep["metrics"][k]} != {refm[k]}'); ok = False
if rep['validation_errors']: print('FAIL: validation errors', rep['validation_errors'][:5]); ok = False
print('SELFTEST', 'PASSED' if ok else 'FAILED'); sys.exit(0 if ok else 1)
