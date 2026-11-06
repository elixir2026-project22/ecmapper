# escher_handlayout_kit

A reproducible toolkit to turn an **Escher draft map** into a hand-laid-out, smoothed, interconnected and packed map.
The creative step (where things go) is written as a **JSON layout spec** - by a person, by a chat-based AI agent
(`agent.py`, any model, **no code execution by the model**), or by a sandboxed coding agent following `PROMPT.md`.
It ships with ecmapper (`toolkits/escher_handlayout_kit`; select this folder as the toolkit in ecmapper's AI layout
panel) but runs on its own: it does not need ecmapper, only numpy, matplotlib and Pillow.

What the procedure produces (worked example: iML1515 central metabolism, 91 reactions, 10 subsystems):
1. each subsystem laid out by hand with textbook shapes (chains, branches, hourglasses, cycles as rings, lanes),
2. subsystems arranged in a biologically plausible flow (glycolysis above TCA, side modules), connected only through a short list of hub metabolites,
3. gentle Bezier curves instead of straight lines (curvature adjustable, 0.6 recommended),
4. labels placed so that lines do not hide them,
5. disconnected pieces packed into free space (rectangle area -24 % in the example).

## Three ways to use it

| you want | do this |
|---|---|
| let an **AI model** do the layout in a chat, safely | `python agent.py --input map.json --out run1 --provider manual` (copy/paste into any chat, no API) or `--provider openai/anthropic` - see **docs/AGENT.md** |
| write the layout yourself | `python run.py inspect ...`, edit a copy of `specs/layout_template.json`, `python run.py check-spec`, `python run.py build` - see **docs/GUIDE.md** |
| let a **coding agent** (with a Python sandbox) do everything | give it the kit and `PROMPT.md` |

## Contents

```
README.md                  this file
agent.py                   chat-based agent: model-agnostic, JSON in/out only (docs/AGENT.md)
run.py                     CLI: inspect | check-spec | build | validate
PROMPT.md                  prompt for a coding agent with a code sandbox (writes the JSON spec, runs run.py, iterates on the metrics)
selftest.py                quick regression test      tests/run_tests.py   full test suite (27 tests, ~100 s, offline)
requirements.txt           numpy, matplotlib, Pillow
schema/layout_spec.schema.json    JSON Schema of the layout spec
escher_handlayout/
  declspec.py              spec schema, strict validator, patch-merge, compiler (spec -> engine); no code execution
  engine.py                Layout engine: node splitting/merging, geometry, Bezier rules, ring arcs, label placement
  pack.py                  components + packing of disconnected pieces (occupancy grid, FFT search)
  render.py                matplotlib preview + Escher JSON export
  inspect_map.py           analysis of the draft map (subsystems, components, merge candidates, traps)
  qa.py                    structural validation + layout-quality metrics
  pipeline.py              build(): spec -> JSON + PNG + QA report
  agent/                   chat harness: llm.py (OpenAI-compatible, Anthropic, manual, scripted), session.py, jsonio.py, prompts/*.md
specs/
  iml1515_central_carbon.layout.json   worked example (complete JSON spec of the bundled map)
  layout_template.json                 empty skeleton
  legacy_python/                       old Python-code specs (execute code; disabled unless --unsafe-python-spec)
docs/GUIDE.md              spec format, layout recipes, curvature, how to read the QA metrics
docs/AGENT.md              the chat agent: modes, conversation, security model, troubleshooting
examples/                  the draft map + expected output (regression reference)
```

## Quick start

```bash
pip install -r requirements.txt
python selftest.py                                                         # SELFTEST PASSED
python run.py build --input examples/central_metabolism_escher.json \
                    --spec specs/iml1515_central_carbon.layout.json --out /tmp/result --curv 0.6
python agent.py --input examples/central_metabolism_escher.json --out /tmp/run1 --provider manual    # chat-based, copy/paste
```
Outputs: `<name>.json` (load it in https://escher.github.io or with the Escher Python package), `<name>_overview.png`,
`qa_report.md/json`. Options of `build`: `--curv 0..1` (0 gentle, 1 strong), `--no-pack`, `--no-png`, `--png-scale`.

## Reference results (bundled example, curvature 0.6)

| item | value |
|---|---|
| reactions / nodes | 91 / 495 -> 510 (duplicated nodes only; no reaction changed) |
| bounding box before / after packing | 81.2 M -> 61.7 M units^2 (about -24 %) |
| QA metrics | crossings 26, labels on lines 5, labels on nodes 0, label overlaps 0, validation errors 0 |
| determinism | byte-identical JSON across runs |
| chat agent, scripted model | reproduces the expected output exactly in 11 turns |

## Limitations (please read)

- The engine expects the **draft-map topology** (substrate -> multimarker -> midmarker -> multimarker -> product, as produced by
  `ecmapper.draft_map`). Maps with other topologies are reported by `inspect` as unsupported.
- The example spec uses the node ids of the bundled map; another map needs its own spec (written by you, by the agent, or by a coding agent).
- Validated structurally and with the kit's metrics, **not** loaded in Escher by the toolkit; label widths are estimated.
- The chat agent has been tested offline (scripted model, fake HTTP server, fuzzing) but **not with real language models yet**;
  prompts may need tuning per model (docs/AGENT.md section 6).
- Curvature 0.75 and higher produced overlapping curves at shared nodes in the example; 0.6 was preferred.
