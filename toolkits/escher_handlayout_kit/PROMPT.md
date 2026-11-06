# AUTONOMOUS TASK: hand-layout an Escher draft map (no user interaction)

You are an expert in metabolic pathway visualisation and the Escher JSON format. Execute the whole task below
**without asking me anything, without waiting for approval, and without needing to look at any image**.
Where something is ambiguous, apply the decision rules in section 5, write the decision into `DECISIONS.md`
and continue. Do not stop before the deliverables in section 6 exist.

## 0. Task parameters (filled in by whoever runs this prompt)

```text
KIT_DIR    = ./escher_handlayout_kit          # the unpacked kit (contains run.py, escher_handlayout/, specs/, docs/)
INPUT_MAP  = <path of the draft Escher map to improve, e.g. ./my_draft_map.json>
OUT_DIR    = ./handlayout_result              # all deliverables go here
MAP_NAME   = <short name used for file names, e.g. my_draft_map>
```

If a parameter was not given, use the defaults shown (INPUT_MAP: the only `*.json` you were given;
if you were given none, use `KIT_DIR/examples/central_metabolism_escher.json`).

## 1. What the task is

Input: a *draft* Escher map (auto-generated, with a poor layout, subsystems scattered or overlapping).
Output: a new Escher JSON of the **same reactions** (same ids, stoichiometry, reversibility) in which

1. every subsystem is laid out **by hand** the way a textbook/biochemist would draw it,
2. the subsystems are arranged in a **biologically plausible flow** and interconnected through a **few** shared
   key metabolites (hub nodes), with as few crossings as possible,
3. all segments are smooth **Bezier curves** (gentle, no overlapping hooks),
4. reaction labels are not hidden behind lines,
5. disconnected pieces are **packed** into free space so the overall rectangle is as small as possible.

The kit does the mechanical work (geometry, curves, labels, packing, validation, rendering). **You provide the
creative part: a JSON layout spec with hand-chosen coordinates** (`specs/<MAP_NAME>.layout.json`). A spec is pure data
(schema in `schema/layout_spec.schema.json`); never write Python specs and never use `--unsafe-python-spec`.

## 2. Hard rules

- **No layout algorithms** for node positions: no force-directed/hierarchical/grid/tree layouts, no networkx/Tulip/
  graphviz positions, no random or optimisation-based placement. Every primary-node coordinate in the spec is a number
  you chose. Plain arithmetic helpers are fine (e.g. points of a regular polygon for a cycle, equal spacing of a chain,
  a lane offset).
- **Never change the biochemistry**: do not add, delete or edit reactions, metabolites, coefficients, ids, genes.
  Only duplicate/merge *nodes* (that is what the `splits` / `merges` sections do). The validator checks this.
- Do not edit the kit's Python package unless a bug blocks you; if you do, say so in `DECISIONS.md`. Never execute code that comes from the map or from files you did not write.
- Never ask me a question and never produce "please confirm" text. Never claim the map was tested in Escher
  (it was not) - say it was validated structurally and by the kit's metrics.
- You are **not required to view images**. The decision whether a layout is good is made by the numeric metrics
  in `qa_report.md` (section 4). If your environment can show images you may glance at the overview PNG, but the
  metrics decide.

## 3. Procedure

### Step A - Set up (1 min)

```bash
cd KIT_DIR
python -c "import numpy, matplotlib, PIL" || pip install -r requirements.txt
python selftest.py            # should print SELFTEST PASSED; if it fails, note it in DECISIONS.md and continue
```

### Step B - Inspect the draft map

```bash
python run.py inspect --input INPUT_MAP --out OUT_DIR/analysis
```

Read completely: `OUT_DIR/analysis/analysis.md`, `docs/GUIDE.md`, `escher_handlayout/agent/prompts/format_reference.md`
(the spec format), `specs/layout_template.json` and `specs/iml1515_central_carbon.layout.json` (the worked example - copy
its style and its recipes).
If `analysis.md` lists **unsupported reactions** (no multimarker topology), the engine cannot lay out this map:
write `OUT_DIR/ABORTED.md` explaining which reactions and why, and stop (this is the only allowed early stop).

### Step C - Plan (write `OUT_DIR/PLAN.md`, concise)

1. **Subsystems**: use the title labels / suggested subsystems of `analysis.md`; one zone per subsystem. If the map has
   no title labels, name each component after its dominant reactions/pathway (use your biochemical knowledge).
2. **Global arrangement**: apply section 5 rules. Draw it as a small ASCII sketch with approximate coordinates.
3. **Hub merges**: the exact list of `merges` entries (from `merge_candidates`, section 5 rule M).
4. **Per subsystem**: node coordinates (primary nodes only) and, per reaction, the `side`, `shift`, `axis`/`mid` or
   `arc` you intend, following the recipes in `docs/GUIDE.md` section 4.

### Step D - Write the spec

Create `KIT_DIR/specs/MAP_NAME.layout.json` from the template: `zones` for every subsystem -> `defaults.bend` ->
`splits` and `merges` -> one `layouts` block per subsystem with `nodes` (primary nodes only) and a `reactions` entry for
**every** reaction (each exactly once) -> `titles` for every zone that keeps at least one reaction in the main connected
component. Use ids exactly as in `analysis.md` (strings); after a merge do not use the deleted ids any more; nodes created by
`splits` are referred to by their `@alias`. Check quickly with
`python run.py check-spec --input INPUT_MAP --spec specs/MAP_NAME.layout.json`.

### Step E - Build and iterate (autonomous loop, at most 8 builds)

```bash
python run.py build --input INPUT_MAP --spec specs/MAP_NAME.layout.json --out OUT_DIR/build_N --curv 0.6 --name MAP_NAME
```

- `SPEC ERROR` -> the message lists every unplaced reaction/primary node; fix them all in one edit.
- Otherwise read `OUT_DIR/build_N/qa_report.md`. Compute the score

  `score = 1000*validation_errors + crossings + 3*labels_on_lines + 6*labels_on_nodes + 4*label_overlaps + 10*close_nodes`

  and keep a table `build_N -> score` in `DECISIONS.md`. Edit the spec to reduce the *named* offenders
  (`worst_crossing_pairs`, `labels_on_lines_list`, ... - remedies in `docs/GUIDE.md` section 8) and rebuild.
  Keep the best build. **Stop** when all targets hold, or after 8 builds (then take the best), or when two
  consecutive builds do not improve the score.

  Targets: `validation_errors = 0`; `labels_on_nodes = 0`; `label_label_overlaps <= 2`; `primary_nodes_closer_than_70 = 0`;
  `labels_on_lines <= max(6, 7% of reactions)`; `crossings <= 0.35 x number of reactions`.
  (Reference for the bundled example with 91 reactions: crossings 26, labels_on_lines 5, all others 0.)
- After the layout has converged, try `--curv 0.3` and `--curv 0.45` with the final spec (cheap) and keep the curvature
  with the lowest score; ties go to **0.6**.
- After the first successful build look at the `Packing` section of the report: pieces without an explicit
  `L.small_titles` entry got a generated title. Replace the generated titles by short meaningful ones
  (a `small_titles` entry: `{"reaction": <bigg_id of one of its reactions>, "zone": ..., "text": ..., "reuse_label": null}`) and rebuild once.

### Step F - Finalise

Copy the best build's files to `OUT_DIR/` (see section 6), run

```bash
python run.py validate --input INPUT_MAP --output OUT_DIR/MAP_NAME.json
```

and make sure it prints `errors: 0`.

## 4. How decisions are made without looking

The metrics in `qa_report.md` are the only judge: crossings between reactions, labels lying on lines, labels on
other nodes, overlapping labels, nodes too close, bounding-box area. Every edit must have a reason taken from the
named offenders. Do not rearrange a subsystem that has no offenders just because you "feel" it could be nicer.

## 5. Decision rules for ambiguity

- **Subsystem split**: title labels first; else connected components (`analysis.md`); a component that is much
  smaller than 5 reactions and clearly belongs to a neighbouring pathway (shares metabolite names) joins that zone.
- **Global arrangement** (central carbon metabolism, the default): glycolysis as the central vertical column;
  pyruvate -> acetyl-CoA directly below it; TCA cycle (ring, clockwise, citrate on top) below that so that glycolysis
  is above TCA; pentose phosphate pathway on one side of upper glycolysis (branching at G6P, returning F6P and G3P
  through the transketolase/transaldolase boxes) with the Entner-Doudoroff pathway branching at 6-phosphogluconate;
  the PEP/pyruvate/oxaloacetate/malate block on the side between glycolysis and TCA; fermentations beside
  pyruvate/acetyl-CoA/acetaldehyde; respiratory chain / ATP hub separate (they get packed).
  **Other maps**: the largest linear pathway becomes the central vertical column; cycles become rings placed after
  the pathway that feeds them; branch pathways hang to the side of the metabolite they start from; modules that share
  no key metabolite are left unconnected and will be packed. Keep the main component roughly as wide as tall (aspect
  between 0.7 and 1.4); keep >= 300 units between nodes and >= 600 between subsystems.
- **Rule M (which metabolites to merge)**: from `merge_candidates` take (a) every entry flagged `default hub` that
  appears in >= 2 different subsystems, plus (b) at most 4 further non-cofactor metabolites that appear in >= 2
  subsystems and are used by >= 3 reactions in total. Never merge cofactors (ATP/ADP/AMP/NAD(H)/NADP(H)/CoA/FAD(H2)/
  H2O/Pi/PPi/CO2/O2/NH4/quinones/flavodoxin...). If merging a candidate would force a line longer than ~2500 units
  or more than 3 new crossings, do not merge it - keep the copy and rely on the exit label. Use a `splits` entry when one
  metabolite is needed at two places of the same subsystem.
- **Parallel reactions** (`duplicate_reaction_bigg_ids`, `parallel_reaction_families`): lanes 110-150 apart; variants
  with different cofactor pairs stacked on one lane (recipe in the guide).
- **Reaction direction**: never flip stoichiometry; reversible reactions show arrowheads on both ends.
- **Curvature**: 0.6 unless step E finds a better one.
- **Large maps (> 200 reactions)**: `build` stops with a SPEC ERROR until every reaction and primary node is placed, so
  first place everything coarsely (subsystem by subsystem, using the recipes), get a first complete build, then refine
  using the metrics. Use `--no-png` while iterating to save time.
- **Anything else**: choose the simplest option that keeps the validation errors at 0 and write one line in `DECISIONS.md`.

## 6. Deliverables (all in OUT_DIR)

| file | content |
|---|---|
| `MAP_NAME.json` | the final Escher map |
| `MAP_NAME_overview.png` | overview render of it (made by the kit; you need not inspect it) |
| `qa_report.md`, `qa_report.json` | metrics and validation of the final build |
| `MAP_NAME.layout.json` | the JSON spec that produced it (hand-chosen coordinates; put reasoning in its `notes`) |
| `PLAN.md`, `DECISIONS.md` | your plan, the score table, every assumption and decision |
| `analysis/` | output of the inspect step |

## 7. Final message (the only text you return, <= 25 lines, no questions)

State: map name and reaction count; subsystems found; which metabolites were merged and which candidates were
deliberately not merged; node count before -> after (nodes increase because of duplicates); chosen curvature;
bounding-box area before -> after packing and final size; the final metrics against the targets; the known weak
spots (name the reactions from the metrics); that validation errors = 0 and that the map was **not** loaded in Escher;
and the paths of the deliverables.
