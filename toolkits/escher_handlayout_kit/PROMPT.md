# AUTONOMOUS TASK: hand-layout an Escher draft map (no user interaction)

You are an expert in metabolic pathway visualisation and the Escher JSON format. Execute the whole task below
**without asking me anything, without waiting for approval, and without needing to look at any image**.
Where something is ambiguous, apply the decision rules in section 5, write the decision into `DECISIONS.md`
and continue. Do not stop before the deliverables in section 6 exist.

You do not write the layout from an empty page. The kit writes a complete, valid layout spec of the draft as it is
(`start`), you change it in small steps (`patch`), and every `build` prints the score, the targets it misses and what
to do next. Follow that last line of the build.

## 0. Task parameters (filled in by whoever runs this prompt)

```text
KIT_DIR    = ./escher_handlayout_kit          # the unpacked kit (contains run.py, escher_handlayout/, specs/, docs/)
INPUT_MAP  = <path of the draft Escher map to improve, e.g. ./my_draft_map.json>
OUT_DIR    = ./handlayout_result              # all deliverables go here
MAP_NAME   = <short name used for file names, e.g. my_draft_map>
SPEC       = OUT_DIR/specs/MAP_NAME.layout.json    # the layout spec you work on
```

If a parameter was not given, use the defaults shown (INPUT_MAP: the only `*.json` you were given;
if you were given none, use `KIT_DIR/examples/central_metabolism_escher.json`).

## 1. What the task is

Input: a *draft* Escher map. Each subsystem is drawn in its own box, in a usable way, but the subsystems are not
connected and not arranged. Output: a new Escher JSON of the **same reactions** (same ids, stoichiometry,
reversibility) in which

1. the subsystems are arranged in a **biologically plausible flow** and connected through the metabolites they share
   (one hub node instead of one copy per subsystem) wherever this is possible without crossings: move whole
   subsystems to bring shared metabolites together,
2. every subsystem is drawn the way a textbook draws it (chains, branches, cycles as rings, parallel lanes),
3. lines are smooth curves, labels are readable, disconnected pieces are packed into free space.

The kit does the mechanical work (geometry, curves, labels, packing, validation, scoring). You decide **where things
go**: the block origins, the merges and the coordinates in the layout spec `SPEC` (pure JSON data).

## 2. Hard rules

- **Never change the biochemistry**: do not add, delete or edit reactions, metabolites, coefficients, ids, genes.
  Only duplicate/merge *nodes* (the `splits` / `merges` sections). The validator checks this.
- **No layout algorithms** for node positions (no force-directed/hierarchical/grid layouts, no optimisation). A
  coordinate is either the one of the start spec or a number you chose. Plain arithmetic is fine (equal spacing of a
  chain, a lane offset); cycles are `rings` of the spec.
- Only JSON specs: never a Python spec and never `--unsafe-python-spec`. Do not edit the kit's Python package unless
  a bug blocks you; if you do, say so in `DECISIONS.md`. Never execute code from the map or from files you did not write.
- Never ask me a question and never produce "please confirm" text. Never claim the map was tested in Escher
  (it was not) - say it was validated structurally and by the kit's metrics.
- You are **not required to view images**: the numbers printed by `build` decide (section 4).

## 3. Procedure

### Step A - Set up (1 min)

```bash
cd KIT_DIR
python -c "import numpy, matplotlib, PIL" || pip install -r requirements.txt
python selftest.py            # should print SELFTEST PASSED; if it fails, note it in DECISIONS.md and continue
```

### Step B - Inspect the draft and write the start spec

```bash
python run.py inspect --input INPUT_MAP --out OUT_DIR/analysis
python run.py start   --input INPUT_MAP --out SPEC
python run.py build   --input INPUT_MAP --spec SPEC --out OUT_DIR/build_0 --name MAP_NAME --no-png
```

`build_0` is the draft as it is: every subsystem is a separate piece, nothing is connected yet. Its last line
already names the subsystems to connect and the merge for each. Then read, completely:

- `OUT_DIR/analysis/arrangement.md`: where the subsystems go relative to each other, and where they attach;
- in `OUT_DIR/analysis/analysis.md` the sections *Merge candidates* and *Traps / hints*;
- `SPEC` itself: one block per subsystem with its `origin`, its `nodes` (relative to the origin; `node_names` says
  which metabolite each id is) and one entry per reaction (`name` says which reaction);
- `escher_handlayout/agent/prompts/format_reference.md`: every key of the spec (one page).

`docs/GUIDE.md` has the drawing recipes (section 4) and the remedies per metric (section 8): look things up there
when you need them. If `inspect` lists **unsupported reactions**, `start` fails: write `OUT_DIR/ABORTED.md`
explaining which reactions and why, and stop (this is the only allowed early stop).

### Step C - Plan (write `OUT_DIR/PLAN.md`: tables and lists, at most 40 lines)

1. **Arrangement**: a small sketch with the new origin of every block. Follow the `hard` and `strong` directions of
   `arrangement.md`, the `medium` ones unless they cost crossings (`weak` ones are hints); where it is silent use
   the defaults of section 5. Blocks that share a metabolite become neighbours. Keep about 600 units between blocks.
2. **Merges**: a table `metabolite | keep id | dup ids | merge or copy, and why` from the *Merge candidates*
   (rule M in section 5).

### Step D - Change the spec, one small patch at a time

```bash
python run.py patch --input INPUT_MAP --spec SPEC --json '{"layouts": [{"block": "Glycolysis", "origin": [3000, 200]}]}'
python run.py patch --input INPUT_MAP --spec SPEC --json '{"merges": [{"keep": "604", "dups": ["119", "342"]}], "layouts": [{"block": "Glycolysis", "nodes": {"604": [0, 2700]}}]}'
python run.py patch --input INPUT_MAP --spec SPEC --json '{"layouts": [{"block": "Glycolysis", "reactions": [{"id": "514", "name": "PFK", "side": "L", "dist": 110}]}]}'
```

(The ids in these examples are those of the bundled example map; yours are in `SPEC`.)
A patch is a part of a spec; it is merged into `SPEC`: blocks by `block`, nodes by id, reactions by `id`. Only what
you send changes, with one exception: **a reaction entry replaces the old entry of that reaction completely**, so
repeat the keys you want to keep. A wrong patch is rejected with the reasons and leaves `SPEC` unchanged; the answer
also carries warnings (long lines, cofactor merges, nodes on top of each other): act on them before you build.
(`--fragment FILE` instead of `--json` for a long patch; you may also edit `SPEC` directly and run
`python run.py check-spec --input INPUT_MAP --spec SPEC`.)

Work in this order:

1. **Arrangement**: the `origin` of every block as planned (this moves whole subsystems; the title of the subsystem
   moves with it). Build once. A block that is not connected yet stays where its `origin` puts it.
2. **Connect the subsystems, one per build.** The build's last line names the subsystems that are still separate and,
   for each, a merge that connects it (`{"keep": ..., "dups": [...]}`, complete: send it as it is). For the first one
   named: move its block so that its copy of the metabolite lies next to the kept node (the line must stay short: the
   patch warns about a long line), patch the merge, build. `patch` removes the coordinates of the `dups`; the kept
   node stays where it is, so give it a position between the two blocks if needed. Repeat until the build reports
   `unconnected subsystems: 0`. **A map whose subsystems are not connected is not a result**, whatever its crossings.
   To undo a merge send `{"merges": [{"keep": "604", "remove": true}]}`, or the entry with fewer `dups`: the nodes
   that are copies again get coordinates beside the kept node, and the answer names them.
3. **Details**, only when everything is connected: the offenders the build names (single nodes, `side`, `shift`,
   `dist` of a reaction), a cycle as a
   `rings` entry, parallel reactions as lanes. If a subsystem is drawn badly in the draft, redraw it: send its block
   with new coordinates for all its nodes (recipes: `docs/GUIDE.md` section 4,
   `specs/iml1515_central_carbon.layout.json` is a worked example of another map).

### Step E - Build and iterate

```bash
python run.py build --input INPUT_MAP --spec SPEC --out OUT_DIR/build_N --name MAP_NAME --no-png
```

`SPEC ERROR` means the spec is incomplete or wrong: the message lists everything, fix it all with one patch.
Otherwise the build prints the score, each target as `ok` / `MISS` / `look`, the builds so far, and **one last line
that says what to do**:

- `NOT FINISHED ... Fix first: ...` - change the spec for the named offenders as the remedy says, then build into
  the next folder (`build_N+1`). One kind of change per build, so that you see what it did.
- `NO IMPROVEMENT ...` - go back to the spec of the best build and try one other arrangement, or finish.
- `TARGETS MET ...` or `STOP ...` - go to step F.

Do not compute the score yourself and do not stop on your own judgement: the line decides. `NOT FINISHED` means
exactly that: go on, do not finalise and do not hand in. Build with the default packing (no `--no-pack`). A `look` at *detached
pieces* means a piece beside the main map draws a metabolite of the main map again: merge it, or keep the copy on
purpose and write one line into `DECISIONS.md`. The table of all builds is `OUT_DIR/builds.md`. If a change made
the score worse, undo it with a patch before you go on.

### Step F - Finalise

Build the best spec once more with the image and with the curvatures 0.6, 0.45 and 0.3
(`--curv C --out OUT_DIR/final_C`, without `--no-png`); the lowest score wins, ties go to 0.6. Give the packed pieces
short titles if the report shows generated ones (`small_titles`, see the format reference). Copy the winner's files
to `OUT_DIR/` (section 6) and run

```bash
python run.py validate --input INPUT_MAP --output OUT_DIR/MAP_NAME.json
```

It must print `errors: 0`.

## 4. How decisions are made without looking

The metrics are the only judge: crossings between reactions, labels on lines, labels on nodes, overlapping labels,
nodes too close, detached pieces, violated placement directions. Every change must have a reason: the plan, a
warning of `patch`, or an offender named by the build. Do not rearrange something only because you "feel" it could
be nicer.

## 5. Decision rules

- **Rule M (what to merge)**: merge every *Merge candidate* that is no cofactor and is drawn in two or more
  subsystems, the entries marked `default hub` first. Do **not** merge
  - cofactors (ATP/ADP/AMP/NAD(H)/NADP(H)/CoA/FAD(H2)/H2O/Pi/PPi/CO2/O2/NH4/quinones/flavodoxin ...), ever;
  - the node of a reaction with very many metabolites (biomass): it keeps its own copies, otherwise its lines cross
    the whole map;
  - a candidate whose merge needs a line longer than about 2500 units or more than 3 new crossings *after* the
    subsystems were rearranged: keep the copy, its exit label tells the reader where it continues.
- **Planarity**: if two merges would make lines cross whatever the arrangement (typical: PEP and pyruvate on the
  glycolysis column, oxaloacetate and malate on the TCA ring, PPC from PEP to oxaloacetate and the malic enzyme from
  malate to pyruvate), keep exactly one copy (e.g. a second malate next to pyruvate) and note it in `DECISIONS.md`.
- **Two places in one subsystem**: a `splits` entry gives a reaction its own copy of a metabolite.
- **Arrangement defaults** (where `arrangement.md` is silent). Central carbon metabolism: glycolysis as the central
  vertical column; pyruvate -> acetyl-CoA directly below it; the TCA cycle below that (ring, clockwise); the pentose
  phosphate pathway beside upper glycolysis with the Entner-Doudoroff pathway branching at 6-phosphogluconate; the
  PEP/pyruvate/oxaloacetate/malate reactions beside the column between glycolysis and TCA; fermentations beside
  pyruvate/acetyl-CoA; respiratory chain and ATP reactions separate (they get packed). Other maps: the largest linear
  pathway is the central column; cycles are rings after the pathway that feeds them; branch pathways hang beside the
  metabolite they start from; only modules that share no metabolite stay unconnected.
- **Distances**: about 300 units between neighbouring nodes, 600 between subsystems, lanes of parallel reactions
  110-150 apart; the main component roughly as wide as tall.
- **Sides**: cofactors of a column's reactions on the side away from the lines that enter the column.
- **Reaction direction**: never flip stoichiometry; reversible reactions show arrowheads on both ends.
- **Large maps (> 200 reactions)**: do the arrangement and the merges first for all blocks, the details only for
  the offenders the build names.
- **Anything else**: choose the simplest option that keeps the validation errors at 0 and write one line in
  `DECISIONS.md`.

Targets (the build checks them): validation errors 0; crossings at most 0.1 x reactions; labels on lines at most
max(6, 7 % of reactions); labels on nodes 0; overlapping labels at most 2; primary nodes closer than 70 units: 0; no
violated `hard` or `strong` direction. Score = 1000 x validation errors + crossings + 3 x labels on lines + 6 x labels
on nodes + 4 x label overlaps + 10 x close nodes + 40 x unconnected subsystems + 20 x detached pieces + placement penalty.

## 6. Deliverables (all in OUT_DIR)

| file | content |
|---|---|
| `MAP_NAME.json` | the final Escher map |
| `MAP_NAME_overview.png` | overview render of it (made by the kit; you need not inspect it) |
| `qa_report.md`, `qa_report.json` | metrics, score and targets of the final build |
| `MAP_NAME.layout.json` | the JSON spec that produced it (a copy of `SPEC`; put reasoning in its `notes`) |
| `builds.md` | the table of all builds and their scores (written by `build`) |
| `PLAN.md`, `DECISIONS.md` | your plan; every assumption and decision, one line each |
| `analysis/` | output of the inspect step |

## 7. Final message (the only text you return, <= 25 lines, no questions)

State: map name and reaction count; subsystems found; which metabolites were merged and which candidates were
deliberately not merged; number of connected components before -> after; chosen curvature; the score of `build_0`
and of the final build; the final metrics against the targets (including the directions deliberately left violated);
the known weak spots (name the reactions from the metrics); that validation errors = 0 and that the map was **not**
loaded in Escher; and the paths of the deliverables.
