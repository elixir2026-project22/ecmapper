# AUTONOMOUS TASK: hand-layout an Escher draft map (no user interaction)

You are an expert in metabolic pathway visualisation and the Escher JSON format. Execute the whole task below
**without asking me anything, without waiting for approval, and without needing to look at any image**.
Where something is ambiguous, apply the decision rules in section 5, write the decision into `DECISIONS.md`
and continue. Do not stop before the deliverables in section 6 exist.

**The design of the map is yours.** You know how this metabolism is drawn in textbooks and in published maps, and
that knowledge decides where the subsystems go and how they connect (Step C). You do not write the layout from an
empty page, though: the kit writes a complete, valid layout spec of the draft as it is (`start`), you change it in
small steps (`patch`), and every `build` reports the defects of the drawing that are left (crossings, labels on lines,
long lines ...) and what to repair next. The score measures defects of the drawing, not whether the map makes
biochemical sense: use it to repair your design, never to choose it.

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
  chain, a lane offset); cycles are `rings` of the spec. The kit's own `tidy` (Step D 3) is allowed: it repairs
  details of the drawing and leaves your arrangement as it is.
- Only JSON specs: never a Python spec and never `--unsafe-python-spec`. Do not edit the kit's Python package unless
  a bug blocks you; if you do, say so in `DECISIONS.md`. Never execute code from the map or from files you did not write.
- Never ask me a question and never produce "please confirm" text. Never claim the map was tested in Escher
  (it was not) - say it was validated structurally and by the kit's metrics.
- You can do the task **without viewing images**: the numbers printed by `build` show the defects of the drawing
  (section 4). If you can view images, look at the overview of a build (build without `--no-png`) when the
  arrangement is done and once more before you finalise, and check that the map reads like your design.

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

`start` does a first arrangement for you: it takes the drawing of every subsystem of the draft (a cycle is redrawn as
a ring, a wobbling chain as a straight column), puts the subsystems where the placement catalogue wants them and
connects them through shared metabolites (`merges`). Its answer lists what it redrew, what it merged and which
metabolites are still drawn more than once. `build_0` is this first arrangement: a connected map to start from.
It is a floor, not a design: the subsystems still have the shapes of the draft, the map is larger than it has to be
and some lines are long. A good map needs subsystems that are **redrawn** (Step D 1); moving blocks alone does not
get there. Then read, completely (`arrangement.md` comes later, in Step C 2):

- in `OUT_DIR/analysis/analysis.md` the sections *Merge candidates* and *Traps / hints*;
- `SPEC` itself: one block per subsystem with its `origin`, its `nodes` (relative to the origin; `node_names` says
  which metabolite each id is) and one entry per reaction (`name` says which reaction);
- `escher_handlayout/agent/prompts/format_reference.md`: every key of the spec (one page).

`docs/GUIDE.md` has the drawing recipes (section 4) and the remedies per metric (section 8): look things up there
when you need them. If `inspect` lists **unsupported reactions**, `start` fails: write `OUT_DIR/ABORTED.md`
explaining which reactions and why, and stop (this is the only allowed early stop).

### Step C - Design (write `OUT_DIR/PLAN.md`: tables and lists, at most 50 lines)

1. **Your design, from what you know.** Do this before you read `arrangement.md`. The subsystems and their reactions
   are in `analysis.md` and `SPEC`. Write down how a biochemist would draw this map: which pathway is the backbone and
   in which direction the flux runs through it; which pathway feeds which, and through which metabolite; which
   metabolites are the hubs (drawn once, with their reactions around them); for every side pathway where it belongs and
   why (what it starts from, ends in, bypasses or reverses); which reactions run in parallel; which subsystems have
   nothing to do with the rest and stay apart. Use the picture you know from textbooks and from published maps of this
   metabolism (section 5 has the usual picture of central carbon metabolism).
2. **Cross-check with the catalogue.** Now read `OUT_DIR/analysis/arrangement.md` completely: where a catalogue of
   placement rules puts the subsystems relative to each other, and where they attach. Take over what you had not
   thought of. Where a direction contradicts your design, decide: follow it, or keep your design and write one line
   with the biochemical reason into `DECISIONS.md` (the build will report the direction as violated; that is accepted,
   see Step E). Overrule a `hard` or `strong` direction only if you are sure; `medium` and `weak` ones are hints.
3. **Geometry**: a table `block | keep, move or redraw | where | why`: which blocks stay as `start` drew them, which
   you move (new `origin`), which you redraw (new coordinates for their nodes). The backbone is one straight column
   with 350 to 450 units between consecutive metabolites; cycles are rings; side pathways are short rows or columns
   next to the metabolite they start from. Keep about 300 units free between neighbouring blocks. **Give hubs room**:
   where many reactions meet (pep, pyruvate, acetyl-CoA, oxaloacetate) the drawing needs more space, not less. Two
   metabolites that several parallel reactions join are about 150 units per reaction apart (pep and pyruvate: 800 or
   more), the next metabolites of the backbone move along, and the pathways around a hub leave it in different
   directions. A small map is no aim: the build counts `crowded reactions`.
4. **Merges**: a table `metabolite | keep id | dup ids | merge or copy, and why` for the metabolites that `start`
   reports as still drawn more than once (rule M in section 5).

If you cannot design the map from your own knowledge, take the arrangement of `arrangement.md` and of section 5 as
your design and follow the last line of every build.

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

1. **Redraw and arrange, one block per build.** Redraw a subsystem whose drawing is sprawling, irregular or far
   from the metabolites it shares: send its block with new coordinates for all its nodes (and `side` / `shift` of its
   reactions where they change), so that it has the shape of your design and its shared metabolites lie on the nodes
   they are merged into. Recipes: `docs/GUIDE.md` section 4; `specs/iml1515_central_carbon.layout.json` is a worked
   example of another map. Move a block with its `origin` (the title moves with it). The build's `long lines` and
   `area excess` say where this is needed. If you cannot work out the coordinates of a whole block, do not redraw it:
   move blocks and nodes instead. A block that is not connected stays where its `origin` puts it.
2. **Connect what is still separate, one subsystem per build.** `start` normally leaves nothing separate; if the
   build reports unconnected subsystems, its last line names them and,
   for each, a merge that connects it (`{"keep": ..., "dups": [...]}`, complete: send it as it is). For the first one
   named: move its block so that its copy of the metabolite lies next to the kept node (the line must stay short: the
   patch warns about a long line), patch the merge, build. `patch` removes the coordinates of the `dups`; the kept
   node stays where it is, so give it a position between the two blocks if needed. Repeat until the build reports
   `unconnected subsystems: 0`. **A map whose subsystems are not connected is not a result**, whatever its crossings.
   To undo a merge send `{"merges": [{"keep": "604", "remove": true}]}`, or the entry with fewer `dups`: the nodes
   that are copies again get coordinates beside the kept node, and the answer names them.
3. **Details**, only when everything is connected. First let the kit tidy the drawing:

   ```bash
   python run.py tidy --input INPUT_MAP --spec SPEC
   ```

   It puts metabolites that are almost above or beside each other exactly into a column or row, turns a metabolite
   that only one reaction uses around that reaction, puts a lane of parallel reactions on the other side of its line,
   and chooses `side` and `dist` of the reactions so that the cofactors lie in free space. It keeps only what makes
   the metrics better, changes neither the place of a block nor a merge, and takes about a minute and a half
   (`--block NAME` for one block). Run it when the arrangement is done and again after you redrew or moved a block,
   then build. Do by hand only what the build still names after that: single nodes, `side`, `shift`, `dist` of a
   reaction, a cycle as a `rings` entry, parallel reactions as lanes.

   `tidy --copies` goes one step further: a reaction that reaches a metabolite of three or more reactions over a long
   line gets its own copy of it beside its other metabolites, where that saves more than 3 points. A copy is a
   decision about the design (rule M in section 5): use it when you would keep such copies anyway, read which ones
   it made, and note them in `DECISIONS.md` or undo the ones you do not want.

### Step E - Build and iterate

```bash
python run.py build --input INPUT_MAP --spec SPEC --out OUT_DIR/build_N --name MAP_NAME --no-png
```

`SPEC ERROR` means the spec is incomplete or wrong: the message lists everything, fix it all with one patch.
Otherwise the build prints the score, each target as `ok` / `MISS` / `look`, the builds so far, and **one last line
that says what to do**:

- `NOT FINISHED ... Fix first: ...` - change the spec for the named offenders as the remedy says, then build into
  the next folder (`build_N+1`). One kind of change per build, so that you see what it did.
- `NO IMPROVEMENT ...` (three changes since the best build were worse) - go back to the spec of the best build and
  try one other arrangement, or finish.
- `TARGETS MET ...` or `STOP ...` - go to step F.

Do not compute the score yourself. The last line is advice on the drawing: `NOT FINISHED` names a defect that is
still there, and normally you repair it and build again; do not stop while the build names defects that you can repair
within your design. You may leave a target missed when its remedy would spoil the design: a direction you overruled
in Step C, a crossing that only a detour would remove. Write one line into `DECISIONS.md` and go on with the next
`MISS`. **Never meet a target with a trick**: no copy of a metabolite (`splits`, an undone merge) whose only purpose
is to satisfy a direction, and no pathway moved away from where it belongs to save a point of the score. The map you
hand in is the one you judge the best map of this metabolism; if it is not the build the last line names, say why in
`DECISIONS.md`. Build with the default packing (no `--no-pack`). A `look` at *detached
pieces* means a piece beside the main map draws a metabolite of the main map again: merge it, or keep the copy on
purpose and write one line into `DECISIONS.md`. The table of all builds is `OUT_DIR/builds.md`. If a change made
the score worse, undo it with a patch before you go on.

### Step F - Finalise

Build the spec of the map you judge best once more with the image and with the curvatures 0.6, 0.45 and 0.3
(`--curv C --out OUT_DIR/final_C`, without `--no-png`); the lowest score wins, ties go to 0.6. Give the packed pieces
short titles if the report shows generated ones (`small_titles`, see the format reference). Copy the winner's files
to `OUT_DIR/` (section 6) and run

```bash
python run.py validate --input INPUT_MAP --output OUT_DIR/MAP_NAME.json
```

It must print `errors: 0`.

## 4. Design and metrics

Two things decide, each in its own field. **Where the subsystems go and how they connect** is decided by your design
(Step C), that is by biochemistry: what feeds what, what a reader expects next to what. **How well it is drawn** is
measured by the metrics: crossings between reactions, labels on lines, labels on nodes, overlapping labels, nodes too
close, long lines, the area, detached pieces. The metrics cannot see whether a map makes sense: a map with a low
score whose pathways sit where no biochemist would look for them is a bad map, and a design is not given up for a few
points. So a change of the arrangement needs a biochemical reason (your plan), and a change of the drawing (a node, a
`side`, a lane) needs a warning of `patch` or an offender named by the build. Without images the numbers are your
only eyes on the drawing: do not tidy what they do not report.

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
- **The usual picture** (your own knowledge of the organism and of the map comes first). Central carbon metabolism: glycolysis as the central
  vertical column; pyruvate -> acetyl-CoA directly below it; the TCA cycle below that (ring, clockwise); the pentose
  phosphate pathway beside upper glycolysis with the Entner-Doudoroff pathway branching at 6-phosphogluconate; the
  PEP/pyruvate/oxaloacetate/malate reactions beside the column between glycolysis and TCA; fermentations beside
  pyruvate/acetyl-CoA; respiratory chain and ATP reactions separate (they get packed). Other maps: the largest linear
  pathway is the central column; cycles are rings after the pathway that feeds them; branch pathways hang beside the
  metabolite they start from; only modules that share no metabolite stay unconnected.
- **Distances**: about 300 units between neighbouring nodes, 600 between subsystems, lanes of parallel reactions
  110-150 apart; more around hubs (Step C 3); the main component roughly as wide as tall.
- **Hubs and copies**: a hub is drawn once for the pathways of the main flow. A side pathway that would be the
  seventh or later reaction at a hub may keep its own copy of the metabolite beside it, if the hub gets too crowded
  otherwise; note it in `DECISIONS.md`.
- **Sides**: cofactors of a column's reactions on the side away from the lines that enter the column.
- **Reaction direction**: never flip stoichiometry; reversible reactions show arrowheads on both ends.
- **Large maps (> 200 reactions)**: do the arrangement and the merges first for all blocks, the details only for
  the offenders the build names.
- **Anything else**: choose the simplest option that keeps the validation errors at 0 and write one line in
  `DECISIONS.md`.

Targets (the build checks them): validation errors 0; crossings at most 0.1 x reactions; labels on lines at most
max(6, 7 % of reactions); labels on nodes 0; overlapping labels at most 2; primary nodes closer than 70 units: 0;
long lines (a primary metabolite more than 800 units from the middle of its reaction) at most max(4, 10 % of
reactions); crowded reactions (4 or more other reaction middles within 250 units) at most max(2, 3 % of reactions); area excess at most 40 (the map's area is at most 0.80 million square units per reaction); no violated
`hard` or `strong` direction. Score = 1000 x validation errors + crossings + 3 x labels on lines + 6 x labels on nodes
+ 4 x label overlaps + 10 x close nodes + 3 x long lines + area excess (1 point per 0.01 million square units per
reaction above 0.40) + 3 x crowded reactions + 40 x unconnected subsystems + 20 x detached pieces + placement penalty.

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
