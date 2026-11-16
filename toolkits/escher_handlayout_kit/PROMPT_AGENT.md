# Lay out an Escher draft map with a coding agent

You are a coding agent: you can run commands, write and run scripts, and look at images. You are also an expert in
metabolic pathways and how textbooks and published maps draw them. This file is your procedure. (`PROMPT.md` is the
procedure for a chat model that can do none of that; do not follow it.)

Work without asking questions. Where something is ambiguous, decide, write one line into `DECISIONS.md` and go on.

## 1. The task

Input: a *draft* Escher map. Each subsystem is drawn in its own box, in a usable way, but the subsystems are neither
connected nor arranged. Output: an Escher map of the **same reactions** (same ids, stoichiometry, reversibility) that
a biochemist would call a good drawing of this metabolism:

1. the pathways sit where a reader expects them and flow into each other through the metabolites they share;
2. every pathway has the shape a textbook gives it: a straight chain, a ring, a branch, parallel lanes;
3. there is room where many reactions meet, the lines are smooth, and every label can be read.

**You design the map and write its layout yourself.** The kit turns your layout (a JSON *spec*: where each primary
metabolite is, on which side each reaction hangs its cofactors) into the map: curves, cofactors, labels, packing,
validation, metrics, pictures.

## 2. Rules

- **Never change the biochemistry**: no reaction, metabolite, coefficient, id or gene is added, deleted or edited.
  Only *nodes* are duplicated or merged (`splits`, `merges`). The validator checks this.
- **You compute coordinates, you do not search for them.** Write a script that produces the spec from your design:
  a column with a spacing, a ring with a centre and a radius, a fan around a hub, a mirrored block
  (`escher_handlayout/shapes.py`). No force-directed, hierarchical or other global layout algorithm, and no optimiser
  over node positions: where things go is your decision, the arithmetic is the script's.
- Only JSON specs for the kit (never `--unsafe-python-spec`). Do not edit the kit; if a bug in it blocks you, work
  around it and say so in `DECISIONS.md`. Do not edit the draft map.
- Never claim the map was tested in Escher: it is validated structurally and by the kit's metrics, and you looked at
  its pictures.

## 3. The commands

`RUN` stands for the command your task description gives for the kit (e.g. `python run_toolkit.py`; in the kit's own
folder `python run.py`).

```text
RUN inspect --input MAP --out analysis          facts about the draft (analysis.md) and where a catalogue of placement
                                                rules puts the subsystems (arrangement.md)
RUN start --input MAP --out reference.json      a complete, valid spec of the draft: a reference for the format, the
                                                ids, the merges and the reaction entries. Not your layout.
RUN script make_spec.py                         runs your script with the kit importable
RUN check-spec --input MAP --spec SPEC          fast check without a build
RUN build --input MAP --spec SPEC --out build_N --max-builds 100     the map, an overview picture, score and targets
RUN view --input MAP --spec SPEC --out views --around pyr_c accoa_c  pictures of 1800 x 1800 units around metabolites
                                                or reactions, the defects marked (--radius, --no-marks; without
                                                --around: the overview with marks)
RUN tidy --input MAP --spec SPEC                mechanical repairs that keep your design: rows and columns, loose ends,
                                                lanes, cofactors; only what makes the metrics better (about 90 s)
RUN validate --input MAP --output FINAL.json    must print errors: 0
```

Read once: `escher_handlayout/agent/prompts/format_reference.md` (every key of a spec, one page),
`escher_handlayout/shapes.py` (the helpers for your script), and `docs/GUIDE.md` section 4 (how to draw a chain, a
ring, lanes, an hourglass reaction). `specs/iml1515_central_carbon.layout.json` is a worked example of another map.

## 4. Procedure

### A. Understand the draft

Run `inspect` and `start`. Read `analysis/analysis.md` (subsystems, reactions with their ids, *Merge candidates*,
*Traps / hints*) and skim `reference.json`. Build `build_0` from `reference.json` and look at its overview: this is
what the kit makes without you. It is connected and valid, and its shapes are those of the draft.

### B. Design (`PLAN.md`, at most 60 lines)

1. **From what you know**, before reading `arrangement.md`: how is this metabolism drawn? Which pathway is the
   backbone and in which direction does it run; what feeds what, through which metabolite; which metabolites are hubs;
   where does each side pathway belong and why (what it starts from, ends in, bypasses or reverses); which reactions
   run in parallel; what has nothing to do with the rest and stays apart.
2. **Cross-check** with `analysis/arrangement.md`. Take over what you had not thought of. Where a direction
   contradicts your design, follow it or keep your design with one line of biochemical reason in `DECISIONS.md`.
3. **A sketch with numbers**: for every block its shape (column, row, ring, fan), its place, its size. Rules of thumb:
   - consecutive metabolites of a chain 350 to 450 units apart, neighbouring blocks 300 or more apart;
   - **hubs get room**: two metabolites joined by several parallel reactions are about 150 units per reaction apart
     (`shapes.hub_gap`; pep and pyruvate: 800 or more), and the pathways around a hub leave it in different
     directions. A small map is no aim;
   - a ring's radius from its number of nodes (`shapes.ring_radius`), reactions as arcs;
   - a metabolite of the main flow is drawn once. A side pathway may keep its own copy of a hub metabolite beside it
     when the hub would get crowded or a line would cross the map; never merge cofactors.
4. **Merges and copies**: which shared metabolites become one node (`merges`), which stay copies and why.

### C. Write the layout as a script (`make_spec.py`)

The script builds the whole spec and writes it (`shapes.save`): zones, merges, splits, one block per subsystem with
the coordinates of its primary metabolites and one entry per reaction, the titles. Take ids, zone names, merges and
reaction entries from `reference.json` where they fit (`shapes.load`), look ids up by name with `shapes.Draft`, and
compute every coordinate from your sketch with `column`, `row`, `line`, `ring`, `fan`, `lanes`, `mirror`, `rotate`,
`move`. Keep the numbers of the design (origins, spacings, radii) as named constants at the top: you will change them.

Run it, then `check-spec`, then `build`. A `SPEC ERROR` lists everything that is wrong; fix the script, not the JSON.

### D. Look, then change the script (repeat)

After **every** build look at the overview, and at a `view` of every hub and of every place the build names. Judge in
this order:

1. Does the map read like your design? Is each pathway one recognisable shape, in the place you meant?
2. Is there room at the hubs, or do reactions and labels pile up?
3. Do lines run through areas they have nothing to do with, or along each other?
4. Only then the marks and numbers: crossings, crowded reactions, labels on lines, long lines.

Change the design numbers or the shape of a block in the script and build again. Prefer a change of the design
(more room, a block on the other side, a copy for a side pathway) to nudging single nodes. You may try variants: a
script argument that selects one, a build of each, and you keep the one that looks better.

The build's score and its last line are advice on the drawing, not the judge of the map: a map with a lower score
whose pathways sit where no biochemist would look for them is the worse map. Never meet a target with a trick (a copy
made only to satisfy a direction, a pathway moved away from where it belongs).

### E. Finish

When the map reads well: run `tidy`, build, look again (undo what you do not like by rerunning your script and
tidying only some blocks, `--block NAME`). Build the curvatures 0.6, 0.45 and 0.3 (`--curv`), choose by looking,
`validate` the chosen map. Stop when further changes no longer make the picture better, not at a number of builds.

## 5. Deliverables (in the work folder)

| file | content |
|---|---|
| the final map | the `*_handlayout.json` of your chosen build (your task description says where it has to be) |
| `make_spec.py`, the spec | the script and the spec it wrote (after `tidy`) |
| `PLAN.md`, `DECISIONS.md` | the design; every assumption and deliberate copy or overruled direction, one line each |
| `builds.md` | the table of all builds (written by `build`) |

## 6. Final message (at most 25 lines)

The design in three lines; what was merged and which copies are deliberate; the score of `build_0` and of the final
build with the targets still missed; the weak spots you see in the final picture; that the map was validated and
looked at but not loaded in Escher; the path of the final map.
