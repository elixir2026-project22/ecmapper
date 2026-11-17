# Lay out the map by the rule set

This toolkit lays out an Escher draft map with **escher_handlayout_kit** (the engine: spec, build, pictures) and
judges the result by a **rule set** (`rules/*.yaml`: typed rules of the form *when* a trigger holds, *do* an action).
The engine's procedure follows this text and stays valid. This text adds what the rules change. Where the two
disagree, a **must** rule wins; otherwise the engine's procedure wins.

## What the rules add to the procedure

1. **After `inspect`** read `analysis/rules_brief.md` (written by `inspect`): section 1 lists the rules that apply to
   this map and the places where, section 2 the decisions to take once, section 3 what the draft fixes. Write the
   decisions of section 2 into `DECISIONS.md`, one line each, before you design.
2. **When you design and change the spec**, use the brief as the list of things to do. How a rule becomes a change of
   the spec:

   | rule says | change of the spec |
   |---|---|
   | a currency metabolite is shared by reactions (E10) | `splits`: one copy of the node per reaction |
   | a metabolite is drawn several times (E11) | `merges`: keep one node, the others become `dups` (not for currency metabolites) |
   | a cycle is drawn as a ring (M12), two cycles share a metabolite (M15) | `rings` in the block, the ring reactions with `arc`; two rings touch at the shared node |
   | antiparallel reactions on one axis (M07) | both reactions between the same two nodes, the second with a small `shift` |
   | paired cofactors on the same side (M06) | the reaction's `side` |
   | side participants short, backbone straight (M05) | the reaction's `dist`, `bend`; move the primary nodes into a line |
   | flow top to bottom (A12), straight axis-aligned runs (M23) | node coordinates of the block: a column, equal x |
   | which units are neighbours, what is above what (A15, A17) | the blocks' `origin` |
   | nodes overlap, a line crosses a node (E30, E31) | `tidy`; then move the node or change `side` / `dist` |

3. **After every `build`** two verdicts are printed. The engine's line (`NOT FINISHED`, `TARGETS MET`, `STOP`, ...)
   is about the drawing. The lines that start with `rules:` and the **last line** are about the rule set:
   - `RULES NOT MET`: a must rule is broken. Repair it before anything else; `build_N/rules_report.md` lists every
     place with its node and reaction ids.
   - `MUST RULES MET; ... should rules not met`: go on with the engine's `Fix first`. Repair the should rules you can
     repair without spoiling the design; for each one you leave, write one line with the reason into `DECISIONS.md`.
   - `ALL CHECKED RULES MET`: only the engine's line counts.
4. **Before you hand in the map** read the section *To judge by eye* of the last `rules_report.md` and look at the
   overview picture with it. These rules have no automatic check; change what clearly breaks one.

A map is finished when the last line of its build does not say `RULES NOT MET` and the engine's line says so
(`TARGETS MET`, `STOP`, or only targets left that you overrule for a reason).

## Limits

- The rule files are the user's: never change `rules/*.yaml`. If a rule cannot be met on this map, say so in
  `DECISIONS.md`; do not look for a trick that satisfies the check.
- Rules reported under *Given by the draft* (which node is primary, the size of a unit) cannot be changed by a
  layout: do not try.
- Read a document of the engine with the command `doc`: `doc format_reference` (every key of a spec), `doc guide`
  (recipes; long, in parts: `doc guide --part 2`).
- `rules --id E10 M12` prints single rules; `check --input MAP --draft draft.escher.json` checks any map again.
