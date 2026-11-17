# Lay out an Escher draft map by a rule set (coding agent)

You are a coding agent. This kit, **escher_rules_kit**, lays out an Escher draft map with **escher_handlayout_kit**
(the engine) and judges it by a **rule set** (`rules/*.yaml`). Your procedure is the engine's `PROMPT_AGENT.md`
(the second instruction file; `RUN doc prompt_agent` prints it) with the additions below. `RUN` is the command your
task description gives for this kit; every engine command is run through it (`RUN inspect ...`, `RUN build ...`).

## The rule set

`rules/NAME.yaml` holds typed rules: *when* a trigger holds, *do* an action, with a strength (`must`: breaking it
makes the map wrong or unreadable; `should`: the default convention; `may`: an option; `open`: undecided).
`rules/NAME.checks.yaml` binds rules to probes of the kit. `RUN rules` lists the rules, `RUN rules --id M07` one.

```text
RUN inspect --input MAP --out analysis      also writes analysis/rules_brief.md: the rules that apply to this map,
                                            and the places where (reactions, metabolites, units with their ids)
RUN build ...                               also checks the built map; the lines `rules: ...` and the last line are
                                            the verdict of the rule set; build_N/rules_report.md lists every place
RUN check --input MAP.json --draft DRAFT.json [--out DIR] [--id E30 E31]    the same check for any map
RUN doc [NAME]                              a document of the engine (format_reference, guide, prompt_agent, ...)
```

## What changes in the engine's procedure

- **A (understand).** Read `analysis/rules_brief.md` completely. Section 2 lists decisions to take once (e.g. whether
  shared metabolites are duplicated at unit boundaries): take them now, one line each in `DECISIONS.md`.
- **B (design).** `PLAN.md` gets one more table: for every rule of section 1 of the brief with places on this map,
  what your design does about it (a ring for each listed cycle, one axis for each antiparallel pair, which units are
  neighbours and why, which metabolites stay copies). The rules to judge by eye (section 4) are design constraints
  as well: read them before the sketch. Where a rule contradicts the engine's placement catalogue
  (`analysis/arrangement.md`), the rule wins if it is `must` or `should`.
- **D (look, then change).** Judge in the engine's order, with one step first: **a broken must rule is repaired
  before anything else**. Then the picture, then the should rules, then the engine's marks and numbers. A should rule
  you leave unmet needs one line of reason in `DECISIONS.md`.
- **E (finish).** Go through *To judge by eye* in the last `rules_report.md` with the overview and the views of the
  hubs. The final message also names: the must rules (all met, or which are not and why), the should rules left
  unmet with their reasons, and the decisions of section 2.

## Limits

- Never change `rules/*.yaml` or `*.checks.yaml`, and never meet a rule with a trick (a copy made only to pass a
  check, a pathway moved away from where it belongs). The rules describe a good map; the checks only approximate
  them.
- Rules reported under *Given by the draft* cannot be changed by a layout.
- The report of the rules is advice about conventions, as the engine's score is advice about the drawing: you look
  at the pictures and decide, and you say what you decided.
