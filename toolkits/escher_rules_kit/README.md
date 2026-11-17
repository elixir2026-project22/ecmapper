# escher_rules_kit

A toolkit that lays out and judges Escher maps **by a rule set**: a human-readable YAML file of typed layout rules
(*when* a trigger holds, *do* an action; strength must / should / may / open). The rules are data, not code: change
the YAML and the kit's briefs, checks and verdicts change with it.

It builds on `escher_handlayout_kit` (the *engine*, in the folder beside this one): the engine turns a JSON layout
spec into a map; this kit says, from the rule set, what a layout of *this* map has to do and which rules a built
map breaks. Select this folder as the toolkit in ecmapper's AI layout panel, or use it from the command line.

| | escher_handlayout_kit | escher_rules_kit |
|---|---|---|
| guidance | placement catalogue measured from published maps, built-in score | a rule set you can read and edit (`rules/*.yaml`) |
| before the layout | `analysis.md`, `arrangement.md` | the same, plus `rules_brief.md`: the rules that apply to this map and where |
| after a build | score, targets, what to fix first | the same, plus the rules the map breaks, by strength, and a verdict |
| needs | numpy, matplotlib, Pillow | PyYAML; the engine for the layout commands |

## The rule set

```
rules/metabolic_map_rules.yaml          the rules (66), premises, conflicts, vocabulary, references
rules/metabolic_map_rules.checks.yaml   settings (currency metabolites, ...) and which probe looks at which rule
rules/metabolic_map_rules.md            the authors' text of the rule set (not read by the kit)
```

A rule:

```yaml
- id: M07
  type: motif
  group: 'Reaction: what one reaction looks like'
  when: Two antiparallel reactions between the same metabolites (PFK / FBPase)
  do: Draw both along one axis, not as a two-node cycle.
  strength: should
  basis: stated
  sources: [17]
  checked: ok
```

and its binding in the checks file:

```yaml
checks:
  M07: {probe: antiparallel_one_axis, max_offset: 350}
```

A *probe* finds the places of a map where a rule's trigger holds (PFK and FBP between f6p and fdp) and, where it
can, the places where the rule is broken. 25 of the 66 rules have a probe; 17 are marked as outside the layout of one
map (viewer, data overlay); the others are listed in every brief and report as rules to judge by eye. **docs/RULES.md**
says how to change rules, strengths, parameters and settings, and how to add a rule set or a probe.

## Commands

```bash
pip install -r requirements.txt
python selftest.py                                            # SELFTEST PASSED
python run.py rules --strength must                           # list rules (filters: --type --strength --group --id)
python run.py rules --validate                                # check the rule files after editing them
python run.py render --out rules.md                           # the rule set as Markdown tables
python run.py probes                                          # the probes, their parameters and their rules

MAP=../escher_handlayout_kit/examples/central_metabolism_escher.json
python run.py brief --input $MAP --out work/analysis          # which rules apply where (rules_brief.md)
python run.py check --input FINAL.json --draft $MAP --out work/check   # which rules a map breaks (rules_report.md)
python run.py check --input NEW.json --previous OLD.json      # also: no placed node has moved (A30)

# the layout, through the engine
python run.py inspect --input $MAP --out work/analysis        # engine analysis + rules_brief.md
python run.py start   --input $MAP --out work/spec.json
python run.py build   --input $MAP --spec work/spec.json --out work/build_0
```

`build` prints the engine's lines and then the verdict of the rules:

```text
rules: 13 checked on this map, 4 met | must: 2 broken (E10 8x, E30 5x) | should: 6 not met (M05 6x, ...) | 30 to judge by eye
  must E10 Currency metabolite: 8x, e.g. mqn8_c (node 672) is shared by NADH10 (reaction 685), NADH17pp (reaction 707): one copy per reaction
  ...
report: work/build_0/rules_report.md
RULES NOT MET: 2 must rules are broken (E10, E30). Repair these first; rules_report.md lists every place.
```

Every command takes `--rules FILE [FILE ...]`; without it all rule-set files of `rules/` are used together.
`brief`, `check`, `rules` and `render` take the filters `--id`, `--skip`, `--type`, `--strength`, `--group`.
The engine is found beside this kit, in a folder `engine` inside it, or through `ESCHER_HANDLAYOUT_KIT`.

## Who uses it how

| you | do this |
|---|---|
| maintain the rules | edit `rules/*.yaml`, run `python run.py rules --validate`, look at `render` |
| judge an existing map | `python run.py check --input MAP.json [--draft DRAFT.json]` |
| lay out with a chat model in ecmapper | choose this folder as the toolkit; the model gets `PROMPT.md` and the engine's `PROMPT.md` |
| lay out with a coding agent | `PROMPT_AGENT.md` (it builds on the engine's `PROMPT_AGENT.md`) |

## Contents

```
run.py                     command line (above)
PROMPT.md, PROMPT_AGENT.md what the rule set adds to the engine's procedure, for a chat model / a coding agent
ecmapper_toolkit.json      how ecmapper runs the kit
rules/                     the rule sets and their checks files
escher_rules/ruleset.py    reads and validates rule sets and checks files; renders Markdown
escher_rules/mapmodel.py   an Escher map as the probes see it: backbone, side participants, units, geometry
escher_rules/probes.py     the probes
escher_rules/report.py     brief and check: texts, JSON, verdict
docs/RULES.md              how to adapt the rule set
selftest.py, tests/        quick test; test suite (python tests/run_tests.py)
```

## Limitations

- A probe approximates a rule. A rule without violations is not proven met, and the lists *to judge by eye* are part
  of every report. Thresholds (how round is a ring, how far is adjacent) are parameters in the checks file, set by
  looking at one map (iML1515 central metabolism); expect to tune them.
- Units are the subsystems of the draft: the pieces of the draft map, named by its title labels. A laid-out map is
  checked with `--draft`; without it the map's own pieces are the units.
- What the engine cannot draw cannot be repaired by a layout: which nodes are primary (M01, M02; `scope: draft`),
  compartment enclosures (A20 to A22), clone markers (E12).
- Flux-guided orientation (M20 to M22) is not computed: a draft map carries no fluxes.
- Not loaded in Escher by the kit; the engine's limitations apply to the layout commands.
