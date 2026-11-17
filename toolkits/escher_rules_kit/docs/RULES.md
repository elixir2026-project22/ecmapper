# Adapting the rule set

Two files per rule set, both YAML, both meant to be edited by hand:

- `rules/NAME.yaml`: the rules themselves. Nothing in it refers to the kit.
- `rules/NAME.checks.yaml`: what the kit can find out about a map for these rules. Optional.

After every change run `python run.py rules --validate`: it names every problem (unknown strength, an id used
twice, a source that is not in the references, a probe or parameter that does not exist).

## 1. The rule-set file

```yaml
schema:                       # the allowed values; extend them here if you need another type or strength
  type: [motif, arrangement, encoding]
  strength: [must, should, may, open]
  basis: [stated, inferred, review]
  checked: {ok: confirmed against source text, fix: 'confirmed, wording corrected', '-': source not opened}
rules:
- id: M07                     # unique; the kit sorts nothing by it, the order of the file is the order of the texts
  type: motif
  group: 'Reaction: what one reaction looks like'      # heading in briefs and in `render`
  when: Two antiparallel reactions between the same metabolites (PFK / FBPase)
  do: Draw both along one axis, not as a two-node cycle.
  strength: should
  basis: stated
  sources: [17]               # numbers of `references`
  checked: ok
  note: optional remark
premises: [...]               # id, text, sources, checked, note
conflicts: [...]              # id, conflict, rules, suggested_resolution
vocabulary: [...]             # term, meaning, synonyms
references: {17: '...'}
```

Required per rule: `id`, `type`, `when`, `do`, `strength`. What the strength means to the kit:

| strength | in a brief | in a check |
|---|---|---|
| `must` | listed with its places | a broken must rule gives the verdict `RULES NOT MET` |
| `should` | listed with its places | reported; to be repaired or left with a written reason |
| `may` | listed with its places | reported as a note, never against the map |
| `open` | listed under *Decisions to take once* | not checked |

Common changes:

- **Make a rule stricter or looser:** change its `strength`.
- **Reword a rule:** change `when` / `do`; the texts go into briefs and reports as they are.
- **Add a rule:** append an entry with a new id. Without a binding it is a rule to judge by eye.
- **Drop a rule:** delete the entry (and its line in the checks file, which `--validate` asks for).
- **Decide an open rule:** change `strength: open` to `must` or `should` and write the decision into `do`.
- **A second rule set** (another organism group, a house style): put `OTHER.yaml` into `rules/`. All rule-set files
  of the folder are used together, so ids must differ; or name the files to use with `--rules`.

## 2. The checks file

```yaml
settings:
  currency: [atp, adp, nad, nadh, h2o, ...]     # base ids (without compartment) of the currency metabolites
  carriers: [coa, nad, nadh, atp, ...]          # must not carry the backbone if another pair can (M02)
  cofactor_pairs: [[atp, adp], [nad, nadh]]     # drawn on the same side of the backbone (M06)
  compartment_suffix: '_[a-z][a-z0-9]?$'        # what is cut off a metabolite id to get its base id
  exchange_prefixes: [EX_, DM_, SK_, sink_]     # exchange, sink and demand reactions
  biomass_pattern: '(?i)biomass'
  exit_label_prefix: '→'                        # text labels that are no unit titles
  primary_radius: 20                            # circle radii in map units
  secondary_radius: 12
checks:
  M07: {probe: antiparallel_one_axis, max_offset: 350}
  A05: {probe: boundary_reactions, edge_fraction: 0.25, report_as: should}
  M01: {probe: one_primary_pair, scope: draft}
  E02: {scope: outside}
```

A binding names a probe and may set its parameters (`python run.py probes` lists every probe with what it looks at,
its parameters and their defaults). Besides the parameters:

| key | meaning |
|---|---|
| `report_as: should` | report a broken rule with this strength instead of the rule's own; for a probe that is only a heuristic for its rule |
| `scope: draft` | the layout cannot change this (decided where the draft is made): reported apart, never against a layout |
| `scope: outside` | the rule is not about the layout of one map (viewer, data overlay, other maps): only named at the end; needs no probe |
| `enabled: false` | keep the binding, do not use it: the rule is judged by eye again |
| `comment: ...` | for the reader |

A rule may also carry its binding itself, as `check: {probe: ...}` in the rule-set file; it then wins over the
checks file.

Common changes:

- **A metabolite is currency in your models** (a quinone, ferredoxin, ACP): add its base id to `currency`. The probe
  `currency_candidates` (rule M04) lists the frequent metabolites of a map that are not in the list.
- **A threshold is too strict on your maps:** change the parameter in the binding (e.g. `min_roundness: 0.4`).
- **Your template says what sits where:** add relations to the binding of A15. Unit names are regular expressions,
  matched without regard to case; a relation whose units are not on a map is skipped.

  ```yaml
  A15:
    probe: relative_position
    relations:
      - {unit: 'glycolysis', is: above, other: 'citric acid cycle|TCA'}        # above, below, left_of, right_of
      - {unit: 'pentose phosphate', is: near, metabolite: g6p, max_distance: 1500}
  ```

  The same probe can serve A11 (biosynthesis left, degradation right) once you have decided conflict C4.
- **A rule should no longer be checked automatically:** `enabled: false`, or delete its binding.

## 3. What a probe returns, and adding one

A probe is a function in `escher_rules/probes.py`:

```python
@probe(max_offset=350)                      # the parameters and their defaults
def antiparallel_one_axis(model, max_offset):
    """One sentence that `run.py probes` prints."""
    ...
    return Finding(applies, violations)     # lists of texts; Finding(applies, checked=False) if it only finds places
```

`applies` are the places where the rule's trigger holds (shown in the brief), `violations` the places where the rule
is broken (shown in the check). Name reactions and nodes with their ids (`model.tag(node)`, `str(reaction)`): the
layout spec addresses them by these ids. `model` is a `MapModel` (`escher_rules/mapmodel.py`): `reactions`, `units`,
`main(reaction)` (the backbone nodes), `xy(node)`, `is_currency(node)`, and the geometry helpers of the module.
Bind the new probe to a rule in the checks file and add a test to `tests/run_tests.py`.

## 4. The Markdown text

`python run.py render --out rules.md` writes the rule set as tables (one per group, with notes, premises,
conflicts, vocabulary and references); `--bindings` adds the probe of every rule. `rules/metabolic_map_rules.md` is
the authors' own text with further sections (not yet covered, changes against the original list); the kit does not
read it, and the test suite checks that its rule rows agree with the YAML.
