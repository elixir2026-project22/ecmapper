"""The brief (which rules apply where on a map, before the layout) and the check (which rules a map breaks)."""
import json
import os

from . import probes as P
from .mapmodel import MapModel

SHOWN = 12  # places per rule in a Markdown text (the JSON files have all)


class Entry:
    """One rule looked at on one map."""

    def __init__(self, rule, finding):
        self.rule, self.finding = rule, finding

    @property
    def state(self):
        """no_probe | found (the probe only finds places) | not_here (trigger not on this map) | met | broken |
        not_checked (the probe lacks an input, e.g. the earlier map)."""
        f = self.finding
        if f is None:
            return 'no_probe'
        if not f.checked:
            return 'found' if f.applies else ('not_checked' if f.summary else 'not_here')
        if f.violations:
            return 'broken'
        return 'met' if f.applies or f.summary else 'not_here'

    def as_dict(self):
        r, f = self.rule, self.finding
        d = dict(id=r.id, type=r.type, group=r.group, when=r.when, do=r.do, strength=r.strength, report_as=r.report_as,
                 scope=r.scope, probe=r.probe, state=self.state)
        if f is not None:
            d.update(applies=f.applies, violations=f.violations, summary=f.summary)
        return d


def evaluate(ruleset, model, previous=None, rules=None):
    """An Entry for every rule (of `rules`, default all) on `model`."""
    out = []
    for rule in (ruleset.rules if rules is None else rules):
        finding = P.run(rule.probe, model, rule.params, previous) if rule.probe else None
        out.append(Entry(rule, finding))
    return out


def _head(rule):
    return f'{rule.id} ({rule.strength}) {rule.when}'


def _places(lines, shown=SHOWN):
    out = [f'- {x}' for x in lines[:shown]]
    if len(lines) > shown:
        out.append(f'- ... and {len(lines) - shown} more (all in the .json file beside this text)')
    return out


def _table(entries):
    out = ['| ID | When | Do | Strength |', '|---|---|---|---|']
    cell = lambda s: str(s).replace('|', '\\|')
    return out + [f'| {e.rule.id} | {cell(e.rule.when)} | {cell(e.rule.do)} | {e.rule.strength} |' for e in entries]


def _write(out_dir, name, text, data):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, name + '.md'), 'w', encoding='utf-8') as f:
        f.write(text)
    with open(os.path.join(out_dir, name + '.json'), 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    return os.path.join(out_dir, name + '.md')


def _units_line(model):
    sizes = sorted(((len(rs), u) for u, rs in model.units.items()), reverse=True)
    return f"{len(model.reactions)} reactions in {len(sizes)} units: " + ', '.join(f'{u} ({n})' for n, u in sizes)


def _outside(entries):
    if not entries:
        return []
    return ['', '## Outside the layout of this map', '',
            'Rules about the viewer, data overlays, other maps or how the model is cut into units (`scope: outside` in '
            'the checks file); nothing to do for them here: ' + ', '.join(f'{e.rule.id} ({e.rule.when})' for e in entries)]


# -- the brief ---------------------------------------------------------------------------------------------------------
def brief(ruleset, map_path, out_dir=None, rules=None):
    """What the rule set asks of this map, before it is laid out: (text, data); written to out_dir/rules_brief.md
    and .json if `out_dir` is given."""
    model = MapModel(map_path, ruleset.settings)
    entries = evaluate(ruleset, model, rules=rules)
    here = [e for e in entries if e.finding is not None and (e.finding.applies or e.finding.violations)
            and e.rule.scope != 'draft' and e.rule.strength != 'open']
    given = [e for e in entries if e.rule.scope == 'draft' and e.finding is not None
             and (e.finding.violations or e.finding.applies)]
    undecided = [e for e in entries if e.rule.strength == 'open']
    by_eye = [e for e in entries if e.state == 'no_probe' and e.rule.strength != 'open' and e.rule.scope == 'layout']
    absent = [e for e in entries if e.state == 'not_here' and e.rule.strength != 'open']
    outside = [e for e in entries if e.rule.scope == 'outside']
    L = [f'# Rules for this map: {model.header.get("map_name") or os.path.basename(map_path)}', '',
         f'{_units_line(model)}. A unit is a subsystem of the draft (its title label).', '',
         f'Rule set: {", ".join(os.path.basename(p) for p in ruleset.files)} ({len(entries)} rules). '
         'Strength: **must** = breaking it makes the map wrong or unreadable; **should** = the default convention, '
         'leave it only for a reason you write down; **may** = an option; **open** = undecided, decide once.', '',
         '## 1. Rules with places found on this map', '',
         'For each rule: what to do, and where on this map its trigger holds. Work through these when you design.', '']
    for e in here:
        L += [f'### {_head(e.rule)}', '', f'Do: {e.rule.do}', '']
        if e.rule.note:
            L += [f'Note: {e.rule.note}', '']
        if e.finding.applies:
            L += ['On this map:'] + _places(e.finding.applies, 8) + ['']
        if e.finding.checked and e.finding.violations:
            L += [f'In the draft as it is, {len(e.finding.violations)} place(s) break this rule; the first:']
            L += _places(e.finding.violations, 3) + ['']
    L += ['## 2. Decisions to take once, before the layout', '']
    for e in undecided:
        L += [f'- **{e.rule.id}** {e.rule.when}: {e.rule.do}']
        if e.finding is not None and e.finding.applies:
            L += ['  On this map:'] + ['  ' + x for x in _places(e.finding.applies, 8)]
    for c in ruleset.conflicts:
        L += [f"- **{c.get('id')}** (rules {c.get('rules')}) {c.get('conflict')}. Suggested: {c.get('suggested_resolution')}"]
    if not undecided and not ruleset.conflicts:
        L += ['- none']
    L += ['', 'Write each decision as one line into DECISIONS.md.', '',
          '## 3. Given by the draft', '',
          'The layout cannot change these (which nodes are primary, which reactions a unit has). Know them; do not try to repair them.', '']
    for e in given:
        L += [f'- **{e.rule.id}** {e.rule.when}: {e.rule.do}']
        lines = e.finding.violations or e.finding.applies
        L += ['  ' + x for x in _places(lines, 6)]
    if not given:
        L += ['- nothing']
    L += ['', '## 4. Rules the kit cannot place: keep them in mind and judge by eye', ''] + _table(by_eye)
    if absent:
        L += ['', '## 5. Rules whose trigger the kit did not find on this map', '',
              ', '.join(f'{e.rule.id} ({e.rule.when})' for e in absent)]
    L += _outside(outside)
    text = '\n'.join(L) + '\n'
    data = dict(map=map_path, rule_files=ruleset.files, units={u: [r.bigg for r in rs] for u, rs in model.units.items()},
                rules=[e.as_dict() for e in entries], conflicts=ruleset.conflicts)
    if out_dir:
        _write(out_dir, 'rules_brief', text, data)
    return text, data


# -- the check ---------------------------------------------------------------------------------------------------------
def check(ruleset, map_path, draft=None, previous=None, out_dir=None, rules=None):
    """The rules `map_path` breaks: (text, data, summary lines); written to out_dir/rules_report.md and .json if
    `out_dir` is given. `draft`: the draft the map was laid out from (its subsystems are the units); `previous`: an
    earlier layout whose positions must be kept."""
    model = MapModel(map_path, ruleset.settings, units_from=draft)
    before = MapModel(previous, ruleset.settings) if previous else None
    entries = evaluate(ruleset, model, before, rules)
    strengths = [s for s in ruleset.schema.get('strength', []) if s != 'open']
    layout = [e for e in entries if e.rule.scope == 'layout']
    broken = {s: [e for e in layout if e.state == 'broken' and e.rule.report_as == s] for s in strengths}
    met = [e for e in layout if e.state == 'met']
    given = [e for e in entries if e.rule.scope == 'draft' and e.state == 'broken']
    by_eye = [e for e in entries if e.state in ('no_probe', 'found') and e.rule.scope != 'outside']
    outside = [e for e in entries if e.rule.scope == 'outside' and e.state in ('no_probe', 'found', 'not_here')]
    rest = [e for e in entries if e.state in ('not_here', 'not_checked')]
    words = {'must': 'broken', 'should': 'not met', 'may': 'noted'}
    name = model.header.get('map_name') or os.path.basename(map_path)
    L = [f'# Rule check: {name}', '', f'`{map_path}`: {_units_line(model)}.', '']
    counts = lambda es: ', '.join(f'{e.rule.id} {len(e.finding.violations)}x' for e in es)
    n_checked = len(met) + sum(len(v) for v in broken.values())
    head = [f'rules: {n_checked} checked on this map, {len(met)} met']
    head += [f'{s}: {len(es)} {words.get(s, "not met")} ({counts(es)})' for s, es in broken.items() if es]
    head += [f'{len(by_eye)} to judge by eye']
    summary = [' | '.join(head)]
    for s, es in broken.items():
        if not es:
            continue
        L += [f'## {s.capitalize()} rules {words.get(s, "not met")}', '']
        for e in es:
            own = '' if e.rule.report_as == e.rule.strength else f' (a {e.rule.strength} rule, reported as {s})'
            L += [f'### {_head(e.rule)}{own}', '', f'Do: {e.rule.do}', ''] + _places(e.finding.violations) + ['']
            if s != 'may':
                first = e.finding.violations[0]
                summary.append(f'  {s} {e.rule.id} {e.rule.when}: {len(e.finding.violations)}x, e.g. {first[:200]}')
    if given:
        L += ['## Given by the draft (not counted: the layout cannot change these)', '']
        for e in given:
            L += [f'- **{e.rule.id}** {e.rule.when}: {len(e.finding.violations)} place(s), e.g. {e.finding.violations[0]}']
        L += ['']
    L += ['## Rules met', '']
    L += [f'- **{e.rule.id}** {e.rule.when}' + (f': {e.finding.summary}' if e.finding.summary else '') for e in met] or ['- none']
    L += ['', '## To judge by eye', '',
          'No probe can decide these. Look at the pictures of the map and go through them; where the kit found places, '
          'their number is given (rules_brief.md lists them).', '', '| ID | When | Do | Strength | Places |', '|---|---|---|---|---|']
    cell = lambda s: str(s).replace('|', '\\|')
    L += [f"| {e.rule.id} | {cell(e.rule.when)} | {cell(e.rule.do)} | {e.rule.strength} | "
          f"{len(e.finding.applies) if e.finding is not None else ''} |" for e in by_eye]
    if rest:
        L += ['', '## Not on this map, or not checked', '']
        L += [f'- **{e.rule.id}** {e.rule.when}' + (f': {e.finding.summary}' if e.finding is not None and e.finding.summary else '')
              for e in rest]
    L += _outside(outside)
    must, should = broken.get('must', []), broken.get('should', [])
    if must:
        verdict = (f"RULES NOT MET: {len(must)} must rule{'s are' if len(must) > 1 else ' is'} broken "
                   f"({', '.join(e.rule.id for e in must)}). Repair these first; rules_report.md lists every place.")
    elif should:
        verdict = (f"MUST RULES MET; {len(should)} should rule{'s' if len(should) > 1 else ''} not met "
                   f"({', '.join(e.rule.id for e in should)}). Repair those you can without spoiling the design; "
                   'for each one you leave, write one line with the reason into DECISIONS.md.')
    else:
        verdict = 'ALL CHECKED RULES MET. Go through "To judge by eye" in rules_report.md with the pictures of the map.'
    state = 'not_met' if must else ('should_open' if should else 'met')
    L += ['', '## Verdict', '', verdict]
    text = '\n'.join(L) + '\n'
    data = dict(map=map_path, draft=draft, previous=previous, rule_files=ruleset.files, state=state, verdict=verdict,
                broken={s: [e.rule.id for e in es] for s, es in broken.items()},
                violations=sum(len(e.finding.violations) for es in broken.values() for e in es),
                rules=[e.as_dict() for e in entries])
    if out_dir:
        summary.append('report: ' + _write(out_dir, 'rules_report', text, data))
    summary.append(verdict)
    return text, data, summary
