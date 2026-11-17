"""Rule sets: human-readable YAML files of typed layout rules, and the companion files that bind rules to probes.

A rule-set file (rules/NAME.yaml) has `schema` (the allowed values of type, strength, basis, checked), `rules`
(id, type, group, when, do, strength, basis, sources, checked, note) and optionally `premises`, `conflicts`,
`vocabulary`, `references`. Its companion rules/NAME.checks.yaml has `settings` (which metabolites are currency, ...)
and `checks` (rule id -> probe and its parameters). See docs/RULES.md.
"""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RULES_DIR = os.path.join(os.path.dirname(HERE), 'rules')
CHECKS_SUFFIX = '.checks.yaml'
#: used when a rule-set file has no `schema` section
DEFAULT_SCHEMA = dict(type=['motif', 'arrangement', 'encoding'], strength=['must', 'should', 'may', 'open'],
                      basis=['stated', 'inferred', 'review'], checked={'ok': '', 'fix': '', '-': ''})
RULE_KEYS = {'id', 'type', 'group', 'when', 'do', 'strength', 'basis', 'sources', 'checked', 'note', 'check'}
REQUIRED = ('id', 'type', 'when', 'do', 'strength')
#: keys of a binding that are not parameters of its probe
BINDING_KEYS = {'probe', 'report_as', 'enabled', 'comment', 'scope'}
SCOPES = ('layout', 'draft', 'outside')


class RuleError(Exception):
    """A rule-set or checks file that cannot be used; the message lists every problem."""


def load_yaml(path):
    """The content of a YAML file (PyYAML, or ruamel.yaml where only that is installed)."""
    with open(path, encoding='utf-8') as f:
        text = f.read()
    try:
        import yaml
    except ImportError:
        yaml = None
    try:
        if yaml is not None:
            return yaml.safe_load(text)
        from ruamel.yaml import YAML
        return YAML(typ='safe', pure=True).load(text)
    except ImportError:
        raise RuleError('reading rule sets needs PyYAML: pip install pyyaml') from None
    except Exception as e:  # a syntax error of either library
        raise RuleError(f'{path} is not valid YAML: {e}') from None


class Rule:
    def __init__(self, data, source):
        self.data, self.source = data, source
        self.id = str(data.get('id', '?'))
        self.type = data.get('type', '')
        self.group = data.get('group', '')
        self.when = str(data.get('when', '')).strip()
        self.do = str(data.get('do', '')).strip()
        self.strength = data.get('strength', '')
        self.basis = data.get('basis', '')
        self.sources = list(data.get('sources') or [])
        self.checked = str(data.get('checked', '-'))
        self.note = str(data.get('note', '') or '').strip()
        self.binding = None  # dict(probe=..., parameters) once a checks file or the rule's own `check` binds it

    @property
    def probe(self):
        return self.binding.get('probe') if self.binding and self.binding.get('enabled', True) else None

    @property
    def params(self):
        return {k: v for k, v in (self.binding or {}).items() if k not in BINDING_KEYS}

    @property
    def report_as(self):
        """The strength a broken rule is reported with (a heuristic probe may report a `must` as `should`)."""
        return (self.binding or {}).get('report_as', self.strength)

    @property
    def scope(self):
        """`layout`; `draft` for a rule the layout cannot change (decided where the draft map is made); `outside`
        for a rule about something else than the layout of one map (the viewer, data overlays, other maps)."""
        return (self.binding or {}).get('scope', 'layout')

    def src(self):
        return ', '.join(f'[{s}]' for s in self.sources)


class RuleSet:
    """The rules of one or more rule-set files, with the settings and bindings of their checks files."""

    def __init__(self):
        self.files, self.rules, self.settings = [], [], {}
        self.schema = {k: (dict(v) if isinstance(v, dict) else list(v)) for k, v in DEFAULT_SCHEMA.items()}
        self.premises, self.conflicts, self.vocabulary, self.references = [], [], [], {}
        self.problems = []

    def __iter__(self):
        return iter(self.rules)

    def __len__(self):
        return len(self.rules)

    def get(self, rule_id):
        return next((r for r in self.rules if r.id == rule_id), None)

    def select(self, only=None, skip=None, types=None, strengths=None, groups=None):
        """The rules with these ids / types / strengths, or whose group contains one of `groups` (None: no filter)."""
        out = []
        for r in self.rules:
            if only and r.id not in only or skip and r.id in skip:
                continue
            if types and r.type not in types or strengths and r.strength not in strengths:
                continue
            if groups and not any(g.lower() in r.group.lower() for g in groups):
                continue
            out.append(r)
        return out

    # -- reading ---------------------------------------------------------------------------------------------------
    def add_file(self, path):
        data = load_yaml(path)
        name = os.path.basename(path)
        if not isinstance(data, dict) or not isinstance(data.get('rules'), list):
            raise RuleError(f"{path}: a rule-set file needs a list `rules`")
        self.files.append(path)
        for key, values in (data.get('schema') or {}).items():
            self.schema[key] = values
        allowed = {k: [str(x) for x in v] for k, v in self.schema.items()}
        refs = {str(k) for k in (data.get('references') or {})}
        for i, raw in enumerate(data['rules']):
            if not isinstance(raw, dict):
                self.problems.append(f'{name}: rule {i + 1} is not a mapping')
                continue
            rule = Rule(raw, name)
            where = f'{name}: rule {rule.id}'
            for key in REQUIRED:
                if not raw.get(key):
                    self.problems.append(f'{where}: `{key}` is missing')
            for key in sorted(set(raw) - RULE_KEYS):
                self.problems.append(f'{where}: unknown key `{key}` (known: {", ".join(sorted(RULE_KEYS))})')
            for key in ('type', 'strength', 'basis', 'checked'):
                if key in raw and str(raw[key]) not in allowed.get(key, []):
                    self.problems.append(f'{where}: {key} `{raw[key]}` is not one of {", ".join(allowed.get(key, []))}')
            if self.get(rule.id):
                self.problems.append(f'{where}: the id is used twice')
            if 'sources' in raw and not isinstance(raw['sources'], list):
                self.problems.append(f'{where}: `sources` must be a list of reference numbers')
            elif refs:
                for s in rule.sources:
                    if str(s) not in refs:
                        self.problems.append(f'{where}: source [{s}] is not in `references`')
            if isinstance(raw.get('check'), dict):
                rule.binding = dict(raw['check'])
            elif 'check' in raw:
                self.problems.append(f'{where}: `check` must be a mapping with `probe`')
            self.rules.append(rule)
        self.premises += data.get('premises') or []
        self.conflicts += data.get('conflicts') or []
        self.vocabulary += data.get('vocabulary') or []
        self.references.update(data.get('references') or {})
        known = {r.id for r in self.rules}
        for c in data.get('conflicts') or []:
            for rid in expand_ids(str(c.get('rules', '')), known):
                if rid not in known:
                    self.problems.append(f"{name}: conflict {c.get('id')}: rule {rid} does not exist")

    def add_checks(self, path):
        data = load_yaml(path) or {}
        name = os.path.basename(path)
        if not isinstance(data, dict):
            raise RuleError(f'{path}: a checks file has `settings` and `checks`')
        self.files.append(path)
        self.settings.update(data.get('settings') or {})
        for rule_id, binding in (data.get('checks') or {}).items():
            rule = self.get(str(rule_id))
            if rule is None:
                self.problems.append(f'{name}: checks: rule {rule_id} is not in the rule set')
            elif not isinstance(binding, dict) or not (binding.get('probe') or binding.get('scope')):
                self.problems.append(f'{name}: checks: {rule_id} needs `probe: NAME` (or `scope: outside`)')
            elif rule.binding is None:  # the rule's own `check` wins
                rule.binding = dict(binding)

    def check_bindings(self, probes):
        """Problems of the bindings against the registry of probes (unknown probe, unknown parameter, report_as)."""
        strengths = [str(x) for x in self.schema.get('strength', [])]
        for rule in self.rules:
            b = rule.binding
            if not b:
                continue
            if b.get('scope', 'layout') not in SCOPES:
                self.problems.append(f"rule {rule.id}: scope `{b['scope']}` is not one of {', '.join(SCOPES)}")
            if not b.get('probe'):
                continue
            probe = probes.get(b.get('probe'))
            if probe is None:
                self.problems.append(f"rule {rule.id}: unknown probe `{b.get('probe')}` (known: {', '.join(sorted(probes))})")
                continue
            for key in sorted(set(rule.params) - set(probe.defaults)):
                self.problems.append(f"rule {rule.id}: probe {probe.name} has no parameter `{key}` "
                                     f"(it has: {', '.join(sorted(probe.defaults)) or 'none'})")
            if 'report_as' in b and str(b['report_as']) not in strengths:
                self.problems.append(f"rule {rule.id}: report_as `{b['report_as']}` is not one of {', '.join(strengths)}")


def expand_ids(text, known):
    """The rule ids of a conflict's `rules` text: 'A04, E10-E14' -> A04 and the known ids from E10 to E14."""
    out = []
    for part in re.split(r'[,;\s]+', text.strip()):
        m = re.fullmatch(r'([A-Za-z]+)(\d+)-(?:\1)?(\d+)', part)
        if m:
            lo, hi = int(m.group(2)), int(m.group(3))
            out += sorted(k for k in known if re.fullmatch(rf'{m.group(1)}\d+', k) and lo <= int(k[len(m.group(1)):]) <= hi)
        elif part:
            out.append(part)
    return out


def default_files():
    """The rule-set files of the kit's rules folder (every *.yaml that is no checks file)."""
    return sorted(p for p in glob.glob(os.path.join(RULES_DIR, '*.yaml')) if not p.endswith(CHECKS_SUFFIX))


def load(paths=None, probes=None, strict=True):
    """The rule set of `paths` (default: all rule-set files of rules/), each with its companion NAME.checks.yaml if
    there is one. `probes`: the registry to check the bindings against. RuleError lists all problems (`strict`)."""
    paths = list(paths or default_files())
    if not paths:
        raise RuleError(f'no rule-set file given and none in {RULES_DIR}')
    rs = RuleSet()
    for path in paths:
        if not os.path.isfile(path):
            raise RuleError(f'{path}: no such rule-set file')
        rs.add_file(path)
    for path in paths:
        checks = re.sub(r'\.ya?ml$', '', path) + CHECKS_SUFFIX
        if os.path.isfile(checks):
            rs.add_checks(checks)
    if probes is not None:
        rs.check_bindings(probes)
    if strict and rs.problems:
        raise RuleError('the rule set has problems:\n  ' + '\n  '.join(rs.problems))
    return rs


# -- the rule set as Markdown ----------------------------------------------------------------------------------------
TYPE_HEADINGS = {'motif': 'Motif rules', 'arrangement': 'Arrangement rules', 'encoding': 'Encoding rules'}


def rule_row(rule):
    """One table row of a rule, as in the Markdown companion of a rule set."""
    cell = lambda s: str(s).replace('|', '\\|').replace('\n', ' ')
    return f'| {rule.id} | {cell(rule.when)} | {cell(rule.do)} | {rule.strength} | {rule.basis} | {rule.src()} | {rule.checked} |'


def to_markdown(rs, rules=None, title='Rule set', bindings=False):
    """The rule set as Markdown tables by type and group, with premises, conflicts, vocabulary and references."""
    rules = list(rs.rules if rules is None else rules)
    out = [f'# {title}', '', 'Written by `run.py render` from ' + ', '.join(f'`{os.path.basename(p)}`' for p in rs.files)
           + '. Change the YAML files, not this text.', '']
    types = list(dict.fromkeys([t for t in rs.schema.get('type', [])] + [r.type for r in rules]))
    for t in types:
        of_type = [r for r in rules if r.type == t]
        if not of_type:
            continue
        out += [f'## {TYPE_HEADINGS.get(t, str(t).capitalize() + " rules")}', '']
        for group in dict.fromkeys(r.group for r in of_type):
            in_group = [r for r in of_type if r.group == group]
            if group:
                out += [f'### {group}', '']
            head = '| ID | When | Do | Strength | Basis | Src | Chk |' + (' Probe |' if bindings else '')
            out += [head, '|---' * head.count('|', 1) + '|']
            out += [rule_row(r) + (f" {'`' + r.probe + '`' if r.probe else ''} |" if bindings else '') for r in in_group]
            notes = [f'- **{r.id}** {r.note}' for r in in_group if r.note]
            out += [''] + (notes + [''] if notes else [])
    if rs.premises and len(rules) == len(rs.rules):
        out += ['## Premises (not rules)', '', '| ID | Premise | Src | Chk |', '|---|---|---|---|']
        for p in rs.premises:
            text = str(p.get('text', '')) + (f" *{p['note']}*" if p.get('note') else '')
            out.append(f"| {p.get('id')} | {text} | {', '.join(f'[{s}]' for s in p.get('sources') or [])} | {p.get('checked', '-')} |")
        out.append('')
    if rs.conflicts and len(rules) == len(rs.rules):
        out += ['## Conflicts between sources', '', '| ID | Conflict | Rules | Suggested resolution |', '|---|---|---|---|']
        out += [f"| {c.get('id')} | {c.get('conflict')} | {c.get('rules')} | {c.get('suggested_resolution')} |" for c in rs.conflicts]
        out.append('')
    if rs.vocabulary and len(rules) == len(rs.rules):
        out += ['## Vocabulary', '', '| Term used here | Meaning | Called elsewhere |', '|---|---|---|']
        out += [f"| {v.get('term')} | {v.get('meaning')} | {v.get('synonyms')} |" for v in rs.vocabulary]
        out.append('')
    cited = {str(s) for r in rules for s in r.sources}
    refs = {k: v for k, v in rs.references.items() if len(rules) == len(rs.rules) or str(k) in cited}
    if refs:
        out += ['## References', ''] + [f'{k}. {v}' for k, v in refs.items()] + ['']
    return '\n'.join(out)
