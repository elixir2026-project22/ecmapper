"""Where the subsystems of a map usually sit relative to each other, from the placement catalogue in ../layout_directions
(see layout_directions/README_v2.md; the same catalogue drives pwl_handlayout_kit/pwl_kit/directions.py).

The subsystems of a map (its title labels / zones) are matched to the catalogue's entities by `brenda_aliases.json`, by
the entities' names and aliases, or by word overlap; the map's primary metabolites are matched to the catalogue's hub
metabolites by their BiGG ids. :class:`Arrangement` then collects for the map

- directions between subsystems: stated rules (`direction_rules`) and the consensus of the measured maps (`consensus`,
  at the level of the matched entities or, if there is none, of their pathways); conflicting maps give no direction;
- directions between hub metabolites on the map (hard and strong textbook rules, e.g. glucose above G6P);
- attachments: relations whose `via` hubs are on the map ("citric acid cycle feeds glutamate metabolism at akg_c") and
  the entry points of cycles - these are merge candidates;
- cycle layouts (order, orientation, first member) of the cycles on the map;
- for subsystems without any of these: their usual region on overview maps (`absolute_placements`) and role.

:func:`check` measures a layout (zone centres and node positions, from a built Layout or from a spec) against the
directions: the violated hard, strong and medium ones and a penalty for the score. Screen coordinates throughout:
x grows to the right, y grows DOWN, so "above" means a smaller y. The catalogue's `direction` always gives the
position of the TARGET relative to the ANCHOR.
"""
import json, math, os, re
from collections import Counter, defaultdict

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'layout_directions')
CATALOGUE = os.path.join(DATA, 'placement_directions_v2.json')
ALIASES = os.path.join(DATA, 'brenda_aliases.json')
COMPASS = ('above', 'above_right', 'right', 'below_right', 'below', 'below_left', 'left', 'above_left')
STRENGTH = {'hard': 4, 'strong': 3, 'medium': 2, 'weak': 1, 'none': 0}
#: penalty of a violated direction in the layout score (by strength); weak directions are hints only
PENALTY = {'hard': 5, 'strong': 5, 'medium': 2}
ATTACHING = ('FEEDS', 'SHARES_HUB', 'ALTERNATIVE_ROUTE', 'REVERSE_OF')
VERB = {'FEEDS': 'feeds', 'SHARES_HUB': 'shares a hub with', 'ALTERNATIVE_ROUTE': 'is an alternative route to', 'REVERSE_OF': 'is the reverse of'}
ROLE_SIDE = {'biosynthesis': 'west (left)', 'degradation': 'east (right)', 'energy': 'middle'}
#: words of map titles -> words of the catalogue (as in pwl_handlayout_kit/tools/make_aliases.py)
WORDS = {'anapleurotic': 'anaplerotic', 'oxalacetate': 'oxaloacetate', 'catabolism': 'degradation', 'synthesis': 'biosynthesis',
         'ketoglutarate': 'oxoglutarate', 'citric': 'tca', 'krebs': 'tca', 'acid': '', 'glycolytic': 'glycolysis', 'isoprenoid': 'terpenoid'}
STOP = {'of', 'the', 'and', 'a', 'an', 'in', 'by', 'to', 'from', 'via', 'with', 'i', 'ii', 'iii', 'pathway', 'metabolism'}
ROLE = {'biosynthesis', 'degradation', 'fermentation', 'cycle'}


def slug(s):
    return re.sub(r'[^a-z0-9]+', '-', str(s).lower()).strip('-')

def words(text):
    out = set()
    for w in re.findall(r'[a-z0-9]+', str(text).lower().replace('_', ' ')):
        w = WORDS.get(w, w)
        if w and w not in STOP:
            plural = len(w) > 4 and w.endswith('s') and not w.endswith(('ss', 'is', 'us'))
            out.add(w[:-1] if plural else w)
    return out

def base_bigg(bigg):
    """BiGG id without the compartment: akg_c -> akg."""
    return re.sub(r'_[a-z]{1,2}$', '', str(bigg))

def compass(direction):
    d = direction[len('adjacent_'):] if direction.startswith('adjacent_') else direction
    return d if d in COMPASS else None

def _taxa_ok(scope, taxa):
    t = scope.get('taxa', ['all'])
    return taxa == 'all' or 'all' in t or taxa in t


class Catalogue:
    """The placement catalogue (loaded once per process). `Catalogue.load()` is None when layout_directions/ is missing."""
    _cache = {}

    def __init__(self, path=CATALOGUE, aliases=ALIASES):
        with open(path, encoding='utf-8') as f: self.data = json.load(f)
        self.path = path
        self.entities = {e['id']: e for e in self.data['entities']}
        self.aliases = {}
        if os.path.isfile(aliases):
            with open(aliases, encoding='utf-8') as f: self.aliases = json.load(f).get('subsystems', {})
        self.hub_of_bigg = {b: e['id'] for e in self.data['entities'] if e['type'] == 'metabolite_hub' for b in e['xrefs'].get('bigg_metabolite', [])}
        # names for matching: entity id -> list of word sets (name, id, aliases without their "SOURCE: " prefix)
        self._names = {}
        for e in self.data['entities']:
            if e['type'] == 'metabolite_hub': continue
            strings = [e['name'], e['id']] + [a.split(':', 1)[-1] for a in e.get('aliases', [])]
            self._names[e['id']] = [(s.strip(), words(s)) for s in strings]
        self._exact = {}
        for eid, lst in self._names.items():
            for s, _ in lst: self._exact.setdefault(' '.join(sorted(words(s))), eid)

    @classmethod
    def load(cls, path=CATALOGUE, aliases=ALIASES):
        key = (path, aliases)
        if key not in cls._cache:
            cls._cache[key] = cls(path, aliases) if os.path.isfile(path) else None
        return cls._cache[key]

    @property
    def version(self): return self.data.get('version', '?')

    def pathway(self, entity):
        """The pathway an entity (module, hub) belongs to; itself for pathways, categories, regions and unknown ones."""
        seen, todo = set(), [entity]
        while todo:
            e = todo.pop(0)
            if e in seen or e not in self.entities: continue
            seen.add(e)
            if self.entities[e]['type'] in ('pathway', 'category', 'region'): return e
            todo += self.entities[e].get('part_of', [])
        return entity

    def match(self, title):
        """(entity ids, how) for a subsystem title: 'alias' (brenda_aliases.json), 'name' (an entity name or alias),
        'words' (word overlap >= 0.5, no biosynthesis/degradation mix-up) or ([], None)."""
        a = self.aliases.get(slug(title))
        if a and a.get('entities'): return list(a['entities']), 'alias'
        w = words(title)
        if not w: return [], None
        key = ' '.join(sorted(w))
        if key in self._exact: return [self._exact[key]], 'name'
        best, best_s = None, 0.0
        wr = w & ROLE
        for eid, lst in self._names.items():
            for s, ws in lst:
                if not ws: continue
                roles = ws & ROLE
                if wr and roles and wr != roles and (wr | roles) <= {'biosynthesis', 'degradation'}: continue
                j = len((w - ROLE) & (ws - ROLE)) / max(1, len((w - ROLE) | (ws - ROLE)))
                if j > best_s: best, best_s = eid, j
        return ([best], 'words') if best and best_s >= 0.5 else ([], None)


class Direction:
    __slots__ = ('anchor', 'target', 'direction', 'strength', 'evidence')

    def __init__(self, anchor, target, direction, strength, evidence):
        self.anchor, self.target, self.direction, self.strength, self.evidence = anchor, target, direction, strength, evidence

    def where(self):
        """E.g. "below glycolysis", "above right of glycolysis"."""
        d = self.direction.replace('_', ' ')
        return f'{d} {self.anchor}' if d in ('above', 'below') else f'{d} of {self.anchor}'

    def line(self): return f'| {self.target} | {self.where()} | {self.strength} | {self.evidence} |'

    def to_dict(self): return dict(anchor=self.anchor, target=self.target, direction=self.direction, strength=self.strength, evidence=self.evidence)


class Arrangement:
    """The placement brief of a map. `subsystems`: ordered dict  key -> dict(title=<text of the title label>,
    metabolites=[BiGG ids of the primary metabolites]); the keys are the names used in the output (zone names)."""

    def __init__(self, subsystems, catalogue=None, taxa='all'):
        self.cat = catalogue or Catalogue.load()
        if self.cat is None: raise FileNotFoundError(f'placement catalogue not found: {CATALOGUE}')
        self.taxa = taxa; self.keys = list(subsystems)
        self.titles = {k: v.get('title') or k for k, v in subsystems.items()}
        self.entities, self.how = {}, {}
        for k in self.keys: self.entities[k], self.how[k] = self.cat.match(self.titles[k])
        self.hubs = defaultdict(set)              # hub entity -> BiGG ids (with compartment) on the map
        for k, v in subsystems.items():
            for b in v.get('metabolites', []):
                h = self.cat.hub_of_bigg.get(base_bigg(b))
                if h: self.hubs[h].add(b)
        self.subsystem_directions = self._subsystem_directions()
        self.metabolite_directions = self._metabolite_directions()

    # ------------------------------------------------------------------ directions
    def _owner(self, entity, level):
        if level == 'exact': return [k for k in self.keys if entity in self.entities[k]]
        return [k for k in self.keys if any(self.cat.pathway(e) == entity for e in self.entities[k])]

    def _subsystem_directions(self):
        best = {}
        def offer(a, b, direction, strength, evidence):
            d = compass(direction)
            if a == b or d is None or STRENGTH.get(strength, 0) == 0: return
            key = frozenset((a, b)); old = best.get(key)
            if old is None or STRENGTH[strength] > STRENGTH[old.strength]: best[key] = Direction(a, b, d, strength, evidence)
        for r in self.cat.data['direction_rules']:
            if not _taxa_ok(r.get('scope', {}), self.taxa): continue
            d = r['direction']
            a_ids, b_ids = r['anchor'], r['target']
            if r.get('flow') == 'target_to_anchor': a_ids, b_ids = b_ids, a_ids
            for a in self._owner(a_ids, 'exact'):
                for b in self._owner(b_ids, 'exact'):
                    offer(a, b, d, r['strength'], (f"rule {r['id']}: {r.get('note', '')[:70]}").rstrip(': '))
        exact_pairs = set()
        for level in ('exact', 'pathway'):
            for c in self.cat.data['consensus']:
                # one map is no evidence for a convention; conflicting maps give no direction
                if c['level'] != level or c['status'] in ('conflict', 'single_source') or c['n_sources'] < 2: continue
                for a in self._owner(c['anchor'], level):
                    for b in self._owner(c['target'], level):
                        if level == 'pathway' and frozenset((a, b)) in exact_pairs: continue
                        if level == 'exact': exact_pairs.add(frozenset((a, b)))
                        ev = f"{c['status']}, {c['n_sources']} sources (R {c['agreement_R']:.2f})" + (' for their pathways' if level == 'pathway' else '')
                        offer(a, b, c['direction'], c['suggested_strength'], ev)
        return sorted(best.values(), key=lambda d: (-STRENGTH[d.strength], d.anchor, d.target))

    def _metabolite_directions(self):
        out = {}
        for r in self.cat.data['direction_rules']:
            d = compass(r['direction'])
            if d is None or r['strength'] not in ('hard', 'strong') or not _taxa_ok(r.get('scope', {}), self.taxa): continue
            for a in sorted(self.hubs.get(r['anchor'], ())):
                for b in sorted(self.hubs.get(r['target'], ())):
                    if (a, b) not in out or STRENGTH[r['strength']] > STRENGTH[out[a, b].strength]:
                        out[a, b] = Direction(a, b, d, r['strength'], r.get('note', '')[:70])
        return list(out.values())

    # ------------------------------------------------------------------ attachments, cycles, fallback
    def attachments(self):
        lines = []
        for r in self.cat.data['relations']:
            if r['relationship'] not in ATTACHING or not _taxa_ok(r.get('scope', {}), self.taxa): continue
            for a in self._owner(r['anchor'], 'exact'):
                for b in self._owner(r['target'], 'exact'):
                    if a == b: continue
                    via = sorted(m for h in r.get('via', []) for m in self.hubs.get(h, ()))
                    at = f" at {', '.join(via)}" if via else ''
                    note = f" ({r['note']})" if r.get('note') else ''
                    lines.append(f"- {a} {VERB[r['relationship']]} {b}{at}{note}")
        return sorted(set(lines))

    def cycles(self):
        lines = []
        on_map = {e for es in self.entities.values() for e in es}
        for c in self.cat.data['cycle_layouts']:
            if c['cycle'] not in on_map: continue
            members = [sorted(self.hubs.get(h, {h}))[0] for h in c['order']]
            k = next(k for k in self.keys if c['cycle'] in self.entities[k])
            first = f", first member at the {c['first_member_position'].replace('_', ' ')}" if c.get('first_member_position') else ''
            lines.append(f"- {k}: {' -> '.join(members)}, {c['orientation'].replace('_', ' ')}{first} ({c['strength']})")
            for entity, member in sorted(c.get('entry_points', {}).items()):
                where = sorted(self.hubs.get(member, {member}))[0]
                if entity in self.hubs: lines.append(f"  - {sorted(self.hubs[entity])[0]} enters at {where}")
                for other in self._owner(entity, 'exact'):
                    if other != k: lines.append(f"  - {other} attaches at {where}")
        return lines

    def fallback(self):
        """Usual region and role of the subsystems that have no direction of at least medium strength."""
        placed = {x for d in self.subsystem_directions if STRENGTH[d.strength] >= 2 for x in (d.anchor, d.target)}
        lines = []
        for k in self.keys:
            if k in placed or not self.entities[k]: continue
            ids = set(self.entities[k]) | {self.cat.pathway(e) for e in self.entities[k]}
            regions = Counter(a['region'] for a in self.cat.data['absolute_placements'] if a['entity'] in ids)
            roles = [self.cat.entities[e].get('process_role') for e in self.entities[k] if e in self.cat.entities]
            parts = []
            if regions: parts.append('usually ' + ', '.join(f'{r} ({n})' for r, n in regions.most_common(2)) + ' on overview maps')
            role = next((r for r in roles if r), None)
            if role in ROLE_SIDE: parts.append(f'{role}: {ROLE_SIDE[role]} (MetaCyc convention)')
            if parts: lines.append(f"- {k}: {'; '.join(parts)}")
        return lines

    # ------------------------------------------------------------------ output
    def matches(self):
        """Lines "key: entities (how)" for every subsystem."""
        out = []
        for k in self.keys:
            e = ', '.join(self.entities[k]) or '(no catalogue entity)'
            how = {'alias': '', 'name': '', 'words': ' (matched by word overlap)', None: ''}[self.how[k]]
            out.append(f'- {k}: {e}{how}')
        return out

    def text(self, max_chars=None, intro=True):
        """Markdown brief. `max_chars` keeps the most important parts (directions first)."""
        out = []
        if intro:
            out += ['# Arrangement of the subsystems (placement catalogue, layout_directions/)', '',
                    'Directions: where the subsystem in the first column is, seen from the subsystem in the second (screen coordinates: '
                    '"above" = smaller y). Strength hard/strong: follow it; medium: follow it unless it costs crossings; weak: a hint. '
                    'The build report names the violated hard, strong and medium ones (placement penalty).', '']
        out += ['## Subsystems and catalogue entities', ''] + self.matches()
        out += ['', '## Directions between subsystems', '', '| subsystem | is | strength | evidence |', '|---|---|---|---|']
        out += [d.line() for d in self.subsystem_directions] or ['| (none) | | | |']
        cycles = self.cycles()
        if cycles: out += ['', '## Cycles', ''] + cycles
        attach = self.attachments()
        if attach: out += ['', '## Attachments (where a subsystem joins another: merge candidates)', ''] + attach
        fb = self.fallback()
        if fb: out += ['', '## Subsystems without a direction', ''] + fb
        if self.metabolite_directions:
            out += ['', '## Directions between metabolites (textbook conventions)', '', '| metabolite | is | strength | note |', '|---|---|---|---|']
            out += [d.line() for d in self.metabolite_directions]
        s = '\n'.join(out) + '\n'
        if max_chars and len(s) > max_chars:
            cut = s.rfind('\n', 0, max_chars - 40)
            s = s[:cut] + '\n(... shortened)\n'
        return s

    def to_dict(self):
        return dict(catalogue_version=self.cat.version, taxa=self.taxa,
                    subsystems={k: dict(title=self.titles[k], entities=self.entities[k], matched_by=self.how[k]) for k in self.keys},
                    subsystem_directions=[d.to_dict() for d in self.subsystem_directions],
                    metabolite_directions=[d.to_dict() for d in self.metabolite_directions],
                    cycles=self.cycles(), attachments=self.attachments(), without_direction=self.fallback())


# ---------------------------------------------------------------------- subsystems of a map, in its three forms

def subsystems_of_analysis(analysis, key=None):
    """From `inspect_map.analyse(...)`: the title labels (or "Module n") -> their primary metabolites. `key(title)` renames."""
    out = {}
    for c in analysis['components']:
        title = c['suggested_subsystem'] or f"Module {c['index'] + 1}"
        k = key(title) if key else title
        e = out.setdefault(k, dict(title=title, metabolites=[]))
        e['metabolites'] = sorted(set(e['metabolites']) | set(c.get('primary_metabolites', [])))
    return out

def subsystems_of_layout(L):
    """From a Layout after a spec was applied: the zones (their map title label, else their title) -> primary metabolites."""
    mets = defaultdict(set)
    for rid, z in L.rxn_zone.items():
        for _, n in L.info[rid]['subs'] + L.info[rid]['prods']:
            if L.N[n]['node_is_primary']: mets[z].add(L.N[n]['bigg_id'])
    out = {}
    for z, st in L.zone_style.items():
        if z not in mets: continue
        out[z] = dict(title=st.get('label') or st.get('title') or z, metabolites=sorted(mets[z]))
    return out

def subsystems_of_spec(spec, L):
    """From a spec's `zones`/`plan` (before any layout exists): zone name -> primary metabolites of its planned reactions."""
    labels = {z['name']: z.get('label') or z.get('title') or z['name'] for z in spec.get('zones', [])}
    out = {}
    for z, e in spec.get('plan', {}).items():
        mets = set()
        for rid in e.get('reactions', []):
            if rid not in L.info: continue
            for _, n in L.info[rid]['subs'] + L.info[rid]['prods']:
                if L.N[n]['node_is_primary']: mets.add(L.N[n]['bigg_id'])
        out[z] = dict(title=labels.get(z, z), metabolites=sorted(mets))
    return out


# ---------------------------------------------------------------------- geometry and check

def layout_geometry(L):
    """(zone -> centre, BiGG id -> positions of its primary nodes) of a finished Layout (after L.finish())."""
    pts = defaultdict(list); pos = defaultdict(list)
    for rid, z in L.rxn_zone.items():
        I = L.info[rid]
        for n in [I['mid']] + [n for _, n in I['subs'] + I['prods'] if L.N[n]['node_is_primary']]:
            if n in L.N and 'x' in L.N[n]: pts[z].append((L.N[n]['x'], L.N[n]['y']))
    for n, v in L.N.items():
        if v['node_type'] == 'metabolite' and v['node_is_primary'] and 'x' in v: pos[v['bigg_id']].append((v['x'], v['y']))
    centres = {z: (sum(p[0] for p in v) / len(v), sum(p[1] for p in v) / len(v)) for z, v in pts.items() if v}
    return centres, dict(pos)

def spec_geometry(spec, L=None):
    """(zone -> centre of its blocks, BiGG id -> node positions) straight from a spec (origin + nodes + rings), without a build."""
    pts = defaultdict(list); pos = defaultdict(list)
    for b in spec.get('layouts', []):
        z = b.get('zone') or b.get('block'); ox, oy = b.get('origin', [0, 0])
        for n, (x, y) in b.get('nodes', {}).items():
            p = (ox + x, oy + y); pts[z].append(p)
            if L is not None and n in L.N: pos[L.N[n]['bigg_id']].append(p)
        for ring in b.get('rings', []):
            cx, cy = ring.get('center', (0, 0)); r = ring.get('radius', 0)
            pts[z].append((ox + cx, oy + cy))
            for n, ang in ring.get('nodes', {}).items():
                p = (ox + cx + r * math.cos(math.radians(ang)), oy + cy + r * math.sin(math.radians(ang)))
                if L is not None and n in L.N: pos[L.N[n]['bigg_id']].append(p)
    centres = {z: (sum(p[0] for p in v) / len(v), sum(p[1] for p in v) / len(v)) for z, v in pts.items() if v}
    return centres, dict(pos)

def violated(d, a, b, tolerance=50.0):
    """Whether point b is NOT in direction d.direction of point a (screen coordinates, y down)."""
    dx, north = b[0] - a[0], a[1] - b[1]
    parts = d.direction.split('_')
    return (('above' in parts and north <= tolerance) or ('below' in parts and north >= -tolerance)
            or ('right' in parts and dx <= tolerance) or ('left' in parts and dx >= -tolerance))

def check(arr, centres, positions):
    """The violated hard, strong and medium directions of a layout, and their penalty.
    Returns dict(violations=[str], penalty=int, checked=int, directions=int)."""
    problems, penalty, checked = [], 0, 0
    for d in arr.subsystem_directions:
        if d.strength in PENALTY and d.anchor in centres and d.target in centres:
            checked += 1
            if violated(d, centres[d.anchor], centres[d.target]):
                problems.append(f'{d.target} should be {d.where()} ({d.strength}; {d.evidence})'); penalty += PENALTY[d.strength]
    for d in arr.metabolite_directions:
        a = positions.get(d.anchor, []); b = positions.get(d.target, [])
        if a and b:
            checked += 1
            if all(violated(d, p, q) for p in a for q in b):
                problems.append(f'{d.target} should be {d.where()} ({d.strength}; textbook convention)'); penalty += PENALTY[d.strength]
    return dict(violations=problems, penalty=penalty, checked=checked,
                directions=len([d for d in arr.subsystem_directions if d.strength in PENALTY]) + len(arr.metabolite_directions))

def check_layout(L, taxa='all', catalogue=None):
    """Placement check of a finished Layout, or None when there is no catalogue or no subsystem matches it."""
    cat = catalogue or Catalogue.load()
    if cat is None: return None
    subs = subsystems_of_layout(L)
    if not subs: return None
    arr = Arrangement(subs, cat, taxa)
    if not any(arr.entities.values()): return None
    out = check(arr, *layout_geometry(L))
    out['subsystems'] = {k: arr.entities[k] for k in arr.keys}
    return out
