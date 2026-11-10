"""Declarative layout specification (pure JSON data - no code is ever executed).

A *spec* describes where things go: subsystems, node coordinates, per-reaction drawing options, hub merges, titles.
This module provides
  * SCHEMA                 JSON-Schema (draft 2020-12 subset) of a spec / spec fragment - single source of truth
  * validate_schema()      strict validator for that subset (no third-party packages)
  * parse_json_strict()    JSON parser that rejects duplicate keys, NaN/Infinity and oversized input
  * merge_fragment()       apply a partial spec ("patch") produced in a chat turn to the current spec
  * compile_spec()         semantic checks against the map + application onto an engine.Layout

Security model: the only inputs are JSON values that passed the schema; they are used as numbers, enum values and
ids that are looked up in the map. Nothing from a spec is ever evaluated, imported, or used as a file path.
"""
import copy, json, math, re
import numpy as np
from ._io import read_text, write_json

class SpecError(Exception):
    pass

LIMITS = dict(max_bytes=1_500_000, coord=20000, max_nodes=6000, max_reactions=3000, max_depth=14)

# ------------------------------------------------------------------------------------------------ schema
_TEXT = r'^[^\x00-\x1f\x7f<>]*$'                      # printable, no angle brackets
_NOTES = r'^[^\x00-\x08\x0b-\x1f\x7f<>]*$'            # like _TEXT but newlines/tabs allowed
_ZONE = r"^[A-Za-z0-9][A-Za-z0-9 .,()/+:'_\-]{0,79}$"
_REF = r'^(\d{1,9}|@[A-Za-z][A-Za-z0-9_]{0,30})$'     # node id of the map, or @alias created by "splits"
_BIGG = r'^[A-Za-z0-9_.\-()]{1,60}$'
_C = LIMITS['coord']
# what the patterns ask for, in the error messages (models write e.g. BiGG ids where node ids belong)
_FORMATS = {_REF: 'a node id of analysis.md (digits, e.g. "1234", not a BiGG id) or an @alias of "splits"',
            r'^\d{1,9}$': 'an id of analysis.md (digits, e.g. "1234", not a BiGG id)',
            _ZONE: 'a name of letters, digits, spaces and .,()/+:\'_-', _BIGG: 'a BiGG id',
            _TEXT: 'one line without < and >', _NOTES: 'text without < and >',
            r'^#[0-9a-fA-F]{6}$': 'a colour #rrggbb with 6 hex digits, e.g. "#eef3fb"'}

def _bad_format(path, what, value, pattern):
    hint = _FORMATS.get(pattern)
    return f'{path}: {what}{value[:40]!r} has a forbidden format' + (f' (expected {hint})' if hint else '')

def _num(lo, hi): return {'type': 'number', 'minimum': lo, 'maximum': hi}
_POINT = {'type': 'array', 'minItems': 2, 'maxItems': 2, 'items': _num(-_C, _C)}
_VEC = {'type': 'array', 'minItems': 2, 'maxItems': 2, 'items': _num(-5000, 5000)}
_REFS = lambda n: {'type': 'array', 'minItems': 1, 'maxItems': n, 'items': {'$ref': '#/$defs/nodeRef'}}

SCHEMA = {
    '$schema': 'https://json-schema.org/draft/2020-12/schema',
    'title': 'escher_handlayout layout spec v1 (also valid as a partial "fragment" in chat turns)',
    'type': 'object', 'additionalProperties': False,
    'properties': {
        'spec_version': {'const': 1},
        'notes': {'type': 'string', 'maxLength': 4000, 'pattern': _NOTES},
        'defaults': {'type': 'object', 'additionalProperties': False,
                     'properties': {'bend': _num(0, 400)}},
        'zones': {'type': 'array', 'maxItems': 40, 'items': {'$ref': '#/$defs/zone'}},
        'plan': {'type': 'object', 'maxProperties': 40, 'propertyNames': {'pattern': _ZONE},
                 'additionalProperties': {'$ref': '#/$defs/planEntry'}},
        'splits': {'type': 'array', 'maxItems': 200, 'items': {'$ref': '#/$defs/split'}},
        'merges': {'type': 'array', 'maxItems': 200, 'items': {'$ref': '#/$defs/merge'}},
        'layouts': {'type': 'array', 'maxItems': 200, 'items': {'$ref': '#/$defs/block'}},
        'titles': {'type': 'object', 'maxProperties': 40, 'propertyNames': {'pattern': _ZONE},
                   'additionalProperties': _POINT},
        'small_titles': {'type': 'array', 'maxItems': 100, 'items': {'$ref': '#/$defs/smallTitle'}},
    },
    '$defs': {
        'nodeRef': {'type': 'string', 'pattern': _REF},
        'zone': {'type': 'object', 'additionalProperties': False, 'required': ['name'],
                 'properties': {'name': {'type': 'string', 'pattern': _ZONE},
                                'color': {'type': 'string', 'pattern': r'^#[0-9a-fA-F]{6}$'},
                                'title': {'type': 'string', 'maxLength': 100, 'pattern': _TEXT},
                                'label': {'type': ['string', 'null'], 'maxLength': 100, 'pattern': _TEXT},
                                'bend': _num(0, 400)}},
        'planEntry': {'type': 'object', 'additionalProperties': False,
                      'properties': {'origin': _POINT, 'size': {'type': 'array', 'minItems': 2, 'maxItems': 2, 'items': _num(100, 40000)},
                                     'reactions': {'type': 'array', 'maxItems': 3000, 'items': {'type': 'string', 'pattern': r'^\d{1,9}$'}},
                                     'description': {'type': 'string', 'maxLength': 300, 'pattern': _TEXT}}},
        'split': {'type': 'object', 'additionalProperties': False, 'required': ['alias', 'node', 'reaction'],
                  'properties': {'alias': {'type': 'string', 'pattern': r'^@[A-Za-z][A-Za-z0-9_]{0,30}$'},
                                 'node': {'type': 'string', 'pattern': r'^\d{1,9}$'},
                                 'reaction': {'type': 'string', 'pattern': r'^\d{1,9}$'}}},
        'merge': {'type': 'object', 'additionalProperties': False, 'required': ['keep', 'dups'],
                  'properties': {'keep': {'$ref': '#/$defs/nodeRef'}, 'dups': _REFS(30)}},
        'smallTitle': {'type': 'object', 'additionalProperties': False, 'required': ['reaction', 'zone', 'text'],
                       'properties': {'reaction': {'type': 'string', 'pattern': _BIGG}, 'zone': {'type': 'string', 'pattern': _ZONE},
                                      'text': {'type': 'string', 'maxLength': 100, 'pattern': _TEXT},
                                      'reuse_label': {'type': ['string', 'null'], 'maxLength': 100, 'pattern': _TEXT}}},
        'pointExpr': {'oneOf': [
            _POINT,
            {'type': 'object', 'additionalProperties': False, 'required': ['node'], 'properties': {'node': {'$ref': '#/$defs/nodeRef'}}},
            {'type': 'object', 'additionalProperties': False, 'required': ['mean'], 'properties': {'mean': _REFS(8)}},
            {'type': 'object', 'additionalProperties': False, 'required': ['lerp'],
             'properties': {'lerp': {'type': 'array', 'minItems': 3, 'maxItems': 3,
                                     'prefixItems': [{'$ref': '#/$defs/pointExpr'}, {'$ref': '#/$defs/pointExpr'}, _num(-2, 3)]}}}]},
        'vecExpr': {'oneOf': [
            _VEC,
            {'type': 'object', 'additionalProperties': False, 'required': ['perp', 'offset'],
             'properties': {'perp': {'type': 'array', 'minItems': 2, 'maxItems': 2, 'items': {'$ref': '#/$defs/nodeRef'}},
                            'offset': _num(-5000, 5000)}}]},
        'axisExpr': {'oneOf': [
            _VEC,
            {'type': 'object', 'additionalProperties': False, 'required': ['from', 'to'],
             'properties': {'from': {'$ref': '#/$defs/pointExpr'}, 'to': {'$ref': '#/$defs/pointExpr'}}}]},
        'arc': {'oneOf': [
            {'type': 'object', 'additionalProperties': False, 'required': ['ring', 'offset'],
             'properties': {'ring': {'type': 'string', 'pattern': _ZONE}, 'offset': _num(-2000, 2000)}},
            {'type': 'object', 'additionalProperties': False, 'required': ['center', 'offset'],
             'properties': {'center': _POINT, 'offset': _num(-2000, 2000)}}]},
        'ring': {'type': 'object', 'additionalProperties': False, 'required': ['name', 'center', 'radius', 'nodes'],
                 'properties': {'name': {'type': 'string', 'pattern': _ZONE}, 'center': _POINT, 'radius': _num(50, 8000),
                                'nodes': {'type': 'object', 'minProperties': 2, 'maxProperties': 60, 'propertyNames': {'pattern': _REF},
                                          'additionalProperties': _num(-720, 720)}}},
        'reaction': {'type': 'object', 'additionalProperties': False, 'required': ['id'],
                     'properties': {
                         'id': {'type': 'string', 'pattern': r'^\d{1,9}$'},
                         'remove': {'type': 'boolean'},
                         'zone': {'type': 'string', 'pattern': _ZONE},
                         'side': {'enum': ['L', 'R', 'U', 'D', 'out', 'in']},
                         'dist': _num(20, 800), 'secdist': _num(20, 1200), 'gap': _num(5, 300), 'step': _num(20, 400),
                         'bend': _num(0, 400), 'bow': _num(-0.5, 0.5), 'curve': {'type': 'boolean'},
                         'shift': {'$ref': '#/$defs/vecExpr'}, 'axis': {'$ref': '#/$defs/axisExpr'}, 'mid': {'$ref': '#/$defs/pointExpr'},
                         'center': _POINT, 'arc': {'$ref': '#/$defs/arc'},
                         'arc_nodes': {'type': 'array', 'minItems': 2, 'maxItems': 2, 'items': {'$ref': '#/$defs/nodeRef'}}}},
        'block': {'type': 'object', 'additionalProperties': False, 'required': ['block'],
                  'properties': {'block': {'type': 'string', 'pattern': _ZONE},
                                 'zone': {'type': 'string', 'pattern': _ZONE},
                                 'origin': _POINT,
                                 'nodes': {'type': 'object', 'maxProperties': 2000, 'propertyNames': {'pattern': _REF}, 'additionalProperties': _POINT},
                                 'remove_nodes': {'type': 'array', 'maxItems': 2000, 'items': {'$ref': '#/$defs/nodeRef'}},
                                 'rings': {'type': 'array', 'maxItems': 20, 'items': {'$ref': '#/$defs/ring'}},
                                 'reactions': {'type': 'array', 'maxItems': 3000, 'items': {'$ref': '#/$defs/reaction'}}}},
    }}


# ------------------------------------------------------------------------------------------------ validator
def _is(v, t):
    if t == 'object': return isinstance(v, dict)
    if t == 'array': return isinstance(v, list)
    if t == 'string': return isinstance(v, str)
    if t == 'boolean': return isinstance(v, bool)
    if t == 'number': return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
    if t == 'integer': return isinstance(v, int) and not isinstance(v, bool)
    if t == 'null': return v is None
    return False

def _resolve(ref, root):
    node = root
    for part in ref.lstrip('#/').split('/'): node = node[part]
    return node

def validate_schema(inst, schema=SCHEMA, root=SCHEMA, path='$', errs=None, depth=0):
    """Validate `inst` against the supported JSON-Schema subset. Returns the list of error strings."""
    errs = [] if errs is None else errs
    if depth > LIMITS['max_depth'] + 6: errs.append(f'{path}: nesting too deep'); return errs
    if '$ref' in schema: return validate_schema(inst, _resolve(schema['$ref'], root), root, path, errs, depth + 1)
    if 'oneOf' in schema:
        ok = []
        for i, s in enumerate(schema['oneOf']):
            e = validate_schema(inst, s, root, path, [], depth + 1)
            if not e: ok.append(i)
        if len(ok) != 1: errs.append(f'{path}: does not match exactly one allowed form ({len(ok)} matched)')
        return errs
    if 'const' in schema and inst != schema['const']: errs.append(f"{path}: must be {schema['const']!r}"); return errs
    if 'enum' in schema and inst not in schema['enum']: errs.append(f"{path}: {inst!r} not in {schema['enum']}"); return errs
    t = schema.get('type')
    if t is not None:
        ts = t if isinstance(t, list) else [t]
        if not any(_is(inst, x) for x in ts): errs.append(f'{path}: expected {"/".join(ts)}, got {type(inst).__name__}'); return errs
    if isinstance(inst, bool): return errs
    if isinstance(inst, (int, float)):
        if 'minimum' in schema and inst < schema['minimum']: errs.append(f"{path}: {inst} < minimum {schema['minimum']}")
        if 'maximum' in schema and inst > schema['maximum']: errs.append(f"{path}: {inst} > maximum {schema['maximum']}")
    elif isinstance(inst, str):
        if 'maxLength' in schema and len(inst) > schema['maxLength']: errs.append(f"{path}: longer than {schema['maxLength']} characters")
        if 'pattern' in schema and not re.search(schema['pattern'], inst): errs.append(_bad_format(path, '', inst, schema['pattern']))
    elif isinstance(inst, list):
        if 'minItems' in schema and len(inst) < schema['minItems']: errs.append(f"{path}: fewer than {schema['minItems']} items")
        if 'maxItems' in schema and len(inst) > schema['maxItems']: errs.append(f"{path}: more than {schema['maxItems']} items")
        for i, x in enumerate(inst):
            if 'prefixItems' in schema and i < len(schema['prefixItems']):
                validate_schema(x, schema['prefixItems'][i], root, f'{path}[{i}]', errs, depth + 1)
            elif 'items' in schema:
                validate_schema(x, schema['items'], root, f'{path}[{i}]', errs, depth + 1)
    elif isinstance(inst, dict):
        props = schema.get('properties', {})
        for k in schema.get('required', []):
            if k not in inst: errs.append(f'{path}: missing required "{k}"')
        if 'minProperties' in schema and len(inst) < schema['minProperties']: errs.append(f"{path}: fewer than {schema['minProperties']} entries")
        if 'maxProperties' in schema and len(inst) > schema['maxProperties']: errs.append(f"{path}: more than {schema['maxProperties']} entries")
        for k, v in inst.items():
            if 'propertyNames' in schema and not re.search(schema['propertyNames']['pattern'], k):
                errs.append(_bad_format(path, 'key ', k, schema['propertyNames']['pattern'])); continue
            if k in props: validate_schema(v, props[k], root, f'{path}.{k}', errs, depth + 1)
            elif schema.get('additionalProperties') is False: errs.append(f'{path}: unknown key "{k[:40]}"')
            elif isinstance(schema.get('additionalProperties'), dict):
                validate_schema(v, schema['additionalProperties'], root, f'{path}.{k}', errs, depth + 1)
    return errs


def parse_json_strict(text, max_bytes=None):
    """json.loads that rejects duplicate keys, NaN/Infinity and oversized input. Raises SpecError."""
    max_bytes = max_bytes or LIMITS['max_bytes']
    if len(text.encode('utf-8', 'ignore')) > max_bytes: raise SpecError(f'input larger than {max_bytes} bytes')
    def pairs(items):
        d = {}
        for k, v in items:
            if k in d: raise SpecError(f'duplicate key {k[:40]!r}')
            d[k] = v
        return d
    def bad(c): raise SpecError(f'non-finite number {c} is not allowed')
    try:
        return json.loads(text, object_pairs_hook=pairs, parse_constant=bad)
    except SpecError:
        raise
    except (ValueError, RecursionError) as e:
        raise SpecError(f'not valid JSON: {e}')


def load_spec_file(path):
    text = read_text(path)
    spec = parse_json_strict(text)
    errs = validate_schema(spec)
    if errs: raise SpecError('spec violates the schema:\n  ' + '\n  '.join(errs[:60]))
    return spec


def write_schema_file(path):
    write_json(path, SCHEMA)


# ------------------------------------------------------------------------------------------------ patching
def merge_fragment(spec, frag):
    """Return a NEW spec = `spec` updated by the partial spec `frag` (both must already satisfy the schema).

    zones / splits / merges / small_titles / layouts / rings / reactions are merged by key (name, alias, keep, reaction,
    block, ring name, reaction id); dictionaries (plan, titles, nodes, defaults) are updated key by key;
    {"id": "..", "remove": true} deletes a reaction entry, "remove_nodes" deletes node coordinates."""
    out = copy.deepcopy(spec)
    out.setdefault('spec_version', 1)
    for k in ('notes',):
        if k in frag: out[k] = frag[k]
    if 'defaults' in frag: out.setdefault('defaults', {}).update(frag['defaults'])
    def by_key(section, key):
        if section not in frag: return
        cur = out.setdefault(section, []); idx = {e[key]: i for i, e in enumerate(cur)}
        for e in frag[section]:
            if e[key] in idx: cur[idx[e[key]]] = copy.deepcopy(e)
            else: idx[e[key]] = len(cur); cur.append(copy.deepcopy(e))
    by_key('zones', 'name'); by_key('splits', 'alias'); by_key('small_titles', 'reaction')
    if 'merges' in frag:
        cur = out.setdefault('merges', []); idx = {e['keep']: i for i, e in enumerate(cur)}
        for e in frag['merges']:
            if e['keep'] in idx: cur[idx[e['keep']]] = copy.deepcopy(e)
            else: idx[e['keep']] = len(cur); cur.append(copy.deepcopy(e))
    for sect in ('plan', 'titles'):
        if sect in frag:
            cur = out.setdefault(sect, {})
            for k, v in frag[sect].items():
                if isinstance(v, dict) and isinstance(cur.get(k), dict): cur[k].update(copy.deepcopy(v))
                else: cur[k] = copy.deepcopy(v)
    if 'layouts' in frag:
        blocks = out.setdefault('layouts', []); idx = {b['block']: i for i, b in enumerate(blocks)}
        for fb in frag['layouts']:
            if fb['block'] not in idx:
                nb = {k: copy.deepcopy(v) for k, v in fb.items() if k not in ('remove_nodes',)}
                nb['reactions'] = [r for r in nb.get('reactions', []) if not r.get('remove')]
                idx[fb['block']] = len(blocks); blocks.append(nb); continue
            b = blocks[idx[fb['block']]]
            for k in ('zone', 'origin'):
                if k in fb: b[k] = copy.deepcopy(fb[k])
            nodes = b.setdefault('nodes', {})
            for k, v in fb.get('nodes', {}).items(): nodes[k] = copy.deepcopy(v)
            for k in fb.get('remove_nodes', []): nodes.pop(k, None)
            if 'rings' in fb:
                rings = b.setdefault('rings', []); ri = {r['name']: i for i, r in enumerate(rings)}
                for r in fb['rings']:
                    if r['name'] in ri: rings[ri[r['name']]] = copy.deepcopy(r)
                    else: ri[r['name']] = len(rings); rings.append(copy.deepcopy(r))
            if 'reactions' in fb:
                rx = b.setdefault('reactions', []); xi = {r['id']: i for i, r in enumerate(rx)}
                for r in fb['reactions']:
                    if r.get('remove'):
                        if r['id'] in xi: rx[xi[r['id']]] = None
                        continue
                    if r['id'] in xi and rx[xi[r['id']]] is not None: rx[xi[r['id']]] = copy.deepcopy(r)
                    else: xi[r['id']] = len(rx); rx.append(copy.deepcopy(r))
                b['reactions'] = [r for r in rx if r is not None]
    return out


# ------------------------------------------------------------------------------------------------ compilation
def compile_spec(spec, L, partial=False):
    """Check `spec` against the map held by Layout `L` and apply it to `L`.

    Returns (problems, warnings). With partial=False a non-empty problem list raises SpecError (the layout is then
    incomplete); with partial=True (used between chat turns) completeness is not required and unresolved
    expressions are skipped silently. Always call on a FRESH Layout."""
    P, W = [], []
    N, R = L.N, L.R
    C = LIMITS['coord']
    nprim = lambda r: [n for _, n in L.info[r]['subs'] + L.info[r]['prods'] if N[n]['node_is_primary']]

    # ---- zones
    zones = {}
    for z in spec.get('zones', []):
        if z['name'] in zones: P.append(f"zone '{z['name']}' declared twice")
        zones[z['name']] = z
    base_bend = spec.get('defaults', {}).get('bend', 45)
    for z in zones.values():
        L.zone(z['name'], z.get('color'), z.get('title') or z['name'], z.get('label'))

    # ---- splits then merges (ids of the map; aliases defined by splits)
    alias, merged_into = {}, {}
    def missing(i, what):
        hint = (f" ({i} is the id of reaction {R[i]['bigg_id']}: give the node ids of its metabolites, "
                f"name[id] in the reactions of analysis.md)") if i in R else ''
        P.append(f'{what}: node {i} does not exist{hint}')
    def prim_node(i, what):
        if i not in N: missing(i, what); return False
        if N[i]['node_type'] != 'metabolite' or not N[i]['node_is_primary']: P.append(f'{what}: node {i} is not a primary metabolite node'); return False
        return True
    def ref(r, what):
        i = alias.get(r) if r.startswith('@') else r
        if i is None: P.append(f'{what}: unknown alias {r}'); return None
        if i in merged_into: P.append(f'{what}: node {i} was merged into {merged_into[i]}; use that id'); return None
        if i not in N: missing(i, what); return None
        return i
    for s in spec.get('splits', []):
        if s['alias'] in alias: P.append(f"split alias {s['alias']} defined twice"); continue
        if not prim_node(s['node'], f"split {s['alias']}"): continue
        if s['reaction'] not in R: P.append(f"split {s['alias']}: reaction {s['reaction']} does not exist"); continue
        if s['node'] not in {n for _, n in L.info[s['reaction']]['subs'] + L.info[s['reaction']]['prods']}:
            P.append(f"split {s['alias']}: reaction {s['reaction']} does not use node {s['node']}"); continue
        alias[s['alias']] = L.split(s['node'], s['reaction'])
    seen_merge = set()
    for m in spec.get('merges', []):
        keep = ref(m['keep'], 'merge keep'); dups = [ref(d, f"merge into {m['keep']}") for d in m['dups']]
        if keep is None or any(d is None for d in dups): continue
        ids = [keep] + dups
        if len(set(ids)) != len(ids) or seen_merge & set(ids): P.append(f"merge {m['keep']}: ids repeated or already used in another merge"); continue
        if not all(prim_node(i, f"merge {m['keep']}") for i in ids): continue
        names = {N[i]['bigg_id'] for i in ids}
        if len(names) != 1: P.append(f"merge {m['keep']}: different metabolites cannot be merged: {sorted(names)}"); continue
        L.merge(keep, *dups); seen_merge |= set(ids)
        for d in dups: merged_into[d] = keep

    # ---- blocks: node coordinates and rings
    placed = {}          # node id -> block name
    ring_center = {}     # ring name -> global centre
    block_nodes, block_zone, block_origin = {}, {}, {}
    names = set()
    for b in spec.get('layouts', []):
        bn = b['block']
        if bn in names: P.append(f"block '{bn}' appears twice"); continue
        names.add(bn)
        o = b.get('origin', [0, 0]); zone = b.get('zone')
        block_origin[bn] = o; block_zone[bn] = zone
        if zone is None: P.append(f"block '{bn}': missing zone")
        elif zone not in zones: P.append(f"block '{bn}': zone '{zone}' is not declared in zones")
        block_nodes[bn] = []
        def put(i, x, y, what):
            if not (math.isfinite(x) and math.isfinite(y) and abs(x) <= C and abs(y) <= C): P.append(f'{what}: coordinate out of range'); return
            if i in placed: P.append(f"{what}: node {i} is already placed in block '{placed[i]}'"); return
            placed[i] = bn; L.pos[i] = (float(x), float(y)); block_nodes[bn].append(i)
        for k, v in b.get('nodes', {}).items():
            i = ref(k, f"block '{bn}' nodes")
            if i is not None and prim_node(i, f"block '{bn}' nodes"): put(i, o[0] + v[0], o[1] + v[1], f"block '{bn}' node {k}")
        for rg in b.get('rings', []):
            if rg['name'] in ring_center: P.append(f"ring '{rg['name']}' defined twice"); continue
            c = np.array((o[0] + rg['center'][0], o[1] + rg['center'][1])); ring_center[rg['name']] = c
            for k, ang in rg['nodes'].items():
                i = ref(k, f"ring '{rg['name']}'")
                if i is not None and prim_node(i, f"ring '{rg['name']}'"):
                    put(i, c[0] + rg['radius'] * np.cos(np.radians(ang)), c[1] + rg['radius'] * np.sin(np.radians(ang)), f"ring '{rg['name']}' node {k}")
    if len(placed) > LIMITS['max_nodes']: P.append('too many placed nodes')

    # ---- expression evaluation
    class Unresolved(Exception): pass
    def pos_of(r, what):
        i = ref(r, what)
        if i is None: raise Unresolved()
        if i not in L.pos:
            if not partial: P.append(f'{what}: node {r} has no coordinates yet')
            raise Unresolved()
        return np.array(L.pos[i], float)
    def ev_point(e, o, what, depth=0):
        if depth > 4: P.append(f'{what}: expression nested too deeply'); raise Unresolved()
        if isinstance(e, list): return np.array((o[0] + e[0], o[1] + e[1]), float)
        if 'node' in e: return pos_of(e['node'], what)
        if 'mean' in e: return np.mean([pos_of(r, what) for r in e['mean']], 0)
        a, bb, t = e['lerp']
        p1, p2 = ev_point(a, o, what, depth + 1), ev_point(bb, o, what, depth + 1)
        return p1 + t * (p2 - p1)
    def ev_vec(e, what):
        if isinstance(e, list): return tuple(e)
        pa, pb = pos_of(e['perp'][0], what), pos_of(e['perp'][1], what)
        if np.linalg.norm(pb - pa) < 1e-9: P.append(f'{what}: perp reference nodes coincide'); raise Unresolved()
        u = (pb - pa) / np.linalg.norm(pb - pa); n = np.array([-u[1], u[0]])
        return tuple(n * e['offset'])
    def ev_axis(e, o, what):
        if isinstance(e, list): v = np.array(e, float)
        else: v = ev_point(e['to'], o, what) - ev_point(e['from'], o, what)
        if np.linalg.norm(v) < 1e-9: P.append(f'{what}: axis has zero length'); raise Unresolved()
        return tuple(v)

    # ---- reactions
    done = {}
    for b in spec.get('layouts', []):
        bn = b['block']; o = b.get('origin', [0, 0])
        for r in b.get('reactions', []):
            rid = r['id']; what = f"reaction {rid} (block '{bn}')"
            if rid not in R: P.append(f'{what}: reaction does not exist'); continue
            if rid in done: P.append(f"{what}: already placed in block '{done[rid]}'"); continue
            zone = r.get('zone', b.get('zone'))
            if zone not in zones: P.append(f"{what}: zone {zone!r} not declared"); continue
            missing = [n for n in nprim(rid) if n not in L.pos]
            if missing:
                if not partial: P.append(f"{what} ({R[rid]['bigg_id']}): primary nodes without coordinates: {', '.join(sorted(set(N[n]['bigg_id'] + '#' + n for n in missing)))}")
                continue
            kw = {}
            try:
                for k in ('side', 'dist', 'secdist', 'gap', 'step', 'bow', 'curve'):
                    if k in r: kw[k] = r[k]
                kw['bend'] = r.get('bend', zones[zone].get('bend', base_bend))
                if 'shift' in r: kw['shift'] = ev_vec(r['shift'], what)
                if 'axis' in r: kw['axis'] = ev_axis(r['axis'], o, what)
                if 'mid' in r: kw['mid'] = tuple(ev_point(r['mid'], o, what))
                if 'center' in r: kw['center'] = (o[0] + r['center'][0], o[1] + r['center'][1])
                if 'arc' in r:
                    a = r['arc']
                    if 'ring' in a:
                        if a['ring'] not in ring_center: P.append(f"{what}: unknown ring '{a['ring']}'"); continue
                        cen = tuple(ring_center[a['ring']])
                    else: cen = (o[0] + a['center'][0], o[1] + a['center'][1])
                    kw['arc'] = (cen, a['offset'])
                    if 'arc_nodes' in r:
                        an = [ref(x, what) for x in r['arc_nodes']]
                        if None in an: continue
                        if any(x not in nprim(rid) for x in an): P.append(f'{what}: arc_nodes must be primary nodes of this reaction'); continue
                        kw['arc_nodes'] = tuple(an)
                elif 'arc_nodes' in r: P.append(f'{what}: arc_nodes requires arc'); continue
                if kw.get('side') in ('out', 'in') and 'arc' not in kw and 'center' not in kw: P.append(f"{what}: side '{kw['side']}' needs arc or center"); continue
                subP = [n for _, n in L.info[rid]['subs'] if N[n]['node_is_primary']]; prodP = [n for _, n in L.info[rid]['prods'] if N[n]['node_is_primary']]
                if 'arc' not in kw and (not subP or not prodP) and not ('axis' in kw and 'mid' in kw):
                    P.append(f"{what} ({R[rid]['bigg_id']}): one side has no primary node, so both axis and mid are required"); continue
                if 'arc' in kw and 'arc_nodes' not in kw and (not subP or not prodP): P.append(f'{what}: arc needs arc_nodes here'); continue
            except Unresolved:
                continue
            try:
                L.reaction(rid, zone, **kw)
            except Exception as e:                         # engine refused (e.g. topology) - report instead of crashing
                P.append(f'{what}: {e}'); continue
            done[rid] = bn
    # ---- titles / small titles
    for z, xy in spec.get('titles', {}).items():
        if z not in zones: P.append(f"title for undeclared zone '{z}'")
        else: L.titles[z] = (float(xy[0]), float(xy[1]))
    for s in spec.get('small_titles', []):
        if s['zone'] not in zones: P.append(f"small_title for undeclared zone '{s['zone']}'")
        else: L.small_titles[s['reaction']] = (s['zone'], s['text'], s.get('reuse_label'))
    # ---- completeness (full builds only)
    if not partial:
        un = [f"{rid} ({R[rid]['bigg_id']})" for rid in R if rid not in done]
        if un: P.append(f'{len(un)} reactions are not placed: ' + ', '.join(un[:80]) + (' ...' if len(un) > 80 else ''))
        for z in zones:
            if z not in L.titles and any(L.rxn_zone.get(r) == z for r in done): W.append(f"zone '{z}' has no title position (titles.{z})")
    L.spec_alias, L.spec_placed, L.spec_done = alias, placed, done      # for the chat harness (read-only info)
    # ---- soft checks against the plan
    plan = spec.get('plan', {})
    for bn, ids in block_nodes.items():
        z = block_zone.get(bn); reg = plan.get(z) if z else None
        if reg and 'origin' in reg and 'size' in reg and ids:
            xs = [L.pos[i][0] for i in ids]; ys = [L.pos[i][1] for i in ids]
            x0, y0 = reg['origin']; w, h = reg['size']
            if min(xs) < x0 - 300 or max(xs) > x0 + w + 300 or min(ys) < y0 - 300 or max(ys) > y0 + h + 300:
                W.append(f"block '{bn}': nodes lie outside the planned region of zone '{z}'")
    unused = [i for i in placed if not any(i in {n for _, n in L.info[r]['subs'] + L.info[r]['prods']} for r in R)]
    if unused: W.append(f'{len(unused)} placed nodes are not used by any reaction')
    if len(R) > LIMITS['max_reactions']: P.append('map has too many reactions for this tool')
    if P and not partial: raise SpecError('spec problems:\n  ' + '\n  '.join(P[:100]))
    return P, W
