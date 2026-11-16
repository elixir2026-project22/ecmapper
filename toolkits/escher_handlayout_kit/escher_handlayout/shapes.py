"""Coordinates by arithmetic: helpers for a script that writes a layout spec (PROMPT_AGENT.md).

A coding agent decides the design (which pathway is a column, which a ring, what hangs where) and lets these
functions turn it into exact numbers. They are plain geometry: nothing here looks at the map and decides a position.

    import sys; sys.path.insert(0, KIT_DIR)
    from escher_handlayout import shapes as S
    d = S.Draft('draft.escher.json')                       # ids by name
    gly = S.column(d.chain('HEX1', 'PGI', 'PFK', 'FBA'), step=400)
    spec = S.load('specs/draft.layout.json'); S.block(spec, 'Glycolysis')['nodes'].update(gly); S.save(spec, 'specs/draft.layout.json')

All coordinates are lists [x, y] (y points down), relative to the origin of the block they are put into.
"""
import math
from .engine import Layout
from .declspec import load_spec_file, format_spec, validate_schema
from ._io import write_text


# ----------------------------------------------------------------------------------------------- ids by name
class Draft:
    """The ids of a draft map by name. Reactions and metabolites are named by their BiGG ids; a name that several
    reactions have (the same reaction drawn in two subsystems) gives all of them, in the order of their ids."""
    def __init__(self, path):
        self.L = L = Layout(path); self.by_name = {}
        for rid in sorted(L.R, key=int): self.by_name.setdefault(L.R[rid]['bigg_id'], []).append(rid)

    def reactions(self, name):
        """All reaction ids with this BiGG id."""
        if name not in self.by_name: raise KeyError(f'no reaction {name!r} in the map')
        return list(self.by_name[name])

    def reaction(self, name, nth=0):
        """The id of the reaction `name` (the `nth` one if several have the name)."""
        return self.reactions(name)[nth]

    def primary(self, reaction, nth=0):
        """{BiGG id of the metabolite: node id} for the primary metabolites of a reaction (name or id)."""
        rid = reaction if reaction in self.L.R else self.reaction(reaction, nth); I = self.L.info[rid]
        return {self.L.N[n]['bigg_id']: n for _, n in I['subs'] + I['prods'] if self.L.N[n]['node_is_primary']}

    def node(self, reaction, metabolite, nth=0):
        """The node id of primary metabolite `metabolite` (BiGG id) in reaction `reaction`."""
        p = self.primary(reaction, nth)
        if metabolite not in p: raise KeyError(f'{reaction} has no primary metabolite {metabolite!r} (it has {sorted(p)})')
        return p[metabolite]

    def nodes_of(self, metabolite):
        """All primary node ids of a metabolite (BiGG id): the copies the draft draws, candidates for a merge."""
        return sorted((k for k, v in self.L.N.items() if v['node_type'] == 'metabolite' and v['node_is_primary'] and v['bigg_id'] == metabolite), key=int)

    def chain(self, *reactions):
        """The primary node ids along a chain of reactions (names or ids), each node once, in the order met: the
        substrate side of the first reaction first. For `column` / `row`."""
        out = []
        for r in reactions:
            rid = r if r in self.L.R else self.reaction(r); I = self.L.info[rid]
            prim = lambda xs: [n for _, n in xs if self.L.N[n]['node_is_primary']]
            s, p = prim(I['subs']), prim(I['prods'])
            if out and any(n in out for n in p) and not any(n in out for n in s): s, p = p, s      # drawn the other way round
            out += [n for n in s + p if n not in out]
        return out


# -------------------------------------------------------------------------------------------------- shapes
def column(ids, x=0, y=0, step=380):
    """The nodes below each other, from (x, y) downwards, `step` apart (a negative step: upwards)."""
    return {i: [x, y + k * step] for k, i in enumerate(ids)}


def row(ids, x=0, y=0, step=380):
    """The nodes beside each other, from (x, y) to the right, `step` apart (a negative step: to the left)."""
    return {i: [x + k * step, y] for k, i in enumerate(ids)}


def line(ids, start, end):
    """The nodes evenly on the straight line from `start` to `end` (both included)."""
    n = max(1, len(ids) - 1)
    return {i: [round(start[0] + (end[0] - start[0]) * k / n), round(start[1] + (end[1] - start[1]) * k / n)] for k, i in enumerate(ids)}


def ring_radius(n, spacing=420):
    """The radius of a ring whose `n` nodes are about `spacing` apart."""
    return round(n * spacing / (2 * math.pi))


def ring(name, ids, center=(0, 0), radius=None, start=-90, clockwise=True, spacing=420):
    """A `rings` entry of a block: the nodes evenly on a circle, the first at angle `start` (degrees; 0 = east, -90 =
    north, clockwise on the screen). Its reactions are drawn with `ring_reactions`."""
    radius = radius or ring_radius(len(ids), spacing); turn = 360.0 / len(ids) * (1 if clockwise else -1)
    return {'name': name, 'center': list(center), 'radius': radius, 'nodes': {i: round(start + k * turn, 2) for k, i in enumerate(ids)}}


def ring_reactions(name, reactions, side='out', dist=150):
    """Reaction entries for the reactions along ring `name`. `reactions`: ids, or (id, [node a, node b]) where the
    reaction touches other primary nodes besides its two ring nodes; an id listed again is a second lane further out."""
    out, seen = [], {}
    for r in reactions:
        rid, nodes = (r, None) if isinstance(r, str) else r
        k = seen.get(rid, 0); seen[rid] = k + 1
        e = {'id': rid, 'arc': {'ring': name, 'offset': 130 * k}, 'side': side, 'dist': dist}
        if nodes: e['arc_nodes'] = list(nodes)
        out.append(e)
    return out


def fan(ids, around, radius, start, end):
    """The nodes on an arc around the point `around`, from angle `start` to angle `end` (degrees, 0 = east, clockwise
    on the screen): the neighbours of a hub, each in its own direction."""
    n = max(1, len(ids) - 1)
    ang = [start + (end - start) * k / n for k in range(len(ids))] if len(ids) > 1 else [(start + end) / 2]
    return {i: [round(around[0] + radius * math.cos(math.radians(a))), round(around[1] + radius * math.sin(math.radians(a)))] for i, a in zip(ids, ang)}


def lanes(reactions, a, b, spacing=130, **keys):
    """Reaction entries for parallel reactions between the nodes `a` and `b`: lanes `spacing` apart, centred on the
    line a - b. `keys` (e.g. dist=90) go into every entry; cofactors hang on the outer side of each lane."""
    n = len(reactions); out = []
    for k, rid in enumerate(reactions):
        off = round((k - (n - 1) / 2) * spacing)
        e = dict({'id': rid}, **keys)
        if off: e['shift'] = {'perp': [a, b], 'offset': off}
        out.append(e)
    return out


def hub_gap(n_lanes, spacing=150, least=450):
    """How far apart two metabolites should be that `n_lanes` parallel reactions join (room for the lanes' bows)."""
    return max(least, round(n_lanes * spacing))


# ------------------------------------------------------------------------------------------- transformations
def move(nodes, dx=0, dy=0):
    return {i: [x + dx, y + dy] for i, (x, y) in nodes.items()}


def mirror(nodes, axis='x', about=0):
    """Mirror at the vertical line x = `about` (axis 'x': left and right swap) or the horizontal line y = `about` ('y')."""
    if axis not in ('x', 'y'): raise ValueError("axis must be 'x' or 'y'")
    return {i: ([2 * about - x, y] if axis == 'x' else [x, 2 * about - y]) for i, (x, y) in nodes.items()}


def rotate(nodes, degrees, about=(0, 0)):
    """Turn clockwise on the screen by `degrees` around the point `about`."""
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    return {i: [round(about[0] + c * (x - about[0]) - s * (y - about[1])), round(about[1] + s * (x - about[0]) + c * (y - about[1]))] for i, (x, y) in nodes.items()}


def scale(nodes, factor, about=(0, 0)):
    return {i: [round(about[0] + factor * (x - about[0])), round(about[1] + factor * (y - about[1]))] for i, (x, y) in nodes.items()}


def bbox(nodes):
    xs = [p[0] for p in nodes.values()]; ys = [p[1] for p in nodes.values()]
    return min(xs), min(ys), max(xs), max(ys)


# -------------------------------------------------------------------------------------------------- spec file
def load(path):
    return load_spec_file(path)


def save(spec, path):
    """Write the spec (checked against the schema first: ValueError lists what is wrong)."""
    errs = validate_schema(spec)
    if errs: raise ValueError('the spec violates the schema:\n  ' + '\n  '.join(errs[:40]))
    write_text(path, format_spec(spec))


def block(spec, name, zone=None, origin=None):
    """The block `name` of the spec; it is added if the spec has none (with `zone`, default the name, and `origin`)."""
    for b in spec.setdefault('layouts', []):
        if b['block'] == name:
            if origin is not None: b['origin'] = list(origin)
            return b
    b = {'block': name, 'zone': zone or name, 'origin': list(origin or (0, 0)), 'nodes': {}, 'reactions': []}
    spec['layouts'].append(b); return b


def set_reactions(b, entries):
    """Put reaction entries into block `b`: an entry replaces the one with the same id."""
    by_id = {r['id']: k for k, r in enumerate(b.setdefault('reactions', []))}
    for e in entries:
        if e['id'] in by_id: b['reactions'][by_id[e['id']]] = dict(e)
        else: by_id[e['id']] = len(b['reactions']); b['reactions'].append(dict(e))
