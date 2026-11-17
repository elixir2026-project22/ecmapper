"""An Escher map as the probes see it: reactions with their primary connection and side participants, metabolite
nodes, units (the subsystems of the draft, named by its title labels) and the settings of the rule set."""
import json
import math
import re
from collections import defaultdict

#: used where a checks file leaves a setting out (rules/*.checks.yaml, `settings`)
DEFAULT_SETTINGS = dict(
    compartment_suffix=r'_[a-z][a-z0-9]?$',
    currency=['atp', 'adp', 'amp', 'gtp', 'gdp', 'nad', 'nadh', 'nadp', 'nadph', 'fad', 'fadh2', 'coa', 'h2o', 'h', 'pi',
              'ppi', 'co2', 'o2', 'nh4', 'hco3', 'h2o2', 'h2'],
    carriers=['coa', 'nad', 'nadh', 'nadp', 'nadph', 'atp', 'adp', 'amp', 'gtp', 'gdp', 'fad', 'fadh2'],
    cofactor_pairs=[['atp', 'adp'], ['nad', 'nadh'], ['nadp', 'nadph'], ['fad', 'fadh2'], ['gtp', 'gdp']],
    exchange_prefixes=['EX_', 'DM_', 'SK_', 'sink_'],
    biomass_pattern='(?i)biomass',
    exit_label_prefix='→',
    primary_radius=20, secondary_radius=12,
)


def read_map(path):
    """(header, body) of an Escher map file."""
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    if not (isinstance(data, list) and len(data) == 2 and isinstance(data[1], dict) and 'reactions' in data[1]):
        raise ValueError(f'{path} is not an Escher map (a JSON list of a header and a body with reactions and nodes)')
    return data[0], data[1]


class Reaction:
    """A reaction of the map. `subs` / `prods`: the metabolite node ids on either side; `primary_subs` /
    `primary_prods`: those drawn as primary nodes; `marker[n]`: the marker node a metabolite node hangs on."""

    def __init__(self, rid, data):
        self.rid, self.data = rid, data
        self.bigg = data.get('bigg_id', rid)
        self.name = data.get('name') or ''
        self.reversible = bool(data.get('reversibility'))
        self.subs, self.prods, self.marker, self.segments = [], [], {}, []
        self.mid, self.unit = None, None

    def __str__(self):
        return f'{self.bigg} (reaction {self.rid})'


class MapModel:
    def __init__(self, path, settings=None, units_from=None):
        """`units_from`: the draft of this map (a MapModel or a path); its subsystems are the units of this map's
        reactions, which keep their ids. Without it the units are taken from this map's own pieces and titles."""
        self.path = path
        self.settings = dict(DEFAULT_SETTINGS, **(settings or {}))
        self.header, body = read_map(path)
        self.nodes, self.labels, self.canvas = body.get('nodes', {}), body.get('text_labels', {}), body.get('canvas', {})
        self._suffix = re.compile(self.settings['compartment_suffix'])
        self.currency = set(self.settings['currency'])
        self.carriers = set(self.settings['carriers'])
        self.reactions = {}
        self.users = defaultdict(list)  # metabolite node id -> reactions that use it
        for rid, data in body.get('reactions', {}).items():
            self.reactions[rid] = r = self._reaction(rid, data)
            for n in r.subs + r.prods:
                self.users[n].append(r)
        if isinstance(units_from, str):
            units_from = MapModel(units_from, settings)
        own = self._own_units()
        if units_from is not None:
            for rid, r in self.reactions.items():
                r.unit = units_from.reactions[rid].unit if rid in units_from.reactions else own.get(rid, 'map')
        else:
            for rid, r in self.reactions.items():
                r.unit = own.get(rid, 'map')
        self.units = defaultdict(list)
        for r in self.reactions.values():
            self.units[r.unit].append(r)

    # -- reading ---------------------------------------------------------------------------------------------------
    def _reaction(self, rid, data):
        r = Reaction(rid, data)
        sign = {m['bigg_id']: m['coefficient'] for m in data.get('metabolites', [])}
        seen = set()
        for seg in data.get('segments', {}).values():
            a, b = seg['from_node_id'], seg['to_node_id']
            if a not in self.nodes or b not in self.nodes:
                continue
            r.segments.append((a, b, seg.get('b1'), seg.get('b2')))
            for n in (a, b):
                if self.nodes[n]['node_type'] == 'midmarker':
                    r.mid = n
            for met, other, is_sub in ((a, b, True), (b, a, False)):
                if self.nodes[met]['node_type'] != 'metabolite' or self.nodes[other]['node_type'] == 'metabolite' or met in seen:
                    continue
                seen.add(met)
                coefficient = sign.get(self.nodes[met].get('bigg_id'))
                if coefficient is not None:  # the stoichiometry decides the side, the direction of the segment otherwise
                    is_sub = coefficient < 0
                (r.subs if is_sub else r.prods).append(met)
                r.marker[met] = other
        if r.mid is None:
            markers = [n for a, b, _, _ in r.segments for n in (a, b) if self.nodes[n]['node_type'] != 'metabolite']
            r.mid = markers[0] if markers else None
        r.primary_subs = [n for n in r.subs if self.nodes[n].get('node_is_primary')]
        r.primary_prods = [n for n in r.prods if self.nodes[n].get('node_is_primary')]
        r.sides = [n for n in r.subs + r.prods if not self.nodes[n].get('node_is_primary')]
        r.one_sided = not (r.subs and r.prods)
        return r

    def _own_units(self):
        """Reaction id -> unit name: the pieces of the map that hang together through primary nodes, each named by the
        title label nearest above it (pieces with the same title are one unit)."""
        parent = {}

        def find(x):
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for rid, r in self.reactions.items():
            find('r' + rid)
            for n in r.primary_subs + r.primary_prods:
                parent[find('n' + n)] = find('r' + rid)
        pieces = defaultdict(list)
        for rid in self.reactions:
            pieces[find('r' + rid)].append(rid)
        prefix = self.settings['exit_label_prefix']
        titles = [t for t in self.labels.values() if not str(t.get('text', '')).startswith(prefix)]
        out = {}
        for i, rids in enumerate(sorted(pieces.values(), key=lambda p: -len(p))):
            pts = [self.xy(n) for rid in rids for n in self.reactions[rid].subs + self.reactions[rid].prods]
            name = f'piece {i + 1}'
            if titles and pts:
                x0, y0, x1, y1 = bbox(pts)
                above = [t for t in titles if t['y'] <= (y0 + y1) / 2] or titles
                name = min(above, key=lambda t: math.hypot(t['x'] - (x0 + x1) / 2, t['y'] - y0))['text'].strip()
            elif len(pieces) == 1:
                name = 'map'
            for rid in rids:
                out[rid] = name
        return out

    # -- metabolites -----------------------------------------------------------------------------------------------
    def bigg(self, n):
        return self.nodes[n].get('bigg_id', '')

    def base(self, bigg_id):
        """A metabolite id without its compartment: pyr_c -> pyr."""
        return self._suffix.sub('', bigg_id)

    def compartment(self, bigg_id):
        m = self._suffix.search(bigg_id)
        return m.group(0)[1:] if m else ''

    def is_currency(self, n):
        return self.base(self.bigg(n)) in self.currency

    def xy(self, n):
        return (self.nodes[n]['x'], self.nodes[n]['y'])

    def radius(self, n):
        return self.settings['primary_radius' if self.nodes[n].get('node_is_primary') else 'secondary_radius']

    def metabolite_nodes(self):
        return [n for n in self.users if self.users[n]]

    def tag(self, n):
        """A metabolite node for a message: pyr_c (node 27)."""
        return f'{self.bigg(n)} (node {n})'

    # -- the backbone ----------------------------------------------------------------------------------------------
    def main(self, r):
        """(substrate nodes, product nodes) of the reaction's backbone: its primary nodes that are no currency
        metabolites; the primary nodes themselves where that leaves a side empty."""
        subs = [n for n in r.primary_subs if not self.is_currency(n)] or r.primary_subs
        prods = [n for n in r.primary_prods if not self.is_currency(n)] or r.primary_prods
        return subs, prods

    def backbone_edges(self, rs):
        """(metabolite id, metabolite id, reaction) for every substrate-product pair on the backbone of reactions `rs`."""
        out = []
        for r in rs:
            subs, prods = self.main(r)
            out += [(self.bigg(a), self.bigg(b), r) for a in subs for b in prods if self.bigg(a) != self.bigg(b)]
        return out

    def is_boundary(self, r):
        """An exchange, sink, demand or biomass reaction."""
        s = self.settings
        return r.one_sided or r.bigg.startswith(tuple(s['exchange_prefixes'])) or self.is_biomass(r)

    def is_biomass(self, r):
        return bool(re.search(self.settings['biomass_pattern'], r.bigg + ' ' + r.name))

    def unit_box(self, unit):
        pts = [self.xy(n) for r in self.units[unit] for n in r.primary_subs + r.primary_prods]
        pts = pts or [self.xy(n) for r in self.units[unit] for n in r.subs + r.prods]
        return bbox(pts) if pts else None

    def extent(self):
        pts = [self.xy(n) for n in self.metabolite_nodes()]
        return bbox(pts) if pts else (0, 0, 1, 1)


# -- geometry ----------------------------------------------------------------------------------------------------------
def bbox(pts):
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def dist(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def point_segment(p, a, b):
    """Distance of point p from the segment a-b."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    d2 = dx * dx + dy * dy
    t = 0.0 if d2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / d2))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


def point_line(p, a, b):
    """Distance of point p from the line through a and b."""
    d = dist(a, b)
    return dist(p, a) if d == 0 else abs(cross(a, b, p)) / d


def box_gap(a, b):
    """The gap between two boxes (0 if they touch or overlap)."""
    dx = max(a[0] - b[2], b[0] - a[2], 0)
    dy = max(a[1] - b[3], b[1] - a[3], 0)
    return math.hypot(dx, dy)


def curve(p, q, b1, b2, pieces=8):
    """The points of a segment's line: straight, or its Bezier curve sampled."""
    if not b1 or not b2:
        return [p, q]
    c1, c2 = (b1['x'], b1['y']), (b2['x'], b2['y'])
    out = []
    for i in range(pieces + 1):
        t = i / pieces
        u = 1 - t
        out.append(tuple(u ** 3 * p[k] + 3 * u * u * t * c1[k] + 3 * u * t * t * c2[k] + t ** 3 * q[k] for k in (0, 1)))
    return out


def roundness(pts):
    """How round the polygon through `pts` is: 4 pi area / perimeter^2 (a circle 1, a square 0.79, a line 0)."""
    n = len(pts)
    area = abs(sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1] for i in range(n))) / 2
    perimeter = sum(dist(pts[i], pts[(i + 1) % n]) for i in range(n))
    return 0.0 if perimeter == 0 else 4 * math.pi * area / perimeter ** 2


def is_convex(pts, tolerance=0.05):
    """Whether the polygon turns one way at every corner (corners that are nearly straight count as either)."""
    n = len(pts)
    signs = set()
    for i in range(n):
        o, a, b = pts[i], pts[(i + 1) % n], pts[(i + 2) % n]
        c = cross(o, a, b)
        scale = dist(o, a) * dist(a, b)
        if scale and abs(c) / scale > tolerance:
            signs.add(c > 0)
    return len(signs) <= 1
