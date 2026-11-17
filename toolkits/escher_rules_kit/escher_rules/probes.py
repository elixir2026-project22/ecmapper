"""Probes: what the kit can find out about a map for a rule.

A probe looks at a map and returns a Finding: `applies` (the places of this map where the rule's trigger holds, for
the brief) and `violations` (the places where the rule is broken, for the check). A probe that can only find the
places returns `checked=False`: the rule stays one to judge by eye. A rule is bound to a probe, with the probe's
parameters, in rules/NAME.checks.yaml; `python run.py probes` lists the probes and their parameters.
"""
import math
import re
from collections import defaultdict

from . import mapmodel as G

PROBES = {}


class Finding:
    def __init__(self, applies=None, violations=None, checked=True, summary=''):
        self.applies = list(applies or [])
        self.violations = list(violations or [])
        self.checked = checked  # False: the probe only finds where the rule applies
        self.summary = summary  # one line with the measured values


class Probe:
    def __init__(self, fn, defaults, needs_previous=False):
        self.fn, self.name, self.defaults, self.needs_previous = fn, fn.__name__, defaults, needs_previous
        self.doc = ' '.join((fn.__doc__ or '').split())


def probe(needs_previous=False, **defaults):
    def register(fn):
        PROBES[fn.__name__] = Probe(fn, defaults, needs_previous)
        return fn
    return register


def run(name, model, params=None, previous=None):
    p = PROBES[name]
    kw = dict(p.defaults, **(params or {}))
    if p.needs_previous:
        return p.fn(model, previous, **kw)
    return p.fn(model, **kw)


def _names(model, r):
    """A reaction as text: PFK (reaction 140): f6p_c -> fdp_c."""
    subs, prods = model.main(r)
    arrow = ' <-> ' if r.reversible else ' -> '
    return f"{r}: {' + '.join(model.bigg(n) for n in subs)}{arrow}{' + '.join(model.bigg(n) for n in prods)}"


def _pathway_reactions(model):
    return [r for r in model.reactions.values() if not model.is_boundary(r)]


# == motif: one reaction ================================================================================================
@probe()
def one_primary_pair(model):
    """Every reaction has exactly one primary substrate node and one primary product node (a one-sided reaction:
    one primary node)."""
    bad = []
    for r in model.reactions.values():
        ns, np_ = len(r.primary_subs), len(r.primary_prods)
        if r.one_sided:
            if ns + np_ != 1:
                bad.append(f'{r}: {ns + np_} primary nodes, one is expected')
        elif (ns, np_) != (1, 1):
            names = ', '.join(model.bigg(n) for n in r.primary_subs + r.primary_prods) or 'none'
            bad.append(f'{r}: {ns} primary substrates and {np_} primary products ({names})')
    return Finding([], bad, summary=f'{len(model.reactions)} reactions')


@probe()
def primary_not_carrier(model):
    """No carrier (setting `carriers`) is drawn as a primary node of a reaction that has another substrate and
    another product to carry the backbone."""
    bad = []
    for r in _pathway_reactions(model):
        free = lambda nodes: [n for n in nodes if model.base(model.bigg(n)) not in model.carriers | model.currency]
        if not (free(r.subs) and free(r.prods)):
            continue  # the reaction converts carriers: nothing else can be primary
        carried = [n for n in r.primary_subs + r.primary_prods if model.base(model.bigg(n)) in model.carriers]
        if carried:
            bad.append(f"{r}: the carrier {', '.join(model.tag(n) for n in carried)} is primary; "
                       f"carrier-free pair: {model.bigg(free(r.subs)[0])} -> {model.bigg(free(r.prods)[0])}")
    return Finding([], bad)


@probe(max_length=500, min_backbone_angle=90)
def side_branches(model, max_length, min_backbone_angle):
    """Side participants hang at most `max_length` map units from their reaction marker, and the backbone does not
    bend at the reaction by more than 180 - `min_backbone_angle` degrees."""
    bad, n_sides = [], 0
    for r in model.reactions.values():
        for n in r.sides:
            n_sides += 1
            d = G.dist(model.xy(n), model.xy(r.marker[n]))
            if d > max_length:
                bad.append(f'{r}: {model.tag(n)} hangs {d:.0f} units from the reaction (at most {max_length})')
        if len(r.primary_subs) == 1 and len(r.primary_prods) == 1 and r.mid:
            a, m, b = model.xy(r.primary_subs[0]), model.xy(r.mid), model.xy(r.primary_prods[0])
            if G.dist(a, m) and G.dist(b, m):
                cos = ((a[0] - m[0]) * (b[0] - m[0]) + (a[1] - m[1]) * (b[1] - m[1])) / (G.dist(a, m) * G.dist(b, m))
                angle = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
                if angle < min_backbone_angle:
                    bad.append(f'{r}: the backbone bends at the reaction ({angle:.0f} degrees between its two halves, '
                               f'at least {min_backbone_angle})')
    return Finding([], bad, summary=f'{n_sides} side participants')


@probe()
def cofactor_pairs_same_side(model):
    """Both members of a cofactor pair (setting `cofactor_pairs`) of a reaction lie on the same side of its
    backbone."""
    pairs = [tuple(p) for p in model.settings['cofactor_pairs']]
    applies, bad = [], []
    for r in model.reactions.values():
        if len(r.primary_subs) != 1 or len(r.primary_prods) != 1:
            continue
        a, b = model.xy(r.primary_subs[0]), model.xy(r.primary_prods[0])
        side = {model.base(model.bigg(n)): n for n in r.sides}
        for x, y in pairs:
            if x in side and y in side and (side[x] in r.subs) != (side[y] in r.subs):
                applies.append(f'{r}: {x} / {y}')
                cx, cy = G.cross(a, b, model.xy(side[x])), G.cross(a, b, model.xy(side[y]))
                if cx * cy < 0:
                    bad.append(f'{r}: {model.tag(side[x])} and {model.tag(side[y])} are on opposite sides of the backbone')
    return Finding(applies, bad)


def _antiparallel(model):
    """Pairs of reactions that convert the same backbone metabolites in opposite directions."""
    by_pair = defaultdict(list)
    for r in _pathway_reactions(model):
        subs, prods = model.main(r)
        key = (frozenset(model.bigg(n) for n in subs), frozenset(model.bigg(n) for n in prods))
        if key[0] and key[1]:
            by_pair[key].append(r)
    out = []
    for (subs, prods), forward in by_pair.items():
        back = by_pair.get((prods, subs), [])
        if back and sorted(subs) < sorted(prods):
            out += [(f, b) for f in forward for b in back if f.unit == b.unit]
    return out


@probe(max_offset=350)
def antiparallel_one_axis(model, max_offset):
    """Two antiparallel reactions between the same metabolites (PFK / FBPase) join the same two nodes, and neither
    bulges more than `max_offset` map units away from the line between them."""
    applies, bad = [], []
    for f, b in _antiparallel(model):
        applies.append(f'{_names(model, f)}  against  {b}  (unit {f.unit})')
        fs, fp = model.main(f)
        bs, bp = model.main(b)
        if set(fs) != set(bp) or set(fp) != set(bs):
            bad.append(f'{f} and {b} are drawn between different copies of their metabolites: let both use the same two nodes')
            continue
        p, q = model.xy(fs[0]), model.xy(fp[0])
        for r in (f, b):
            off = G.point_line(model.xy(r.mid), p, q) if r.mid else 0
            if off > max_offset:
                bad.append(f'{r} bulges {off:.0f} units from the axis {model.bigg(fs[0])} - {model.bigg(fp[0])} '
                           f'(at most {max_offset}): a two-node cycle')
    return Finding(applies, bad)


# == motif: the shape of a unit =========================================================================================
def _cycles_of(edges, max_length, budget=200000):
    """The simple cycles of 3 to `max_length` metabolites in the graph of `edges` (a, b, reaction), as lists of
    metabolite ids."""
    adj = defaultdict(set)
    for a, b, _ in edges:
        adj[a].add(b)
        adj[b].add(a)
    found, steps = [], [0]

    def walk(start, path, seen):
        steps[0] += 1
        if steps[0] > budget:
            return
        for nxt in sorted(adj[path[-1]]):
            if nxt == start and len(path) >= 3 and path[1] < path[-1]:
                found.append(list(path))
            elif nxt > start and nxt not in seen and len(path) < max_length:
                seen.add(nxt)
                walk(start, path + [nxt], seen)
                seen.discard(nxt)

    for start in sorted(adj):
        walk(start, [start], {start})
    return found


def _directed(cycle, edges):
    """Whether the reactions allow going round the cycle in one direction."""
    ways = defaultdict(set)
    for a, b, r in edges:
        ways[(a, b)].add(True)
        if r.reversible:
            ways[(b, a)].add(True)
    n = len(cycle)
    forward = all(ways.get((cycle[i], cycle[(i + 1) % n])) for i in range(n))
    backward = all(ways.get((cycle[(i + 1) % n], cycle[i])) for i in range(n))
    return forward or backward


def unit_cycles(model, max_length=12, max_shared=2):
    """unit -> its rings: the cycles of the unit's backbone, the longest that can be gone round first; a cycle that
    shares more than `max_shared` metabolites with a ring already taken is a bypass of it, not a ring of its own.
    Cycles are searched inside a unit only (rule M11)."""
    out = {}
    for unit, rs in model.units.items():
        edges = [e for e in model.backbone_edges([r for r in rs if not model.is_boundary(r)])
                 if not (model.base(e[0]) in model.currency or model.base(e[1]) in model.currency)]
        cycles = _cycles_of(edges, max_length)
        cycles.sort(key=lambda c: (not _directed(c, edges), -len(c), c))
        rings = []
        for c in cycles:
            if all(len(set(c) & set(k)) <= max_shared for k in rings):
                rings.append(c)
        if rings:
            out[unit] = (rings, edges)
    return out


def _ring_nodes(model, cycle, edges):
    """The node of every metabolite of a cycle as its reactions draw it, or (None, metabolite) where two reactions
    that follow each other on the cycle use different copies of a metabolite."""
    n = len(cycle)
    chosen = {}
    for i, met in enumerate(cycle):
        prev_, next_ = cycle[i - 1], cycle[(i + 1) % n]
        used = []
        for other in (prev_, next_):
            nodes = set()
            for a, b, r in edges:
                if {a, b} == {met, other}:
                    nodes |= {x for x in r.subs + r.prods if model.bigg(x) == met and model.nodes[x].get('node_is_primary')}
            used.append(nodes)
        common = used[0] & used[1]
        if not common:
            return None, met
        chosen[met] = sorted(common)[0]
    return [chosen[m] for m in cycle], None


@probe(min_length=3, max_length=12, max_shared=2, min_roundness=0.5)
def cycles_as_rings(model, min_length, max_length, max_shared, min_roundness):
    """A cycle inside a unit (`min_length` to `max_length` metabolites) is drawn as a ring: closed, convex, and round
    (4 pi area / perimeter^2 of the polygon of its metabolites at least `min_roundness`). Cycles sharing more than
    `max_shared` metabolites with a longer one count as its bypass."""
    applies, bad = [], []
    for unit, (rings, edges) in unit_cycles(model, max_length, max_shared).items():
        for c in rings:
            if len(c) < min_length:
                continue
            text = f"unit {unit}: {' - '.join(c)} ({len(c)} metabolites)"
            applies.append(text)
            nodes, broken = _ring_nodes(model, c, edges)
            if nodes is None:
                bad.append(f'{text}: not closed, {broken} is drawn as separate copies where the cycle passes it')
                continue
            pts = [model.xy(x) for x in nodes]
            rnd = G.roundness(pts)
            if not G.is_convex(pts) or rnd < min_roundness:
                bad.append(f'{text}: not a ring (roundness {rnd:.2f}, at least {min_roundness}'
                           + ('' if G.is_convex(pts) else '; its outline turns inward') + ')')
    return Finding(applies, bad)


@probe(max_length=12)
def cycles_sharing_metabolites(model, max_length):
    """Finds pairs of rings (of one unit or of two) that share one or two metabolites: candidates for two adjacent
    rings joined at the shared element."""
    rings = [(unit, c) for unit, (rs, _) in unit_cycles(model, max_length).items() for c in rs]
    applies = []
    for i, (u1, c1) in enumerate(rings):
        for u2, c2 in rings[i + 1:]:
            shared = sorted(set(c1) & set(c2))
            if 1 <= len(shared) <= 2:
                applies.append(f"{u1}: {' - '.join(c1)}  and  {u2}: {' - '.join(c2)}  share {', '.join(shared)}")
    return Finding(applies, checked=False)


@probe(min_branches=3)
def branch_points(model, min_branches):
    """Finds the metabolites of a unit's backbone with at least `min_branches` neighbours: where a trunk ends in a
    fan."""
    applies = []
    for unit, rs in model.units.items():
        nb = defaultdict(set)
        for a, b, _ in model.backbone_edges([r for r in rs if not model.is_boundary(r)]):
            if model.base(a) not in model.currency and model.base(b) not in model.currency:
                nb[a].add(b)
                nb[b].add(a)
        for met in sorted(nb, key=lambda m: (-len(nb[m]), m)):
            if len(nb[met]) >= min_branches:
                applies.append(f"unit {unit}: {met} branches to {', '.join(sorted(nb[met]))}")
    return Finding(applies, checked=False)


@probe(min_units=2)
def homologous_units(model, min_units):
    """Finds metabolites that exist in several compartments on this map (one pathway present in two compartments:
    candidates for parallel limbs)."""
    comps = defaultdict(set)
    for n in model.metabolite_nodes():
        if not model.is_currency(n):
            comps[model.base(model.bigg(n))].add(model.compartment(model.bigg(n)))
    multi = sorted(b for b, cs in comps.items() if len(cs) >= min_units)
    applies = [f"{len(multi)} metabolites are drawn in several compartments: {', '.join(multi[:25])}"
               + (' ...' if len(multi) > 25 else '')] if multi else []
    return Finding(applies, checked=False)


# == motif: orientation and geometry ====================================================================================
@probe(direction='down', min_fraction=0.6)
def flow_direction(model, direction, min_fraction):
    """At least `min_fraction` of the irreversible reactions outside rings run in `direction` (down, up, right,
    left), measured from the primary substrate to the primary product."""
    vec = {'down': (0, 1), 'up': (0, -1), 'right': (1, 0), 'left': (-1, 0)}[direction]
    in_ring = {m for rings, _ in unit_cycles(model).values() for c in rings for m in c}
    against, total = [], 0
    for r in _pathway_reactions(model):
        subs, prods = model.main(r)
        if r.reversible or len(subs) != 1 or len(prods) != 1 or model.is_currency(subs[0]) or model.is_currency(prods[0]):
            continue
        if model.bigg(subs[0]) in in_ring and model.bigg(prods[0]) in in_ring:
            continue
        a, b = model.xy(subs[0]), model.xy(prods[0])
        total += 1
        if (b[0] - a[0]) * vec[0] + (b[1] - a[1]) * vec[1] <= 0:
            against.append(r)
    if not total:
        return Finding()
    fraction = 1 - len(against) / total
    summary = f'{fraction:.0%} of {total} irreversible reactions outside rings run {direction}'
    bad = []
    if fraction < min_fraction:
        bad = [f"{summary} (at least {min_fraction:.0%}); against it: {', '.join(str(r) for r in against[:12])}"
               + (' ...' if len(against) > 12 else '')]
    return Finding([], bad, summary=summary)


@probe(tolerance_deg=10, min_fraction=0.5)
def axis_aligned(model, tolerance_deg, min_fraction):
    """At least `min_fraction` of the backbone connections outside rings are horizontal or vertical within
    `tolerance_deg` degrees."""
    in_ring = {m for rings, _ in unit_cycles(model).values() for c in rings for m in c}
    aligned = total = 0
    for r in _pathway_reactions(model):
        subs, prods = model.main(r)
        if len(subs) != 1 or len(prods) != 1 or model.is_currency(subs[0]) or model.is_currency(prods[0]):
            continue
        if model.bigg(subs[0]) in in_ring and model.bigg(prods[0]) in in_ring:
            continue
        a, b = model.xy(subs[0]), model.xy(prods[0])
        if a == b:
            continue
        total += 1
        angle = math.degrees(math.atan2(abs(b[1] - a[1]), abs(b[0] - a[0])))
        aligned += angle <= tolerance_deg or angle >= 90 - tolerance_deg
    if not total:
        return Finding()
    fraction = aligned / total
    summary = f'{fraction:.0%} of {total} backbone connections outside rings are axis-aligned'
    return Finding([], [f'{summary} (at least {min_fraction:.0%})'] if fraction < min_fraction else [], summary=summary)


# == arrangement ========================================================================================================
@probe(min_reactions=6, max_reactions=60, max_map=200)
def unit_size(model, min_reactions, max_reactions, max_map):
    """A unit has at most `max_reactions` reactions and the map at most `max_map`; units below `min_reactions` are
    listed (to be drawn with a neighbour), not reported as broken."""
    applies, bad = [], []
    sizes = sorted(((len(rs), u) for u, rs in model.units.items()), reverse=True)
    small = [f'{u} ({n})' for n, u in sizes if n < min_reactions]
    if small:
        applies.append(f"units below {min_reactions} reactions (draw them attached to a neighbour): {', '.join(small)}")
    for n, u in sizes:
        if n > max_reactions:
            bad.append(f'unit {u} has {n} reactions (at most {max_reactions}): split it')
    if len(model.reactions) > max_map:
        bad.append(f'the map has {len(model.reactions)} reactions (at most {max_map} per map): split it into several maps')
    summary = f"{len(model.reactions)} reactions in {len(sizes)} units: {', '.join(f'{u} {n}' for n, u in sizes)}"
    return Finding(applies, bad, summary=summary)


def _interfaces(model):
    """metabolite id -> the units it occurs in on a backbone (currency metabolites left out), for those in several."""
    units = defaultdict(set)
    for r in _pathway_reactions(model):
        for n in r.primary_subs + r.primary_prods:
            if not model.is_currency(n):
                units[model.bigg(n)].add(r.unit)
    return {m: us for m, us in units.items() if len(us) > 1}


def _links(model):
    """(unit, unit) -> the metabolites they share, most first."""
    links = defaultdict(list)
    for met, us in sorted(_interfaces(model).items()):
        us = sorted(us)
        for i, a in enumerate(us):
            for b in us[i + 1:]:
                links[(a, b)].append(met)
    return dict(sorted(links.items(), key=lambda kv: (-len(kv[1]), kv[0])))


@probe()
def interface_metabolites(model):
    """Finds the metabolites shared between units (the interface metabolites), each with its units and the number of
    nodes it is drawn with."""
    count = defaultdict(set)
    for n in model.metabolite_nodes():
        count[model.bigg(n)].add(n)
    applies = [f"{met}: {', '.join(sorted(us))} (drawn as {len(count[met])} node{'s' if len(count[met]) != 1 else ''})"
               for met, us in sorted(_interfaces(model).items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    return Finding(applies, checked=False)


@probe(min_shared=2, max_gap=1500)
def linked_units_adjacent(model, min_shared, max_gap):
    """Units that share at least `min_shared` interface metabolites are neighbours: the gap between their boxes is at
    most `max_gap` map units."""
    applies, bad = [], []
    for (a, b), mets in _links(model).items():
        applies.append(f"{a}  -  {b}: {len(mets)} shared ({', '.join(mets[:8])}{' ...' if len(mets) > 8 else ''})")
        if len(mets) >= min_shared:
            ba, bb = model.unit_box(a), model.unit_box(b)
            if ba and bb and G.box_gap(ba, bb) > max_gap:
                bad.append(f'{a} and {b} share {len(mets)} metabolites but are {G.box_gap(ba, bb):.0f} units apart (at most {max_gap})')
    return Finding(applies, bad)


@probe(edge_fraction=0.25)
def boundary_reactions(model, edge_fraction):
    """Exchange, sink and demand reactions lie in the outer `edge_fraction` of the map (0: not checked); the
    metabolite nodes of a biomass reaction are used by no other reaction (its precursors are clones)."""
    applies, bad = [], []
    x0, y0, x1, y1 = model.extent()
    w, h = max(x1 - x0, 1), max(y1 - y0, 1)
    for r in model.reactions.values():
        if not model.is_boundary(r):
            continue
        kind = 'biomass' if model.is_biomass(r) else 'exchange / sink / demand'
        applies.append(f'{r}: {kind}, unit {r.unit}')
        if model.is_biomass(r):
            shared = [n for n in r.subs + r.prods if len(model.users[n]) > 1]
            if shared:
                bad.append(f"{r}: biomass uses nodes of other reactions ({', '.join(model.tag(n) for n in shared[:6])}"
                           f"{' ...' if len(shared) > 6 else ''}): give it its own copies")
        elif edge_fraction and r.mid:
            x, y = model.xy(r.mid)
            if min(x - x0, x1 - x, (y - y0) * w / h, (y1 - y) * w / h) > edge_fraction * w:
                bad.append(f'{r} lies in the interior of the map: move it to the edge, next to the pathway it feeds or drains')
    return Finding(applies, bad)


def _unit_named(model, pattern):
    return [u for u in model.units if re.search(pattern, u, re.I)]


@probe(relations=[])
def relative_position(model, relations):
    """The units keep the relations listed in `relations`. Each is {unit: PATTERN, is: above | below | left_of |
    right_of, other: PATTERN} (centres compared) or {unit: PATTERN, is: near, metabolite: BASE_ID, max_distance: N}
    (the unit's box is at most N units from a node of that metabolite used outside the unit). PATTERN is a regular
    expression matched against the unit names, ignoring case; a relation whose units are not on the map is skipped."""
    applies, bad = [], []
    centre = lambda box: ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2)
    for rel in relations:
        kind = rel.get('is')
        for unit in _unit_named(model, rel.get('unit', '$^')):
            box = model.unit_box(unit)
            if box is None:
                continue
            if kind == 'near':
                met, limit = rel.get('metabolite', ''), rel.get('max_distance', 1500)
                outside = [n for n in model.metabolite_nodes() if model.base(model.bigg(n)) == met
                           and any(r.unit != unit for r in model.users[n])]
                if not outside:
                    continue
                applies.append(f'{unit} near {met}')
                d = min(G.box_gap(box, model.xy(n) * 2) for n in outside)
                if d > limit:
                    bad.append(f'{unit} is {d:.0f} units from the nearest {met} of another unit (at most {limit})')
                continue
            for other in _unit_named(model, rel.get('other', '$^')):
                obox = model.unit_box(other)
                if other == unit or obox is None:
                    continue
                applies.append(f"{unit} {str(kind).replace('_', ' ')} {other}")
                (x, y), (ox, oy) = centre(box), centre(obox)
                ok = {'above': y < oy, 'below': y > oy, 'left_of': x < ox, 'right_of': x > ox}.get(kind)
                if ok is None:
                    bad.append(f"relation `{kind}` is unknown (above, below, left_of, right_of, near)")
                elif not ok:
                    bad.append(f"{unit} is not {str(kind).replace('_', ' ')} {other}")
    return Finding(applies, bad)


@probe()
def compartments(model):
    """Finds the compartments of the map's metabolites (by id suffix)."""
    comps = defaultdict(int)
    for n in model.metabolite_nodes():
        comps[model.compartment(model.bigg(n)) or '(none)'] += 1
    if len(comps) < 2:
        return Finding(checked=False)
    return Finding([f"{len(comps)} compartments: " + ', '.join(f'{c} ({k} nodes)' for c, k in sorted(comps.items()))], checked=False)


@probe()
def transport_reactions(model):
    """Finds the reactions whose backbone crosses a compartment boundary."""
    applies = []
    for r in _pathway_reactions(model):
        cs = {model.compartment(model.bigg(n)) for n in r.subs + r.prods}
        if len(cs) > 1:
            applies.append(f"{r}: {' / '.join(sorted(c or '(none)' for c in cs))}, unit {r.unit}")
    return Finding(applies, checked=False)


@probe(needs_previous=True, tolerance=1.0)
def positions_stable(model, previous, tolerance):
    """No node that an earlier version of the map (`--previous MAP`) already placed has moved by more than
    `tolerance` map units."""
    if previous is None:
        return Finding(checked=False, summary='not checked: give the earlier map with --previous')
    moved = []
    for n, node in model.nodes.items():
        old = previous.nodes.get(n)
        if old and old.get('node_type') == node.get('node_type') and old.get('bigg_id') == node.get('bigg_id'):
            d = G.dist((old['x'], old['y']), (node['x'], node['y']))
            if d > tolerance:
                moved.append((d, n))
    moved.sort(reverse=True)
    bad = [f"{model.bigg(n) or model.nodes[n]['node_type']} (node {n}) moved {d:.0f} units" for d, n in moved]
    common = sum(1 for n in model.nodes if n in previous.nodes)
    return Finding([], bad, summary=f'{len(moved)} of {common} nodes of the earlier map moved')


# == encoding ===========================================================================================================
@probe(max_reactions=1)
def currency_not_shared(model, max_reactions):
    """A node of a currency metabolite (setting `currency`) is used by at most `max_reactions` reaction."""
    bad, total = [], 0
    for n in model.metabolite_nodes():
        if model.is_currency(n):
            total += 1
            if len(model.users[n]) > max_reactions:
                bad.append(f"{model.tag(n)} is shared by {', '.join(str(r) for r in model.users[n][:6])}"
                           f"{' ...' if len(model.users[n]) > 6 else ''}: one copy per reaction")
    return Finding([], bad, summary=f'{total} nodes of currency metabolites')


@probe(min_degree=8)
def currency_candidates(model, min_degree):
    """Finds metabolites in at least `min_degree` reactions that the setting `currency` does not list: decide for each
    whether it is currency here (add it to the setting) or a hub of the backbone (pyruvate, glutamate, acetyl-CoA)."""
    degree = defaultdict(set)
    for r in model.reactions.values():
        for n in r.subs + r.prods:
            degree[model.base(model.bigg(n))].add(r.rid)
    found = sorted(((len(rs), b) for b, rs in degree.items() if len(rs) >= min_degree and b not in model.currency), reverse=True)
    return Finding([f'{b}: in {k} reactions' for k, b in found], checked=False)


@probe(max_nodes=1, per_unit=False)
def drawn_once(model, max_nodes, per_unit):
    """A metabolite that is no currency metabolite is drawn as at most `max_nodes` node per compartment
    (`per_unit: true`: per unit, so that copies at unit boundaries are allowed)."""
    nodes = defaultdict(set)
    for n in model.metabolite_nodes():
        if model.is_currency(n):
            continue
        if per_unit:
            for r in model.users[n]:
                nodes[(model.bigg(n), r.unit)].add(n)
        else:
            nodes[(model.bigg(n), None)].add(n)
    bad = []
    for (met, unit), ns in sorted(nodes.items(), key=lambda kv: (-len(kv[1]), kv[0][0], str(kv[0][1]))):
        if len(ns) > max_nodes:
            units = sorted({r.unit for n in ns for r in model.users[n]})
            where = f'in unit {unit}' if per_unit else f"units: {', '.join(units)}"
            bad.append(f"{met} is drawn {len(ns)} times (nodes {', '.join(sorted(ns, key=_num))}; {where})")
    return Finding([], bad, summary=f'{len({m for m, _ in nodes})} metabolites that are no currency')


def _num(s):
    return (0, int(s)) if str(s).isdigit() else (1, str(s))


@probe(clearance=0)
def no_node_overlap(model, clearance):
    """No two metabolite nodes overlap: their centres are at least the sum of their radii (settings `primary_radius`,
    `secondary_radius`) plus `clearance` apart."""
    ns = sorted(model.metabolite_nodes(), key=lambda n: model.nodes[n]['x'])
    reach = 2 * max(model.settings['primary_radius'], model.settings['secondary_radius']) + clearance
    bad = []
    for i, a in enumerate(ns):
        ax, ay = model.xy(a)
        for b in ns[i + 1:]:
            bx, by = model.xy(b)
            if bx - ax > reach:
                break
            if math.hypot(bx - ax, by - ay) < model.radius(a) + model.radius(b) + clearance:
                bad.append(f'{model.tag(a)} and {model.tag(b)} overlap at ({ax:.0f}, {ay:.0f})')
    return Finding([], bad, summary=f'{len(ns)} metabolite nodes')


@probe(clearance=0)
def edges_avoid_nodes(model, clearance):
    """No line of a reaction passes through a metabolite node that it does not end at (closer to its centre than its
    radius plus `clearance`)."""
    ns = [(n, model.xy(n), model.radius(n) + clearance) for n in model.metabolite_nodes()]
    bad = []
    for r in model.reactions.values():
        hit, own = set(), set(r.subs + r.prods)
        for a, b, b1, b2 in r.segments:
            pts = G.curve(model.xy(a), model.xy(b), b1, b2)
            x0, y0, x1, y1 = G.bbox(pts)
            for n, (x, y), rad in ns:
                if n in own and min(G.dist((x, y), pts[0]), G.dist((x, y), pts[-1])) < rad:
                    continue  # the line starts at a marker drawn on a node of its own reaction
                if n in (a, b) or n in hit or not (x0 - rad <= x <= x1 + rad and y0 - rad <= y <= y1 + rad):
                    continue
                if any(G.point_segment((x, y), pts[i], pts[i + 1]) < rad for i in range(len(pts) - 1)):
                    hit.add(n)
        bad += [f'a line of {r} passes through {model.tag(n)} at ({model.xy(n)[0]:.0f}, {model.xy(n)[1]:.0f})' for n in sorted(hit, key=_num)]
    return Finding([], bad, summary=f'{sum(len(r.segments) for r in model.reactions.values())} lines')


@probe(target_width_mm=180, min_pt=5, label_units=20)
def label_legible(model, target_width_mm, min_pt, label_units):
    """With the map scaled to `target_width_mm`, a metabolite label (`label_units` map units high) is at least
    `min_pt` points high."""
    x0, _, x1, _ = model.extent()
    width = float(model.canvas.get('width') or (x1 - x0) or 1)
    pt = label_units * target_width_mm / width * 72 / 25.4
    need = min_pt * width * 25.4 / (72 * label_units)
    summary = f'at {target_width_mm} mm width the metabolite labels are {pt:.1f} pt high'
    bad = [] if pt >= min_pt else [f'{summary} (at least {min_pt} pt): legible from {need:.0f} mm width, or split the map']
    return Finding([], bad, summary=summary)
