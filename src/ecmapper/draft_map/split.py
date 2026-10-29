"""Splitting a draft map into a central map and auxiliary maps, and an overview of the maps."""

from __future__ import annotations

import copy
import html
import math
import re
from collections import Counter, defaultdict

import numpy as np

from .core import DraftMap, metabolite_key, text_width
from .layout import Ports, SequencePair, anneal_arrangement

CENTRAL = "Central metabolism"

#: subsystems of the BRENDA map that make up the central map by default
DEFAULT_CENTRAL = [
    "glycolysis",
    "gluconeogenesis",
    "pentose-phosphate-pathway",
    "Entner-Doudoroff-pathway",
    "oxidative-decarboxylation-of-pyruvate",
    "citric-acid-cycle",
    "oxidative-phosphorylation",
    "acetate-fermentation",
    "ethanol-fermentation",
    "lactate-fermentation",
]

#: names of the auxiliary maps: the colour groups of the BRENDA map, named after the subsystems in them
GROUP_NAMES = {
    "pw1": "Energy and C1 metabolism",
    "pw2": "Lipids",
    "pw3": "Amino acids",
    "pw4": "Nucleotides and cofactors",
    "pw5": "Carbohydrates",
    "pw6": "Fermentation and degradation",
    "pw7": "Other pathways",
    "pw0": "Other pathways",
}


def slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def reactions_of(dm: DraftMap, subsystems) -> set[str]:
    subsystems = set(subsystems)
    return {r for n in dm.nodes.values() if n.subsystem in subsystems for r in n.reactions}


def assign_maps(dm: DraftMap, central=DEFAULT_CENTRAL, explicit: dict[str, list[str]] | None = None,
                min_reactions: int = 25) -> tuple[list[str], dict[str, str]]:
    """The map of each subsystem: map names in order and {subsystem ID: map name}.

    Subsystems listed in `explicit` ({map name: [subsystem IDs]}) go to their
    map, those in `central` to the central map, all others to the map of
    their colour group. Colour group maps with fewer than `min_reactions`
    model reactions are merged into the map they have most inter-pathway
    links with (ties and no links: the largest map).
    """
    explicit = explicit or {}
    assignment = {}
    for name, subsystems in explicit.items():
        for s in subsystems:
            if s in dm.subsystems:
                assignment.setdefault(s, name)
    for s in central:
        if s in dm.subsystems:
            assignment.setdefault(s, CENTRAL)
    groups = list(dict.fromkeys(GROUP_NAMES.values()))
    for k, s in dm.subsystems.items():
        assignment.setdefault(k, GROUP_NAMES.get(s.color_group, groups[-1]))
    fixed = set(explicit) | {CENTRAL}
    while True:
        members = defaultdict(list)
        for k, name in assignment.items():
            members[name].append(k)
        size = {name: len(reactions_of(dm, subs)) for name, subs in members.items()}
        small = sorted((name for name in members if name not in fixed and size[name] < min_reactions),
                       key=lambda name: (size[name], name))
        if not small or len(members) == 1:
            break
        name = small[0]
        links = Counter()
        for l in dm.links:
            ma, mb = assignment[dm.nodes[l.a].subsystem], assignment[dm.nodes[l.b].subsystem]
            if ma != mb and name in (ma, mb):
                links[mb if ma == name else ma] += 1
        target = max((m for m in members if m != name), key=lambda m: (links[m], size[m], m))
        dm.log.append(f"split: map {name!r} ({size[name]} reactions) merged into {target!r}")
        for k in members[name]:
            assignment[k] = target
    order = [name for name in [*explicit, CENTRAL, *groups] if name in assignment.values()]
    return order, assignment


def split_map(dm: DraftMap, order: list[str], assignment: dict[str, str]) -> list[tuple[str, DraftMap]]:
    """One draft map per map name with the subsystems assigned to it.

    Inter-pathway links between maps are dropped; their end nodes get
    `stubs`, the other maps on which the metabolite continues. Needs edges
    inside subsystems (:func:`.core.localise_edges`).
    """
    map_of = {n.id: assignment[n.subsystem] for n in dm.nodes.values()}
    rank = {name: i for i, name in enumerate(order)}
    stubs = defaultdict(set)
    for l in dm.links:
        if map_of[l.a] != map_of[l.b]:
            stubs[l.a].add(map_of[l.b])
            stubs[l.b].add(map_of[l.a])
    result = []
    for name in order:
        sub = DraftMap()
        sub.subsystems = {k: copy.deepcopy(s) for k, s in dm.subsystems.items() if assignment[k] == name}
        sub.nodes = {i: copy.deepcopy(n) for i, n in dm.nodes.items() if map_of[i] == name}
        for e in dm.edges:
            if (e.enzyme in sub.nodes) != (e.met in sub.nodes):
                raise ValueError(f"edge {e.enzyme} - {e.met} between maps: localise the edges first")
        sub.edges = [copy.copy(e) for e in dm.edges if e.enzyme in sub.nodes]
        sub.links = [copy.copy(l) for l in dm.links if l.a in sub.nodes and l.b in sub.nodes]
        for i, n in sub.nodes.items():
            if stubs.get(i):
                n.stubs = sorted(stubs[i], key=rank.get)
        sub.orig_centres = {k: dm.orig_centres[k] for k in sub.subsystems}
        sub.log = list(dm.log)
        result.append((name, sub))
    n_stubs = sum(1 for v in stubs.values() if v)
    dm.log.append(f"split: {len(order)} maps; {sum(map_of[l.a] != map_of[l.b] for l in dm.links)} inter-pathway links "
                  f"between maps replaced by references on {n_stubs} metabolite nodes")
    return result


def cross_map_metabolites(dm: DraftMap, assignment: dict[str, str]) -> list[dict]:
    """Metabolites linked across maps: name, model metabolites and the maps."""
    maps, names, mets = defaultdict(set), {}, defaultdict(set)
    for l in dm.links:
        a, b = dm.nodes[l.a], dm.nodes[l.b]
        ma, mb = assignment[a.subsystem], assignment[b.subsystem]
        if ma != mb:
            k = metabolite_key(a)
            maps[k] |= {ma, mb}
            names.setdefault(k, getattr(a, "base_label", a.label))
            mets[k] |= a.model_mets | b.model_mets
    rows = [dict(name=names[k], model_metabolites=sorted(mets[k]), maps=sorted(v)) for k, v in maps.items()]
    rows.sort(key=lambda r: (-len(r["maps"]), str(r["name"])))
    return rows


# --- overview ---------------------------------------------------------------------------------------

BOX_H = 46.0
BOX_GAP = 24.0
REGION_PAD = 30.0
REGION_TITLE = 50.0
OVERVIEW_STYLE = (
    "text{font-family:Liberation Sans, Arial, Helvetica, sans-serif;text-anchor:middle;}"
    "rect.region{fill:rgb(248,248,248);stroke:rgb(150,150,150);stroke-width:2;}"
    "text.region{font-size:22px;font-weight:bold;fill:rgb(40,40,40);}"
    "rect.sub{stroke:rgb(120,120,120);}text.sub{font-size:14px;}text.count{font-size:11px;fill:rgb(80,80,80);}"
    "path.inside{stroke:rgb(150,150,150);fill:none;opacity:0.7;}"
    "path.between{stroke:rgb(150,60,170);fill:none;opacity:0.6;}"
)


def overview(maps: list[tuple[str, DraftMap]], dm: DraftMap, assignment: dict[str, str], gap: float = 120.0,
             iterations: int = 20000, seed: int = 1, hubs=frozenset()) -> dict:
    """Layout of an overview: each subsystem collapsed to a small box, one region per map.

    `maps` are the arranged maps (absolute positions). Inside a region the
    collapsed boxes keep the relative order of the subsystems on the map
    (sequence pair of their centres, packed). The regions are arranged with
    few crossings of the links between maps (annealing as for the
    subsystems). `dm` is the map before splitting, with all links; the
    links of the metabolites in `hubs` (keys of :func:`.core.hub_table`)
    are left out.
    """
    regions = []
    boxes = {}  # subsystem -> [x, y, w, h] relative to its region
    for name, m in maps:
        keys = list(m.subsystems)
        w = np.array([text_width(m.subsystems[k].label, 14) + 16 for k in keys])
        h = np.full(len(keys), BOX_H)
        cx = np.array([m.subsystems[k].x + m.subsystems[k].w / 2 for k in keys])
        cy = np.array([m.subsystems[k].y + m.subsystems[k].h / 2 for k in keys])
        x, y = SequencePair(w, h, BOX_GAP).pack(*SequencePair.from_positions(cx, cy))
        title_w = text_width(name, 22) + 2 * REGION_PAD
        content_w = float((x + w).max())
        dx = REGION_PAD + max(0.0, (title_w - 2 * REGION_PAD - content_w) / 2)
        for k, bx, by, bw in zip(keys, x, y, w):
            boxes[k] = [float(bx) + dx, float(by) + REGION_TITLE, float(bw), BOX_H]
        rw = max(content_w + 2 * REGION_PAD, title_w)
        rh = float((y + h).max()) + REGION_TITLE + REGION_PAD
        regions.append(dict(name=name, w=rw, h=rh))
    index = {r["name"]: i for i, r in enumerate(regions)}
    pairs = Counter()
    for l in dm.links:
        sa, sb = dm.nodes[l.a].subsystem, dm.nodes[l.b].subsystem
        if sa != sb and metabolite_key(dm.nodes[l.a]) not in hubs:
            pairs[tuple(sorted((sa, sb)))] += 1
    cross = [(a, b, c) for (a, b), c in sorted(pairs.items()) if assignment[a] != assignment[b]]
    rw = np.array([r["w"] for r in regions])
    rh = np.array([r["h"] for r in regions])
    if len(regions) > 2:
        centre = lambda k: (boxes[k][0] + boxes[k][2] / 2, boxes[k][1] + BOX_H / 2)  # noqa: E731
        ports = Ports(
            [r["name"] for r in regions], rw, rh,
            np.array([index[assignment[a]] for a, _, _ in cross], int),
            np.array([index[assignment[b]] for _, b, _ in cross], int),
            np.array([centre(a) for a, _, _ in cross], float).reshape(-1, 2),
            np.array([centre(b) for _, b, _ in cross], float).reshape(-1, 2),
            [(a, b) for a, b, _ in cross], rh,
        )
        x0, y0 = SequencePair(rw, rh, gap).pack(np.arange(len(regions)), np.arange(len(regions)))
        x, y, *_ = anneal_arrangement(ports, x0, y0, gap=gap, iterations=iterations, seed=seed)
    else:
        x, y = SequencePair(rw, rh, gap).pack(np.arange(len(regions)), np.arange(len(regions)))
    for r, rx, ry in zip(regions, x, y):
        r["x"], r["y"] = float(rx), float(ry)
    for k, b in boxes.items():
        r = regions[index[assignment[k]]]
        b[0] += r["x"]
        b[1] += r["y"]
    inside = [(a, b, c) for (a, b), c in sorted(pairs.items()) if assignment[a] == assignment[b]]
    return dict(regions=regions, boxes=boxes, links_inside=inside, links_between=cross)


def write_overview_svg(layout: dict, maps: list[tuple[str, DraftMap]], path: str, links: dict[str, str],
                       margin: float = 40.0) -> None:
    """SVG of the overview; `links` maps each map name to the file shown when clicking its title."""
    subsystems = {k: s for _, m in maps for k, s in m.subsystems.items()}
    counts = {k: len({r for n in m.members(k) for r in n.reactions}) for _, m in maps for k in m.subsystems}
    regions, boxes = layout["regions"], layout["boxes"]
    x0 = min(r["x"] for r in regions) - margin
    y0 = min(r["y"] for r in regions) - margin
    W = max(r["x"] + r["w"] for r in regions) - x0 + margin
    H = max(r["y"] + r["h"] for r in regions) - y0 + margin
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="{x0:.1f} {y0:.1f} {W:.1f} {H:.1f}" '
        f'width="{W:.0f}" height="{H:.0f}">',
        f"<style>{OVERVIEW_STYLE}</style>",
        f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{W:.1f}" height="{H:.1f}" style="fill:white"/>',
    ]
    for r in regions:
        title = html.escape(r["name"])
        href = links.get(r["name"])
        text = f'<a href="{html.escape(href)}">{title}</a>' if href else title
        out.append(f'<rect class="region" rx="12" x="{r["x"]:.1f}" y="{r["y"]:.1f}" width="{r["w"]:.1f}" '
                   f'height="{r["h"]:.1f}"/><text class="region" x="{r["x"] + r["w"] / 2:.1f}" '
                   f'y="{r["y"] + 32:.1f}">{text}</text>')

    def centre(k):
        b = boxes[k]
        return b[0] + b[2] / 2, b[1] + b[3] / 2

    for css, items in (("inside", layout["links_inside"]), ("between", layout["links_between"])):
        for a, b, c in items:
            (ax, ay), (bx, by) = centre(a), centre(b)
            width = 1.5 + 1.5 * math.log2(c)
            out.append(f'<path class="{css}" style="stroke-width:{width:.1f}" '
                       f'd="M {ax:.1f} {ay:.1f} L {bx:.1f} {by:.1f}"><title>{html.escape(subsystems[a].label)} – '
                       f'{html.escape(subsystems[b].label)}: {c} links</title></path>')
    for k, (bx, by, bw, bh) in boxes.items():
        s = subsystems[k]
        out.append(
            f'<g><rect class="sub" rx="8" x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" height="{bh:.1f}" '
            f'style="fill:{s.fill};fill-opacity:0.85"/><text class="sub" x="{bx + bw / 2:.1f}" y="{by + 19:.1f}">'
            f'{html.escape(s.label)}</text><text class="count" x="{bx + bw / 2:.1f}" y="{by + 36:.1f}">'
            f'{counts[k]} reactions</text></g>'
        )
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))
