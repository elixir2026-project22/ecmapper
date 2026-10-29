"""SVG output of draft maps in the style of the BRENDA map."""

from __future__ import annotations

import html
import json

from .core import DraftMap
from .routes import bezier_path

STYLE = """
polygon.pw, rect.pw{stroke:none;}
text{font-family:Liberation Sans, Arial, Helvetica, sans-serif;text-anchor:middle;}
path.link{stroke-width:2;stroke:black;fill:none;}
path.inter-pw-link{stroke:rgb(90,90,90);stroke-width:4;stroke-dasharray:14 8;opacity:0.55;}
path.inter-pw-link.added{stroke:rgb(150,60,170);}
rect{stroke-linejoin:bevel;}
g.nodeEN rect{fill:rgb(139, 179, 255);stroke:rgb(97, 125, 179);}
g.nodeEN text{font-size:12px;}
g.nodeSM rect, g.nodeHM rect{fill:rgb(255, 183, 139);stroke:rgb(179, 128, 97);}
g.nodeSM text, g.nodeHM text{font-size:12px;}
g.nodePathwayTitle rect{fill:rgb(255, 255, 255);stroke:rgb(179, 179, 179);}
g.nodePathwayTitle text{font-size:18px;font-weight:bolder;}
g.nodeCM rect{fill:rgb(212, 207, 114);stroke:rgb(148, 145, 80);}
g.nodeCM text{font-size:12px;}
g.nodeEN.ecm-on-map rect{fill:rgb(110, 200, 100);stroke:rgb(50, 130, 45);stroke-width:3;}
g.nodeEN.ecm-vicinity rect{fill:rgb(80, 200, 200);stroke:rgb(30, 120, 130);stroke-width:3;}
g.nodeEN.ecm-added rect{fill:rgb(205, 160, 230);stroke:rgb(130, 60, 160);stroke-width:3;stroke-dasharray:6 3;}
g.ecm-added-met rect{stroke-dasharray:5 3;stroke-width:2;}
g.hub rect{stroke:rgb(190, 40, 40);stroke-width:3;}
text.rxn{font-size:9px;fill:rgb(40,40,40);}
text.ref{font-size:11px;font-style:italic;fill:rgb(130,50,160);}
a{fill:inherit;}
"""

MARKER = (
    '<defs><marker style="fill: black;" orient="auto-start-reverse" markerUnits="userSpaceOnUse" '
    'markerHeight="10" markerWidth="10" refY="0" refX="8" viewBox="0 -5 10 10" id="in">'
    '<path d="M0,-5L10,0L0,5"/></marker></defs>'
)

BRENDA = "https://www.brenda-enzymes.org"


def _clip(cx, cy, w, h, tx, ty, margin=1.0):
    """Point where the line from the centre of a box to (tx, ty) leaves the box."""
    dx, dy = tx - cx, ty - cy
    if dx == 0 and dy == 0:
        return cx, cy
    t = min((w / 2 + margin) / abs(dx) if dx else float("inf"), (h / 2 + margin) / abs(dy) if dy else float("inf"))
    t = min(t, 1.0)
    return cx + t * dx, cy + t * dy


def _line(a, b):
    p = _clip(a.x, a.y, a.w, a.h, b.x, b.y)
    q = _clip(b.x, b.y, b.w, b.h, a.x, a.y)
    return f"M {p[0]:.1f} {p[1]:.1f} L {q[0]:.1f} {q[1]:.1f}"


def write_svg(dm: DraftMap, path: str, reaction_labels: bool = True, md=None, margin: float = 60.0,
              show_links: bool = True) -> None:
    nodes = dm.nodes
    xs = [s.x for s in dm.subsystems.values()] + [n.x - n.w / 2 for n in nodes.values()]
    ys = [s.y - s.title_h / 2 for s in dm.subsystems.values()] + [n.y - n.h / 2 for n in nodes.values()]
    xe = [s.x + s.w for s in dm.subsystems.values()] + [n.x + n.w / 2 for n in nodes.values()]
    ye = [s.y + s.h for s in dm.subsystems.values()] + [n.y + n.h / 2 for n in nodes.values()]
    x0, y0 = min(xs) - margin, min(ys) - margin
    W, H = max(xe) - x0 + margin, max(ye) - y0 + margin
    dm.view_origin = (x0, y0)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="{x0:.1f} {y0:.1f} {W:.1f} {H:.1f}" '
        f'width="{W:.0f}" height="{H:.0f}" style="font-family: Arial, Helvetica, sans-serif; font-size: 15px;">',
        f"<style>{STYLE}</style>",
        MARKER,
        f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{W:.1f}" height="{H:.1f}" style="fill:white"/>',
        '<g class="pw-con" id="pw-con">',
    ]
    for s in dm.subsystems.values():
        out.append(
            f'<g class="pw-sub"><rect id="pw_{html.escape(s.id)}" class="pw {s.color_group}" x="{s.x:.1f}" '
            f'y="{s.y:.1f}" width="{s.w:.1f}" height="{s.h:.1f}" style="opacity: 0.5; fill: {s.fill};"/></g>'
        )
    out.append("</g>")
    out.append('<g class="all-link" id="all-link">')
    if show_links:
        for l in dm.links:
            a, b = nodes[l.a], nodes[l.b]
            css = "link inter-pw-link" + (" added" if l.origin == "added" else "")
            out.append(f'<path class="{css}" d="{_line(a, b)}"><title>{html.escape(a.label)}</title></path>')
    for e in dm.edges:
        if e.kind != "main":
            continue
        a, b = nodes[e.enzyme], nodes[e.met]
        markers = (' marker-start="url(#in)"' if e.arrow_enzyme else "") + (' marker-end="url(#in)"' if e.arrow_met else "")
        route = getattr(e, "route", None)  # around other nodes, see .routes
        out.append(f'<path class="link"{markers} d="{bezier_path(route) if route else _line(a, b)}"/>')
    out.append("</g>")
    out.append('<g class="all-linkCof" id="all-linkCof">')
    for e in dm.edges:
        if e.kind == "main":
            continue
        a, b = nodes[e.enzyme], nodes[e.met]
        markers = (' marker-start="url(#in)"' if e.arrow_enzyme else "") + (' marker-end="url(#in)"' if e.arrow_met else "")
        out.append(f'<path class="link"{markers} d="{_line(a, b)}"/>')
    out.append("</g>")
    out.append('<g class="all-node" id="all-node">')
    for s in dm.subsystems.values():
        tx, ty = s.x + s.w / 2, s.y
        out.append(
            f'<g class="node nodePathwayTitle pw_{html.escape(s.id)}" transform="translate({tx:.1f},{ty:.1f})">'
            f'<rect rx="10" height="{s.title_h:.0f}" width="{s.title_w:.0f}" x="{-s.title_w / 2:.1f}" '
            f'y="{-s.title_h / 2:.1f}"/><text dominant-baseline="central">{html.escape(s.label)}</text></g>'
        )
    for kind in ("main", "enzyme", "cofactor"):
        if kind == "cofactor":
            out.append('</g><g class="all-nodeCof" id="all-nodeCof">')
        for n in nodes.values():
            if n.kind != kind:
                continue
            css = f"node {n.svg_class} pw_{html.escape(n.subsystem)}"
            if n.ligand_id is not None:
                css += f" lid{n.ligand_id}"
            if n.kind == "enzyme":
                css += {"on-map": " ecm-on-map", "vicinity": " ecm-vicinity", "added": " ecm-added"}.get(n.status, "")
            elif n.status in ("added", "linked", "new"):
                css += " ecm-added-met"
            if getattr(n, "hub", False):
                css += " hub"
            rx = {"enzyme": ' rx="15"', "cofactor": ' rx="5"'}.get(n.kind, "")
            if n.kind == "enzyme" and n.status != "added" and n.ec_label and "?" not in n.ec_label:
                href = f"{BRENDA}/enzyme.php?ecno={n.ec_label.split('/')[0]}"
            elif n.kind != "enzyme" and n.ligand_id is not None:
                href = f"{BRENDA}/ligand.php?brenda_ligand_id={n.ligand_id}"
            else:
                href = None
            label = html.escape(n.label)
            text = f'<a href="{href}" target="_blank">{label}</a>' if href else label
            title = n.title or ""
            if n.kind == "enzyme" and n.reactions:
                if md is not None:
                    parts = [f"{r}: {md.reactions[r]['name']}" for r in n.reactions]
                else:
                    parts = n.reactions
                title = (title + " | " if title else "") + ("added from model: " if n.status == "added" else "model: ") + "; ".join(parts)
                if n.status == "added" and n.ec_label:
                    title += f" | EC {n.ec_label}"
            if getattr(n, "hub", False):
                title += " | hub metabolite: drawn in several subsystems, links between them not drawn"
            extra = ""
            stubs = getattr(n, "stubs", None)
            if stubs:
                # with --split: the metabolite continues on other maps
                extra = f'<text class="ref" y="{n.h / 2 + 11:.1f}">→ {html.escape(", ".join(stubs))}</text>'
                title += " | also on: " + ", ".join(stubs)
            if reaction_labels and n.kind == "enzyme" and n.reactions and n.status != "added":
                ids = ", ".join(n.reactions[:3]) + (" ..." if len(n.reactions) > 3 else "")
                extra = f'<text class="rxn" y="{n.h / 2 + 9:.1f}">{html.escape(ids)}</text>'
            out.append(
                f'<g id="{n.id}" class="{css}" transform="translate({n.x:.1f},{n.y:.1f})">'
                f'<rect{rx} height="{n.h:.1f}" width="{n.w:.1f}" x="{-n.w / 2:.1f}" y="{-n.h / 2:.1f}"/>'
                f'<text dominant-baseline="central">{text}</text>{extra}<title>{html.escape(title)}</title></g>'
            )
    out.append("</g></svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def write_reaction_positions(dm: DraftMap, path: str) -> int:
    """Positions of the model reactions on the map, {reaction ID: [[x, y], ...]} (enzyme node centres)."""
    positions: dict[str, list] = {}
    for n in dm.nodes.values():
        for r in n.reactions:
            positions.setdefault(r, []).append([round(n.x, 1), round(n.y, 1)])
    with open(path, "w", encoding="utf-8") as f:
        json.dump(positions, f, indent=0, sort_keys=True)
    return len(positions)


def write_cnapy_maps(dm: DraftMap, path: str, box_width: float = 80.0, box_height: float = 24.0) -> int:
    """CNApy box positions ({reaction ID: [x, y]}, a ``.maps`` file) for the SVG written by :func:`write_svg`.

    With the SVG as background (bg-size 1), one map unit is one scene unit
    with the origin at the top left of the SVG. Each reaction box is put on
    the first enzyme node of the reaction (nodes of the map before added ones).
    """
    x0, y0 = getattr(dm, "view_origin", (0.0, 0.0))
    order = {"on-map": 0, "vicinity": 1, "added": 2}
    boxes = {}
    for n in sorted((n for n in dm.nodes.values() if n.kind == "enzyme"), key=lambda n: order.get(n.status, 3)):
        for r in n.reactions:
            boxes.setdefault(r, [round(n.x - x0 - box_width / 2, 1), round(n.y - y0 - box_height / 2, 1)])
    with open(path, "w", encoding="utf-8") as f:
        json.dump(boxes, f)
    return len(boxes)
