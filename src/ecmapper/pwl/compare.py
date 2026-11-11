"""Compare layouts of PWL subsystems with the original drawings of the BRENDA map.

The original positions come from ``BRENDA_map_reactions.json``
(:meth:`ecmapper.draft_map.core.DraftMap.from_brenda`); nodes are matched by ID:
enzyme nodes by reaction ID, main metabolites by ``node_index.json``. Cofactor
nodes are not compared (PWL does not identify them).

Metrics (:func:`similarity`), all invariant to translation, scale, rotation
and reflection of the generated layout:

``procrustes``
    Procrustes disparity of the enzyme and main-metabolite positions (0: same
    shape, 1: unrelated), after the best scaling, rotation and reflection.
``directions``
    Share of the main links (enzyme to metabolite) that point in the same of
    the four directions as in the original, for the best of the eight
    symmetries of the square.
``neighbours``
    Mean share of each node's three nearest nodes that are also among its
    three nearest in the original.
``orthogonal``
    Share of main links within 10 degrees of horizontal or vertical
    (``orthogonal_original`` for the original).
``aspect``
    |log| of the ratio of the aspect ratios (width/height) of the two
    drawings, allowing a rotation by 90 degrees.
"""

from __future__ import annotations

import csv
import html
import json
import math
import re
import statistics
from pathlib import Path

import numpy as np

from ..draft_map.core import DraftMap

SYMMETRIES = [((1, 0), (0, 1)), ((0, -1), (1, 0)), ((-1, 0), (0, -1)), ((0, 1), (-1, 0)),
              ((-1, 0), (0, 1)), ((1, 0), (0, -1)), ((0, 1), (1, 0)), ((0, -1), (-1, 0))]


def original_subsystem(original: DraftMap, sid: str) -> DraftMap:
    """The original drawing of one subsystem."""
    from ..draft_map.relayout import subsystem_map

    return subsystem_map(original, sid)


def procrustes(a: np.ndarray, b: np.ndarray) -> float:
    """Procrustes disparity of two point sets (rows matched), with scaling, rotation and reflection."""
    if len(a) < 3:
        return 0.0
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-9 or nb < 1e-9:
        return 1.0
    a, b = a / na, b / nb
    s = np.linalg.svd(a.T @ b, compute_uv=False)
    return float(max(0.0, 1.0 - s.sum() ** 2))


def _dir4(v) -> tuple[int, int]:
    if abs(v[0]) >= abs(v[1]):
        return (1 if v[0] >= 0 else -1, 0)
    return (0, 1 if v[1] >= 0 else -1)


def _orthogonal(vectors, tolerance=10.0) -> float:
    if not vectors:
        return 1.0
    t = math.radians(tolerance)
    ok = 0
    for x, y in vectors:
        a = math.atan2(abs(y), abs(x))
        ok += a <= t or a >= math.pi / 2 - t
    return ok / len(vectors)


def _aspect(points) -> float:
    p = np.array(points)
    w, h = np.ptp(p[:, 0]) + 1.0, np.ptp(p[:, 1]) + 1.0
    return w / h


def similarity(generated: DraftMap, original: DraftMap) -> dict:
    """The metrics of the module doc for two drawings of a subsystem (nodes matched by ID)."""
    common = [i for i, n in generated.nodes.items() if n.kind != "cofactor" and i in original.nodes]
    out = dict(matched=len(common))
    if len(common) < 2:
        return out
    a = np.array([(generated.nodes[i].x, generated.nodes[i].y) for i in common])
    b = np.array([(original.nodes[i].x, original.nodes[i].y) for i in common])
    out["procrustes"] = round(procrustes(b, a), 4)
    links = [(e.enzyme, e.met) for e in generated.edges if e.kind == "main"
             and e.enzyme in original.nodes and e.met in original.nodes]
    vg = [(generated.nodes[m].x - generated.nodes[e].x, generated.nodes[m].y - generated.nodes[e].y) for e, m in links]
    vo = [(original.nodes[m].x - original.nodes[e].x, original.nodes[m].y - original.nodes[e].y) for e, m in links]
    if links:
        best = 0
        for (r0, r1) in SYMMETRIES:
            same = sum(_dir4((r0[0] * x + r0[1] * y, r1[0] * x + r1[1] * y)) == _dir4(o) for (x, y), o in zip(vg, vo))
            best = max(best, same)
        out["directions"] = round(best / len(links), 4)
    k = min(3, len(common) - 1)
    da = np.linalg.norm(a[:, None] - a[None], axis=2)
    db = np.linalg.norm(b[:, None] - b[None], axis=2)
    np.fill_diagonal(da, np.inf)
    np.fill_diagonal(db, np.inf)
    na = np.argsort(da, axis=1)[:, :k]
    nb = np.argsort(db, axis=1)[:, :k]
    out["neighbours"] = round(float(np.mean([len(set(x) & set(y)) / k for x, y in zip(na, nb)])), 4)
    out["orthogonal"] = round(_orthogonal(vg), 4)
    out["orthogonal_original"] = round(_orthogonal(vo), 4)
    ra, rb = _aspect(a), _aspect(b)
    out["aspect"] = round(min(abs(math.log(ra / rb)), abs(math.log(ra * rb))), 4)
    return out


# ---------------------------------------------------------------------------------- drawings


def side_by_side(left: DraftMap, right: DraftMap, captions=("", ""), gap: float = 80.0) -> str:
    """One SVG with two drawings next to each other (drawn with :func:`ecmapper.draft_map.render.svg_text`)."""
    from ..draft_map.render import svg_text

    parts, x, height = [], 0.0, 0.0
    for dm, caption in zip((left, right), captions):
        text = svg_text(dm, reaction_labels=False, margin=40)
        vx, vy, w, h = dm.view_box
        inner = re.sub(r"^<svg ", f'<svg x="{x:.1f}" y="40" ', text, count=1)
        parts.append(inner)
        parts.append(f'<text x="{x + 10:.1f}" y="28" style="font: bold 22px Arial, sans-serif">'
                     f'{html.escape(caption)}</text>')
        x += w + gap
        height = max(height, h + 40)
    width = x - gap
    return (f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="{width:.0f}" height="{height:.0f}" '
            f'viewBox="0 0 {width:.1f} {height:.1f}"><rect width="100%" height="100%" fill="white"/>'
            + "\n".join(parts) + "</svg>\n")


# ---------------------------------------------------------------------------------- report

COLUMNS = ["subsystem", "method", "chosen", "nodes", "procrustes", "directions", "neighbours", "orthogonal",
           "aspect", "crossings_main", "crossings_all", "overlapping_nodes", "orthogonal_original",
           "crossings_main_original", "crossings_all_original"]


def write_report(rows: list[dict], outdir: Path, failures: list[str] | None = None) -> str:
    """``metrics.csv``, ``metrics.json`` and ``summary.md`` (returned) from the rows of :func:`row`."""
    outdir = Path(outdir)
    with open(outdir / "metrics.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(outdir / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=1)
    methods = list(dict.fromkeys(r["method"] for r in rows))
    by = {m: [r for r in rows if r["method"] == m] for m in methods}
    chosen = [r for r in rows if r["chosen"]]
    by["chosen"] = chosen

    def med(m, key):
        values = [r[key] for r in by[m] if r.get(key) is not None]
        return statistics.median(values) if values else float("nan")

    def mean(m, key):
        values = [r[key] for r in by[m] if r.get(key) is not None]
        return statistics.fmean(values) if values else float("nan")

    subsystems = sorted({r["subsystem"] for r in rows})
    lines = [
        "# Layouts from PWL: comparison with the BRENDA map", "",
        f"{len(subsystems)} subsystems. Median over subsystems (crossings: total). "
        "procrustes: lower is more similar; directions, neighbours, orthogonal: higher is more similar; "
        "aspect: |log| ratio of aspect ratios. See `ecmapper/pwl/compare.py`.", "",
        "| method | subsystems | procrustes | directions | neighbours | orthogonal | aspect | main crossings | "
        "all crossings | overlapping nodes | best procrustes | chosen |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    best_count = {m: 0 for m in methods}
    for sid in subsystems:
        rs = [r for r in rows if r["subsystem"] == sid and r.get("procrustes") is not None]
        if rs:
            best_count[min(rs, key=lambda r: r["procrustes"])["method"]] += 1
    chosen_count = {m: sum(1 for r in chosen if r["method"] == m) for m in methods}
    for m in methods + ["chosen"]:
        total = lambda key: sum(r.get(key) or 0 for r in by[m])  # noqa: E731
        lines.append(
            f"| {m} | {len(by[m])} | {med(m, 'procrustes'):.3f} | {med(m, 'directions'):.3f} | "
            f"{med(m, 'neighbours'):.3f} | {med(m, 'orthogonal'):.3f} | {med(m, 'aspect'):.2f} | "
            f"{total('crossings_main')} | {total('crossings_all')} | {total('overlapping_nodes')} | "
            f"{best_count.get(m, '') if m != 'chosen' else ''} | {chosen_count.get(m, '') if m != 'chosen' else ''} |")
    original_main = sum(r.get("crossings_main_original") or 0 for r in chosen)
    original_all = sum(r.get("crossings_all_original") or 0 for r in chosen)
    lines += ["", f"Original drawings: {original_main} main crossings, {original_all} crossings in all; "
              f"median share of orthogonal main links {med('chosen', 'orthogonal_original'):.3f}.",
              f"Mean procrustes of the chosen layouts: {mean('chosen', 'procrustes'):.3f}."]
    if failures:
        lines += ["", "Failed: " + "; ".join(failures)]
    text = "\n".join(lines) + "\n"
    (outdir / "summary.md").write_text(text, encoding="utf-8")
    return text


def row(sid: str, method: str, chosen: bool, dm: DraftMap, original: DraftMap | None, metrics: dict,
        original_metrics: dict | None) -> dict:
    r = dict(subsystem=sid, method=method, chosen=chosen,
             nodes=sum(1 for n in dm.nodes.values() if n.kind != "cofactor"), **metrics)
    r.pop("drawing_area", None)
    if original is not None:
        r.update(similarity(dm, original))
    if original_metrics:
        r["crossings_main_original"] = original_metrics["crossings_main"]
        r["crossings_all_original"] = original_metrics["crossings_all"]
    return r


__all__ = ["original_subsystem", "procrustes", "row", "side_by_side", "similarity", "write_report"]


# ---------------------------------------------------------------------------------- whole directories


def layout_directory(pwl_dir, outdir, brenda_json=None, brenda_svg=None, methods=None, subsystems=None,
                     routes: bool = True, log=print) -> list[dict]:
    """Lay out the PWL files of `pwl_dir` (``subsystems/*.pwl``, ``node_index.json``) and write to `outdir`:

    - ``layout/<id>.svg``: the chosen layout, in the style of the BRENDA map;
    - with `brenda_json`: ``compare/<id>.svg`` (original and chosen layout side by side), ``metrics.csv``,
      ``metrics.json`` and ``summary.md`` (all candidate layouts, see :func:`similarity`).

    Returns the rows of the metrics.
    """
    from ..draft_map.layout import internal_metrics
    from ..draft_map.render import write_svg
    from .layout import METHODS, choose, layout_candidates, metrics, route
    from .syntax import parse_file

    pwl_dir, outdir = Path(pwl_dir), Path(outdir)
    (outdir / "layout").mkdir(parents=True, exist_ok=True)
    index_path = pwl_dir / "node_index.json"
    node_index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else {}
    original = DraftMap.from_brenda(str(brenda_json), str(brenda_svg) if brenda_svg else None) if brenda_json else None
    if original is not None:
        (outdir / "compare").mkdir(exist_ok=True)
    files = sorted((pwl_dir / "subsystems").glob("*.pwl")) if (pwl_dir / "subsystems").is_dir() else sorted(
        pwl_dir.glob("*.pwl"))
    rows, failures = [], []
    for path in files:
        sid = path.stem
        if subsystems and sid not in subsystems:
            continue
        try:
            doc = parse_file(path)
            orig = original_subsystem(original, doc.id) if original is not None and doc.id in original.subsystems \
                else None
            fill = orig.subsystems[doc.id].fill if orig is not None else "#dedede"
            candidates = layout_candidates(doc, node_index.get(doc.id, {}), methods or METHODS, fill)
            chosen = choose(candidates)
            dm = candidates[chosen]
            if routes:
                route(dm)
            write_svg(dm, str(outdir / "layout" / f"{sid}.svg"), reaction_labels=False)
            if orig is not None:
                om = internal_metrics(orig)
                for m, c in candidates.items():
                    rows.append(row(doc.id, m, m == chosen, c, orig, metrics(c), om))
                (outdir / "compare" / f"{sid}.svg").write_text(
                    side_by_side(orig, dm, ("BRENDA map", f"from PWL ({chosen})")), encoding="utf-8")
            log(f"{sid}: {chosen}")
        except Exception as e:  # noqa: BLE001 - report and go on with the other subsystems
            failures.append(f"{sid}: {e}")
            log(f"{sid}: FAILED {e}")
    if rows:
        write_report(rows, outdir, failures)
    return rows
