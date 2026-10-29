"""Condensing subsystems and arranging them with few crossings of inter-pathway links."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

from .core import DraftMap

# --- 1D compaction (order preserving, LP) -----------------------------------------------------------


def compact_axis(pos, size, opos, osize, gap, ogap, compactness=0.5, max_scale=1.0, keep_order=False):
    """Positions along one axis that keep the order of overlapping boxes and are compact.

    Boxes whose extents in the other axis overlap (with `ogap`) keep their
    order along this axis and are at least `gap` apart. Among such positions,
    the LP minimises the L1 distance to the original positions scaled by a
    common factor s <= `max_scale` plus `compactness` * n * extent: the layout
    is shrunk uniformly where possible and nodes are moved individually only
    to close larger gaps. Returns positions with the minimum lower edge at 0.
    """
    pos, size = np.asarray(pos, float), np.asarray(size, float)
    opos, osize = np.asarray(opos, float), np.asarray(osize, float)
    n = len(pos)
    if n == 0:
        return pos
    origin = pos - (pos - size / 2).min()
    if n == 1:
        return size / 2
    lo, hi = opos - osize / 2 - ogap / 2, opos + osize / 2 + ogap / 2
    overlap = (lo[:, None] < hi[None, :]) & (lo[None, :] < hi[:, None])
    idx = np.arange(n)
    before = (pos[:, None] < pos[None, :]) | ((pos[:, None] == pos[None, :]) & (idx[:, None] < idx[None, :]))
    pairs = np.argwhere(overlap & before)
    if keep_order is True:
        # all nodes keep their order along the axis (consecutive pairs suffice), overlapping ones with gap
        order = np.lexsort((idx, pos))
        chain = np.stack([order[:-1], order[1:]], axis=1)
        pairs = np.concatenate([pairs, chain[~overlap[chain[:, 0], chain[:, 1]]]])
    elif keep_order:
        # nodes closer than `keep_order` in the other axis keep their order
        near = (lo[:, None] - keep_order < hi[None, :]) & (lo[None, :] - keep_order < hi[:, None]) & ~overlap
        pairs = np.concatenate([pairs, np.argwhere(near & before)])
    # variables: p (n), t (n), s, E
    nv = 2 * n + 2
    s_col, e_col = 2 * n, 2 * n + 1
    rows, cols, vals, b = [], [], [], []
    r = 0
    for i, j in pairs:  # p_i - p_j <= -sep
        rows += [r, r]
        cols += [i, j]
        vals += [1.0, -1.0]
        b.append(-((size[i] + size[j]) / 2 + gap) if overlap[i, j] else 0.0)
        r += 1
    for i in range(n):
        # p_i - s o_i - t_i <= 0 ; -p_i + s o_i - t_i <= 0 ; p_i - E <= -size_i/2
        rows += [r, r, r]
        cols += [i, s_col, n + i]
        vals += [1.0, -origin[i], -1.0]
        b.append(0.0)
        r += 1
        rows += [r, r, r]
        cols += [i, s_col, n + i]
        vals += [-1.0, origin[i], -1.0]
        b.append(0.0)
        r += 1
        rows += [r, r]
        cols += [i, e_col]
        vals += [1.0, -1.0]
        b.append(-size[i] / 2)
        r += 1
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, nv))
    c = np.zeros(nv)
    c[n : 2 * n] = 1.0
    c[e_col] = compactness * n
    bounds = [(size[i] / 2, None) for i in range(n)] + [(0, None)] * n + [(0, max_scale), (0, None)]
    res = linprog(c, A_ub=A, b_ub=np.array(b), bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"compaction failed: {res.message}")
    return res.x[:n]


def compact_boxes(x, y, w, h, gap_x, gap_y, compactness=0.5, rounds=2, max_scale=1.0, keep_order=False):
    """Alternating x and y compaction of boxes given by centres and sizes (see :func:`compact_axis`)."""
    x, y = np.asarray(x, float).copy(), np.asarray(y, float).copy()
    for _ in range(rounds):
        x = compact_axis(x, w, y, h, gap_x, gap_y, compactness, max_scale, keep_order)
        y = compact_axis(y, h, x, w, gap_y, gap_x, compactness, max_scale, keep_order)
    return x, y


TITLE_GAP = 22.0
PADDING = 35.0


def condense_subsystems(dm: DraftMap, gap_x=22.0, gap_y=28.0, compactness=0.5, rounds=2, keep_order=150.0) -> None:
    """Step 4: condense the drawing of each subsystem and fit its box (and title) around it."""
    for s in dm.subsystems.values():
        members = dm.members(s.id)
        if not members:
            continue
        x = np.array([n.x for n in members])
        y = np.array([n.y for n in members])
        w = np.array([n.w for n in members])
        h = np.array([n.h for n in members])
        # room for the model reaction IDs below enzyme nodes of the map
        label = np.array([n.kind == "enzyme" and bool(n.reactions) and n.status != "added" for n in members])
        lw = np.array([0.54 * 9 * len(", ".join(n.reactions[:3])) + 4 for n in members])
        we = np.where(label, np.maximum(w, lw), w)
        he = np.where(label, h + 24, h)
        ye = np.where(label, y + 12, y)
        x, ye = compact_boxes(x, ye, we, he, gap_x, gap_y, compactness, rounds, keep_order=keep_order)
        y = np.where(label, ye - 12, ye)
        h = he
        w = we
        x0, x1 = (x - w / 2).min(), (x + w / 2).max()
        y0, y1 = (ye - h / 2).min(), (ye + h / 2).max()
        content_w = x1 - x0
        box_w = max(content_w + 2 * PADDING, s.title_w + 2 * PADDING)
        top = s.title_h / 2 + TITLE_GAP
        box_h = (y1 - y0) + top + PADDING
        # node positions relative to the box (top left at 0, 0), content centred horizontally
        dx = (box_w - content_w) / 2 - x0
        dy = top - y0
        for n, nx, ny in zip(members, x, y):
            n.x, n.y = nx + dx, ny + dy
        s.x, s.y, s.w, s.h = 0.0, 0.0, box_w, box_h


def to_relative(dm: DraftMap) -> None:
    """Express node positions relative to the top left corner of their subsystem box."""
    for n in dm.nodes.values():
        s = dm.subsystems[n.subsystem]
        n.x -= s.x
        n.y -= s.y


def to_absolute(dm: DraftMap) -> None:
    for n in dm.nodes.values():
        s = dm.subsystems[n.subsystem]
        n.x += s.x
        n.y += s.y


def fit_boxes_original(dm: DraftMap) -> None:
    """Boxes of the original map around the (uncondensed) drawing of each subsystem, positions relative."""
    for s in dm.subsystems.values():
        members = dm.members(s.id)
        x0 = min(n.x - n.w / 2 for n in members)
        x1 = max(n.x + n.w / 2 for n in members)
        y0 = min(n.y - n.h / 2 for n in members)
        y1 = max(n.y + n.h / 2 for n in members)
        # keep the original box if the content fits into it
        bx0, by0 = min(s.x, x0 - PADDING), min(s.y, y0 - s.title_h / 2 - TITLE_GAP)
        bx1, by1 = max(s.x + s.w, x1 + PADDING), max(s.y + s.h, y1 + PADDING)
        s.x, s.y, s.w, s.h = bx0, by0, bx1 - bx0, by1 - by0


# --- inter-pathway links and metrics -------------------------------------------------------------


@dataclass
class Ports:
    """The inter-pathway links as pairs of (subsystem index, position relative to the box)."""

    names: list[str]  # subsystem IDs, index = box
    w: np.ndarray
    h: np.ndarray
    a: np.ndarray  # box index of the first end of each link
    b: np.ndarray
    pa: np.ndarray  # (L, 2) relative positions
    pb: np.ndarray
    ends: list[tuple[str, str]]
    my: np.ndarray  # y' = my - y mirrors the drawing of a box vertically below its title

    @classmethod
    def from_map(cls, dm: DraftMap) -> Ports:
        names = list(dm.subsystems)
        index = {k: i for i, k in enumerate(names)}
        a, b, pa, pb, ends = [], [], [], [], []
        for l in dm.links:
            na, nb = dm.nodes[l.a], dm.nodes[l.b]
            if na.subsystem == nb.subsystem:
                continue
            a.append(index[na.subsystem])
            b.append(index[nb.subsystem])
            pa.append((na.x, na.y))
            pb.append((nb.x, nb.y))
            ends.append((l.a, l.b))
        w = np.array([dm.subsystems[k].w for k in names])
        h = np.array([dm.subsystems[k].h for k in names])
        my = np.array([dm.subsystems[k].title_h / 2 + TITLE_GAP + dm.subsystems[k].h - PADDING for k in names])
        return cls(names, w, h, np.array(a, int), np.array(b, int), np.array(pa, float).reshape(-1, 2),
                   np.array(pb, float).reshape(-1, 2), ends, my)

    def endpoints(self, x, y, fx=None, fy=None):
        """Absolute link end points for box positions (top left) and flips."""
        pa, pb = self.pa.copy(), self.pb.copy()
        if fx is not None:
            pa[:, 0] = np.where(fx[self.a], self.w[self.a] - pa[:, 0], pa[:, 0])
            pb[:, 0] = np.where(fx[self.b], self.w[self.b] - pb[:, 0], pb[:, 0])
        if fy is not None:
            pa[:, 1] = np.where(fy[self.a], self.my[self.a] - pa[:, 1], pa[:, 1])
            pb[:, 1] = np.where(fy[self.b], self.my[self.b] - pb[:, 1], pb[:, 1])
        p1 = pa + np.stack([x[self.a], y[self.a]], axis=1)
        p2 = pb + np.stack([x[self.b], y[self.b]], axis=1)
        return p1, p2


def _orient(o, a, b):
    return (a[..., 0] - o[..., 0]) * (b[..., 1] - o[..., 1]) - (a[..., 1] - o[..., 1]) * (b[..., 0] - o[..., 0])


def count_crossings(p1, p2) -> int:
    """Number of proper crossings between line segments (segments sharing an end point do not cross)."""
    if len(p1) < 2:
        return 0
    A1, A2 = p1[:, None, :], p2[:, None, :]
    B1, B2 = p1[None, :, :], p2[None, :, :]
    o1, o2 = _orient(A1, A2, B1), _orient(A1, A2, B2)
    o3, o4 = _orient(B1, B2, A1), _orient(B1, B2, A2)
    x = (o1 * o2 < 0) & (o3 * o4 < 0)
    return int(np.triu(x, 1).sum())


def count_box_crossings(p1, p2, a, b, x, y, w, h, shrink=10.0) -> int:
    """Number of (link, box) pairs where a link passes through a box other than those of its ends."""
    if len(p1) == 0:
        return 0
    x0 = (x + shrink)[None, :]
    x1 = (x + w - shrink)[None, :]
    y0 = (y + shrink)[None, :]
    y1 = (y + h - shrink)[None, :]
    px, py = p1[:, 0:1], p1[:, 1:2]
    dx, dy = (p2[:, 0] - p1[:, 0])[:, None], (p2[:, 1] - p1[:, 1])[:, None]
    t0 = np.zeros((len(p1), len(x)))
    t1 = np.ones((len(p1), len(x)))
    ok = np.ones_like(t0, bool)
    for p, q in ((-dx, px - x0), (dx, x1 - px), (-dy, py - y0), (dy, y1 - py)):
        p = np.broadcast_to(p, t0.shape)
        q = np.broadcast_to(q, t0.shape)
        zero = np.abs(p) < 1e-12
        ok &= ~(zero & (q < 0))
        with np.errstate(divide="ignore", invalid="ignore"):
            r = q / p
        t0 = np.where(~zero & (p < 0), np.maximum(t0, r), t0)
        t1 = np.where(~zero & (p > 0), np.minimum(t1, r), t1)
    hit = ok & (t0 < t1)
    own = (np.arange(len(x))[None, :] == a[:, None]) | (np.arange(len(x))[None, :] == b[:, None])
    return int((hit & ~own).sum())


def evaluate_arrangement(ports: Ports, x, y, fx=None, fy=None) -> dict:
    p1, p2 = ports.endpoints(x, y, fx, fy)
    W = float((x + ports.w).max() - x.min()) if len(x) else 0.0
    H = float((y + ports.h).max() - y.min()) if len(y) else 0.0
    return dict(
        links=len(p1),
        crossings=count_crossings(p1, p2),
        box_crossings=count_box_crossings(p1, p2, ports.a, ports.b, x, y, ports.w, ports.h),
        link_length=float(np.linalg.norm(p2 - p1, axis=1).sum()),
        width=W,
        height=H,
        area_ratio=W * H / float((ports.w * ports.h).sum()),
    )


# --- arrangements --------------------------------------------------------------------------------------


def arrangement_original(dm: DraftMap, gap=80.0, compactness=0.5):
    """Keep the arrangement of the original map, close the gaps (order preserving compaction of the boxes).

    The boxes are placed at their original centres, scaled with the size of
    the (condensed) boxes, then compacted. Needs `orig_centres` on the map.
    """
    names = list(dm.subsystems)
    w = np.array([dm.subsystems[k].w for k in names])
    h = np.array([dm.subsystems[k].h for k in names])
    cx = np.array([dm.orig_centres[k][0] for k in names])
    cy = np.array([dm.orig_centres[k][1] for k in names])
    x, y = compact_boxes(cx, cy, w + gap, h + gap, 0.0, 0.0, compactness, rounds=3, max_scale=10.0)
    return x - w / 2, y - h / 2


class SequencePair:
    """Packing of boxes encoded by a sequence pair (Murata et al. 1995).

    a before b in both sequences: a is left of b; a before b in the first and
    after b in the second sequence: a is above b. The packing places each box
    as far left and up as the relations allow.
    """

    def __init__(self, w, h, gap):
        self.w, self.h = np.asarray(w, float) + gap, np.asarray(h, float) + gap
        self.gap = gap
        self.n = len(w)

    @staticmethod
    def from_positions(cx, cy):
        cx, cy = np.asarray(cx, float), np.asarray(cy, float)
        return list(np.argsort(cx + cy, kind="stable")), list(np.argsort(cx - cy, kind="stable"))

    def pack(self, gp, gn):
        n = self.n
        gp = np.asarray(gp)
        rank = np.empty(n, int)
        rank[np.asarray(gn)] = np.arange(n)
        x, y = np.zeros(n), np.zeros(n)
        for k in range(1, n):
            b = gp[k]
            prev = gp[:k]
            left = rank[prev] < rank[b]
            if left.any():
                pl = prev[left]
                x[b] = (x[pl] + self.w[pl]).max()
            if (~left).any():
                pa = prev[~left]
                y[b] = (y[pa] + self.h[pa]).max()
        return x, y


@dataclass
class AnnealWeights:
    crossing: float = 1.0
    box_crossing: float = 0.5
    length: float = 0.02  # per link length in units of the mean box size
    area: float = 100.0  # per (bounding box area / sum of box areas)
    aspect: float = 200.0  # per log(aspect ratio / target)^2
    target_aspect: float = 1.6


def anneal_arrangement(ports: Ports, init_x, init_y, gap=80.0, iterations=200000, seed=1, flips=False,
                       weights: AnnealWeights | None = None, t_start=None, t_end=0.05, log=None):
    """Arrange boxes with few link crossings: simulated annealing over sequence pairs (and flips).

    `flips`: False, ``"x"`` (mirror subsystems left-right) or ``"xy"`` (also upside down).

    The cost is a weighted sum of the crossings of links, the links passing
    through other boxes, the total link length, the area of the bounding box
    and the deviation from the target aspect ratio. Moves swap two boxes in
    one or both sequences or mirror a box. Starts from the arrangement given
    by the box positions (top left). Returns x, y (top left), flip_x, flip_y
    and the evaluation.
    """
    from .fast import metrics, pack

    weights = weights or AnnealWeights()
    rng = np.random.default_rng(seed)
    n = len(ports.w)
    W_, H_ = ports.w + gap, ports.h + gap
    a, b = ports.a.astype(np.int64), ports.b.astype(np.int64)
    pax, pay = ports.pa[:, 0].copy(), ports.pa[:, 1].copy()
    pbx, pby = ports.pb[:, 0].copy(), ports.pb[:, 1].copy()
    my = ports.my.astype(float)
    gp, gn = SequencePair.from_positions(init_x + ports.w / 2, init_y + ports.h / 2)
    gp, gn = np.array(gp, np.int64), np.array(gn, np.int64)
    fx, fy = np.zeros(n, np.bool_), np.zeros(n, np.bool_)
    unit = math.sqrt(float((ports.w * ports.h).mean()))
    total_area = float((W_ * H_).sum())

    def cost(gp, gn, fx, fy):
        x, y = pack(gp, gn, W_, H_)
        c, bc, length = metrics(x, y, ports.w, ports.h, a, b, pax, pay, pbx, pby, fx, fy, my, 10.0)
        W = float((x + W_).max())
        H = float((y + H_).max())
        aspect = math.log((W / H) / weights.target_aspect) ** 2
        return (weights.crossing * c + weights.box_crossing * bc + weights.length * length / unit
                + weights.area * W * H / total_area + weights.aspect * aspect)

    current = cost(gp, gn, fx, fy)
    best = (current, gp.copy(), gn.copy(), fx.copy(), fy.copy())
    if t_start is None:
        deltas = []
        for _ in range(200):
            i, j = rng.choice(n, 2, replace=False)
            g2 = gp.copy()
            g2[i], g2[j] = g2[j], g2[i]
            deltas.append(abs(cost(g2, gn, fx, fy) - current))
        t_start = (float(np.median(deltas)) or 1.0) / -math.log(0.3)
    accepted = 0
    for it in range(iterations):
        temp = t_start * (t_end / t_start) ** (it / max(1, iterations - 1))
        move = rng.random()
        ngp, ngn, nfx, nfy = gp, gn, fx, fy
        if move < 0.3:
            ngp = gp.copy()
            i, j = rng.choice(n, 2, replace=False)
            ngp[i], ngp[j] = ngp[j], ngp[i]
        elif move < 0.6:
            ngn = gn.copy()
            i, j = rng.choice(n, 2, replace=False)
            ngn[i], ngn[j] = ngn[j], ngn[i]
        elif move < (0.85 if flips else 1.0):
            # exchange the places of two boxes
            u, v = rng.choice(n, 2, replace=False)
            ngp, ngn = gp.copy(), gn.copy()
            iu, iv = np.flatnonzero(ngp == u)[0], np.flatnonzero(ngp == v)[0]
            ngp[iu], ngp[iv] = v, u
            iu, iv = np.flatnonzero(ngn == u)[0], np.flatnonzero(ngn == v)[0]
            ngn[iu], ngn[iv] = v, u
        else:
            k = rng.integers(n)
            if flips == "x" or rng.random() < 0.5:
                nfx = fx.copy()
                nfx[k] = not nfx[k]
            else:
                nfy = fy.copy()
                nfy[k] = not nfy[k]
        value = cost(ngp, ngn, nfx, nfy)
        if value <= current or rng.random() < math.exp((current - value) / temp):
            gp, gn, fx, fy, current = ngp, ngn, nfx, nfy, value
            accepted += 1
            if value < best[0]:
                best = (value, gp.copy(), gn.copy(), fx.copy(), fy.copy())
        if log is not None and it % max(1, iterations // 10) == 0:
            log(f"  anneal {it}: T={temp:.3g} cost={current:.1f} best={best[0]:.1f} accepted={accepted}")
    _, gp, gn, fx, fy = best
    x, y = pack(gp, gn, W_, H_)
    if not flips:
        fx, fy = np.zeros(n, bool), np.zeros(n, bool)
    result = evaluate_arrangement(ports, x, y, fx, fy)
    result["cost"] = best[0]
    return x, y, fx, fy, result


def apply_arrangement(dm: DraftMap, names, x, y, fx=None, fy=None) -> None:
    """Set the box positions (top left) and flip the drawings of subsystems (relative node positions)."""
    index = {k: i for i, k in enumerate(names)}
    for k, s in dm.subsystems.items():
        i = index[k]
        s.x, s.y = float(x[i]), float(y[i])
        s.flip_x = bool(fx[i]) if fx is not None else False
        s.flip_y = bool(fy[i]) if fy is not None else False
    for n in dm.nodes.values():
        s = dm.subsystems[n.subsystem]
        if s.flip_x:
            n.x = s.w - n.x
        if s.flip_y:
            # keep the title on top: mirror the content between title and bottom padding
            top = s.title_h / 2 + TITLE_GAP
            n.y = top + (s.h - PADDING) - n.y


def internal_metrics(dm: DraftMap) -> dict:
    """Crossings of the links inside the subsystems, overlapping nodes and the area of the drawings."""
    edges = defaultdict(list)
    for e in dm.edges:
        a, b = dm.nodes[e.enzyme], dm.nodes[e.met]
        if a.subsystem == b.subsystem:
            edges[a.subsystem].append((e.kind, (a.x, a.y), (b.x, b.y)))
    result = dict(crossings_main=0, crossings_all=0, overlapping_nodes=0, drawing_area=0.0)
    for s in dm.subsystems:
        members = dm.members(s)
        if not members:
            continue
        box = np.array([n.box for n in members])
        result["drawing_area"] += float((box[:, 2].max() - box[:, 0].min()) * (box[:, 3].max() - box[:, 1].min()))
        ov = ((box[:, None, 0] < box[None, :, 2]) & (box[None, :, 0] < box[:, None, 2])
              & (box[:, None, 1] < box[None, :, 3]) & (box[None, :, 1] < box[:, None, 3]))
        result["overlapping_nodes"] += int(np.triu(ov, 1).sum())
        for kinds, key in (({"main"}, "crossings_main"), ({"main", "cofactor"}, "crossings_all")):
            seg = [(p, q) for k, p, q in edges[s] if k in kinds]
            if seg:
                p1 = np.array([p for p, _ in seg], float)
                p2 = np.array([q for _, q in seg], float)
                result[key] += count_crossings(p1, p2)
    return result
