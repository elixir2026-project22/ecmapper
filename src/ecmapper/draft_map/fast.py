"""numba kernels for the annealing of arrangements."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def pack(gp, gn, w, h):
    n = len(gp)
    rank = np.empty(n, np.int64)
    for i in range(n):
        rank[gn[i]] = i
    x = np.zeros(n)
    y = np.zeros(n)
    for k in range(1, n):
        b = gp[k]
        xb = 0.0
        yb = 0.0
        for j in range(k):
            a = gp[j]
            if rank[a] < rank[b]:
                v = x[a] + w[a]
                if v > xb:
                    xb = v
            else:
                v = y[a] + h[a]
                if v > yb:
                    yb = v
        x[b] = xb
        y[b] = yb
    return x, y


@njit(cache=True)
def _orient(ox, oy, ax, ay, bx, by):
    return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)


@njit(cache=True)
def metrics(x, y, w, h, a, b, pax, pay, pbx, pby, fx, fy, my, shrink):
    """crossings, box crossings and total length of the links for box positions and flips."""
    L = len(a)
    n = len(x)
    x1 = np.empty(L)
    y1 = np.empty(L)
    x2 = np.empty(L)
    y2 = np.empty(L)
    for i in range(L):
        px = w[a[i]] - pax[i] if fx[a[i]] else pax[i]
        py = my[a[i]] - pay[i] if fy[a[i]] else pay[i]
        qx = w[b[i]] - pbx[i] if fx[b[i]] else pbx[i]
        qy = my[b[i]] - pby[i] if fy[b[i]] else pby[i]
        x1[i] = x[a[i]] + px
        y1[i] = y[a[i]] + py
        x2[i] = x[b[i]] + qx
        y2[i] = y[b[i]] + qy
    crossings = 0
    for i in range(L):
        for j in range(i + 1, L):
            o1 = _orient(x1[i], y1[i], x2[i], y2[i], x1[j], y1[j])
            o2 = _orient(x1[i], y1[i], x2[i], y2[i], x2[j], y2[j])
            if o1 * o2 >= 0:
                continue
            o3 = _orient(x1[j], y1[j], x2[j], y2[j], x1[i], y1[i])
            o4 = _orient(x1[j], y1[j], x2[j], y2[j], x2[i], y2[i])
            if o3 * o4 < 0:
                crossings += 1
    box_crossings = 0
    length = 0.0
    for i in range(L):
        dx = x2[i] - x1[i]
        dy = y2[i] - y1[i]
        length += np.sqrt(dx * dx + dy * dy)
        for k in range(n):
            if k == a[i] or k == b[i]:
                continue
            bx0 = x[k] + shrink
            bx1 = x[k] + w[k] - shrink
            by0 = y[k] + shrink
            by1 = y[k] + h[k] - shrink
            t0 = 0.0
            t1 = 1.0
            ok = True
            for side in range(4):
                if side == 0:
                    p = -dx
                    q = x1[i] - bx0
                elif side == 1:
                    p = dx
                    q = bx1 - x1[i]
                elif side == 2:
                    p = -dy
                    q = y1[i] - by0
                else:
                    p = dy
                    q = by1 - y1[i]
                if abs(p) < 1e-12:
                    if q < 0:
                        ok = False
                        break
                else:
                    r = q / p
                    if p < 0:
                        if r > t0:
                            t0 = r
                    else:
                        if r < t1:
                            t1 = r
            if ok and t0 < t1:
                box_crossings += 1
    return crossings, box_crossings, length
