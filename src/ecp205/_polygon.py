"""Section properties of a closed polygon outline - the geometry engine behind
:mod:`ecp205.open_sections`.

A rolled section is described by its outline with the root and toe radii
replaced by short chords (:func:`arc`). Green's theorem then gives the area,
centroid and second moments exactly for that polygon. With 32 chords per
quarter circle the chord error is below 1e-4 of the fillet area - far inside
rolling tolerance.

This is pure geometry. It contains no provision of ECP 205.
"""

from __future__ import annotations

import math
from typing import NamedTuple

__all__ = ["Props", "arc", "props"]

Point = tuple[float, float]

N_ARC = 32   #: chords per quarter circle


class Props(NamedTuple):
    """Area and centroidal second moments of a plane figure. cm units."""

    A: float     #: area, cm2
    xc: float    #: centroid, cm
    yc: float
    Ix: float    #: about the horizontal axis through the centroid, cm4
    Iy: float    #: about the vertical axis through the centroid, cm4
    Ixy: float   #: product of inertia about the centroidal axes, cm4


def arc(cx: float, cy: float, r: float, a0: float, a1: float,
        n: int = N_ARC) -> list[Point]:
    """Points on a circular arc from angle ``a0`` to ``a1`` (degrees), inclusive.

    A zero radius returns the single corner point, so a section with sharp
    corners needs no special case.

    >>> arc(0, 0, 0, 0, 90)
    [(0, 0)]
    >>> [tuple(round(v, 9) for v in p) for p in arc(0, 0, 1, 0, 90, n=2)]
    [(1.0, 0.0), (0.707106781, 0.707106781), (0.0, 1.0)]
    """
    if r == 0:
        return [(cx, cy)]
    steps = max(1, round(n * abs(a1 - a0) / 90))
    out = []
    for i in range(steps + 1):
        a = math.radians(a0 + (a1 - a0) * i / steps)
        out.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return out


def props(pts: list[Point]) -> Props:
    """Properties of a simple polygon given counter-clockwise.

    >>> p = props([(0, 0), (4, 0), (4, 2), (0, 2)])     # 4 x 2 rectangle
    >>> p.A, p.xc, p.yc, round(p.Ix, 6), round(p.Iy, 6), round(p.Ixy, 9)
    (8.0, 2.0, 1.0, 2.666667, 10.666667, 0.0)
    """
    A = Sx = Sy = Ixx = Iyy = Pxy = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        c = x0 * y1 - x1 * y0
        A += c
        Sx += (x0 + x1) * c
        Sy += (y0 + y1) * c
        Ixx += (y0 * y0 + y0 * y1 + y1 * y1) * c
        Iyy += (x0 * x0 + x0 * x1 + x1 * x1) * c
        Pxy += (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * c
    A /= 2
    if A <= 0:
        raise ValueError("outline must be a counter-clockwise simple polygon")
    xc, yc = Sx / (6 * A), Sy / (6 * A)
    return Props(A, xc, yc,
                 Ixx / 12 - A * yc * yc,
                 Iyy / 12 - A * xc * xc,
                 Pxy / 24 - A * xc * yc)
