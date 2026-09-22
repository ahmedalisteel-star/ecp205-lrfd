"""Angles, double angles and tees - the open sections whose compressive strength
is governed by flexural-torsional buckling (ECP 205 cl. 4.3).

This module supplies their **geometry only**: area, centroid, principal axes,
the shear centre, J, Cw and the polar radius of gyration about the shear centre.
It implements no provision of the code. The buckling checks that use these
properties belong to cl. 4.3.

Catalogue values are **nominal**. Before a section goes on a drawing, confirm it
against the mill certificate or the supplier's table - rolling tolerances and
national variants exist, and a profile that is not actually available is not a
design.

Where the numbers come from
---------------------------
Angle dimensions (legs, t, root radius r1, toe radius r2) are EN 10056-1. Every
property used in design - A, centroid, Ix, Iy, Ixy, Iu, Iv - is computed from
that outline with the root and toe radii included (:mod:`ecp205._polygon`), so
the set is self-consistent (``Iu + Iv == Ix + Iy`` exactly). The catalogue's own
A, centroid, Ix and Iy are kept in the table only as witnesses:
:func:`ecp205.sections.verify` re-derives them from the geometry and reports any
disagreement over 3 %, so a typo in a dimension cannot pass silently.

J and Cw use the usual thin-walled open-section expressions (Seaburg and
Carter, *Torsional Analysis of Structural Steel Members*, 1997, which gives
section properties only - no design provisions are taken from it).

Axes
----
x is horizontal and y vertical. An :class:`Angle` has its heel at the origin, the
leg ``b`` vertical and the leg ``d`` horizontal - for a catalogue unequal angle
the long leg is vertical. u and v are the major and minor principal axes. For a
:class:`Tee` and a :class:`DoubleAngle` the y axis is the axis of symmetry, as in
cl. 4.3.1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from ._polygon import Props, arc, props
from .core import DENSITY
from .sections import ISection
from .sections import get as get_I

__all__ = ["Angle", "Tee", "DoubleAngle", "EQUAL_ANGLES", "UNEQUAL_ANGLES",
           "ANGLES", "get_angle", "search_angles", "tee_from", "double_angle",
           "verify_angles"]


# ---------------------------------------------------------------- catalogue
# EN 10056-1 equal angles. Dimensions in cm.
#                      a     t     r1    r2    | witnesses: A    e     Ix
# e is the centroid distance from the back of either leg; Ix is about the
# centroidal axis parallel to a leg (Ix = Iy).
EQUAL_ANGLES = {
    "L 25x25x3":      (2.5, 0.3, 0.35, 0.2, 1.42, 0.72, 0.8),
    "L 25x25x4":      (2.5, 0.4, 0.35, 0.2, 1.85, 0.76, 1.01),
    "L 25x25x5":      (2.5, 0.5, 0.35, 0.2, 2.26, 0.8, 1.2),
    "L 40x40x4":      (4, 0.4, 0.6, 0.3, 3.08, 1.12, 4.47),
    "L 40x40x5":      (4, 0.5, 0.6, 0.3, 3.79, 1.16, 5.43),
    "L 40x40x6":      (4, 0.6, 0.6, 0.3, 4.48, 1.2, 6.31),
    "L 50x50x5":      (5, 0.5, 0.7, 0.35, 4.8, 1.4, 11),
    "L 50x50x6":      (5, 0.6, 0.7, 0.35, 5.69, 1.45, 12.8),
    "L 50x50x8":      (5, 0.8, 0.7, 0.35, 7.41, 1.52, 16.3),
    "L 60x60x5":      (6, 0.5, 0.8, 0.4, 5.82, 1.64, 19.4),
    "L 60x60x6":      (6, 0.6, 0.8, 0.4, 6.91, 1.69, 22.8),
    "L 60x60x8":      (6, 0.8, 0.8, 0.4, 9.03, 1.77, 29.2),
    "L 60x60x10":     (6, 1, 0.8, 0.4, 11.1, 1.85, 34.9),
    "L 65x65x5":      (6.5, 0.5, 0.9, 0.45, 6.37, 1.77, 24.7),
    "L 65x65x6":      (6.5, 0.6, 0.9, 0.45, 7.53, 1.81, 29.2),
    "L 65x65x8":      (6.5, 0.8, 0.9, 0.45, 9.76, 1.88, 37.5),
    "L 65x65x10":     (6.5, 1, 0.9, 0.45, 12, 1.96, 45.1),
    "L 70x70x8":      (7, 0.8, 0.9, 0.45, 10.6, 2.01, 47.5),
    "L 70x70x10":     (7, 1, 0.9, 0.45, 13.1, 2.09, 57.2),
    "L 75x75x5":      (7.5, 0.5, 0.9, 0.45, 7.36, 1.99, 38.5),
    "L 75x75x6":      (7.5, 0.6, 0.9, 0.45, 8.75, 2.04, 45.6),
    "L 75x75x7":      (7.5, 0.7, 0.9, 0.45, 10.1, 2.09, 52.4),
    "L 75x75x8":      (7.5, 0.8, 0.9, 0.45, 11.5, 2.13, 58.9),
    "L 75x75x10":     (7.5, 1, 0.9, 0.45, 14.1, 2.21, 71.2),
    "L 75x75x12":     (7.5, 1.2, 0.9, 0.45, 16.7, 2.29, 82.6),
    "L 80x80x6":      (8, 0.6, 1, 0.5, 9.25, 2.17, 55.8),
    "L 80x80x8":      (8, 0.8, 1, 0.5, 12.3, 2.26, 72.2),
    "L 80x80x10":     (8, 1, 1, 0.5, 15.1, 2.34, 87.5),
    "L 90x90x6":      (9, 0.6, 1.1, 0.55, 10.6, 2.41, 80.3),
    "L 90x90x7":      (9, 0.7, 1.1, 0.55, 12.2, 2.45, 92.5),
    "L 90x90x8":      (9, 0.8, 1.1, 0.55, 13.9, 2.5, 104),
    "L 90x90x9":      (9, 0.9, 1.1, 0.55, 15.5, 2.54, 116),
    "L 90x90x10":     (9, 1, 1.1, 0.55, 17.1, 2.58, 127),
    "L 90x90x12":     (9, 1.2, 1.1, 0.55, 20.3, 2.66, 148),
    "L 100x100x6":    (10, 0.6, 1.2, 0.6, 11.8, 2.64, 111),
    "L 100x100x8":    (10, 0.8, 1.2, 0.6, 15.5, 2.74, 145),
    "L 100x100x10":   (10, 1, 1.2, 0.6, 19.2, 2.82, 177),
    "L 100x100x12":   (10, 1.2, 1.2, 0.6, 22.7, 2.9, 207),
    "L 100x100x15":   (10, 1.5, 1.2, 0.6, 27.9, 3.02, 249),
    "L 120x120x8":    (12, 0.8, 1.3, 0.65, 18.7, 3.23, 255),
    "L 120x120x10":   (12, 1, 1.3, 0.65, 23.2, 3.31, 313),
    "L 120x120x12":   (12, 1.2, 1.3, 0.65, 27.5, 3.4, 368),
    "L 120x120x14":   (12, 1.4, 1.3, 0.65, 31.8, 3.48, 420),
    "L 120x120x15":   (12, 1.5, 1.3, 0.65, 33.9, 3.51, 445),
    "L 130x130x10":   (13, 1, 1.4, 0.7, 25.2, 3.55, 401),
    "L 130x130x12":   (13, 1.2, 1.4, 0.7, 30, 3.64, 472),
    "L 130x130x15":   (13, 1.5, 1.4, 0.7, 37, 3.76, 573),
    "L 130x130x16":   (13, 1.6, 1.4, 0.7, 39.3, 3.8, 605),
    "L 150x150x10":   (15, 1, 1.6, 0.8, 29.14, 4.03, 624),
    "L 150x150x12":   (15, 1.2, 1.6, 0.8, 34.8, 4.12, 737),
    "L 150x150x15":   (15, 1.5, 1.6, 0.8, 43, 4.25, 898),
    "L 150x150x18":   (15, 1.8, 1.6, 0.8, 51, 4.37, 1050),
    "L 150x150x19":   (15, 1.9, 1.6, 0.8, 53.7, 4.4, 1100),
    "L 150x150x20":   (15, 2, 1.6, 0.8, 56.3, 4.44, 1150),
    "L 200x200x16":   (20, 1.6, 1.8, 0.9, 61.8, 5.52, 2340),
    "L 200x200x18":   (20, 1.8, 1.8, 0.9, 69.14, 5.6, 2600),
    "L 200x200x20":   (20, 2, 1.8, 0.9, 76.6, 5.68, 2850),
    "L 200x200x24":   (20, 2.4, 1.8, 0.9, 90.8, 5.84, 3330),
}

# EN 10056-1 unequal angles, long leg ``a`` first. Dimensions in cm.
#                      a     b     t     r1    r2    | witnesses: A  ea  eb  Ix  Iy
# ea is the centroid distance measured along the long leg (from the back of the
# short leg), eb along the short leg. Ix is about the centroidal axis parallel to
# the short leg (the larger value), Iy parallel to the long leg.
UNEQUAL_ANGLES = {
    "L 65x50x5":      (6.5, 5, 0.5, 0.7, 0.35, 5.54, 1.99, 1.25, 23.2, 11.9),
    "L 65x50x6":      (6.5, 5, 0.6, 0.7, 0.35, 6.58, 2.04, 1.29, 27.2, 14),
    "L 65x50x8":      (6.5, 5, 0.8, 0.7, 0.35, 8.6, 2.11, 1.37, 34.8, 17.7),
    "L 75x50x6":      (7.5, 5, 0.6, 0.8, 0.4, 7.19, 2.44, 1.21, 40.5, 14.4),
    "L 75x50x8":      (7.5, 5, 0.8, 0.8, 0.4, 9.41, 2.52, 1.29, 52, 18.4),
    "L 75x50x10":     (7.5, 5, 1, 0.8, 0.4, 11.6, 2.6, 1.36, 62.6, 21.9),
    "L 80x60x6":      (8, 6, 0.6, 0.9, 0.45, 8.11, 2.47, 1.48, 51.4, 24.8),
    "L 80x60x7":      (8, 6, 0.7, 0.9, 0.45, 9.38, 2.51, 1.52, 59, 28.4),
    "L 80x60x8":      (8, 6, 0.8, 0.9, 0.45, 10.6, 2.55, 1.56, 66.3, 31.8),
    "L 100x50x6":     (10, 5, 0.6, 1, 0.5, 8.73, 3.49, 1.04, 89.7, 15.3),
    "L 100x50x8":     (10, 5, 0.8, 1, 0.5, 11.4, 3.59, 1.12, 116, 19.5),
    "L 100x50x10":    (10, 5, 1, 1, 0.5, 14.1, 3.67, 1.2, 141, 23.4),
    "L 100x65x7":     (10, 6.5, 0.7, 1, 0.5, 11.2, 3.23, 1.51, 113, 37.6),
    "L 100x65x8":     (10, 6.5, 0.8, 1, 0.5, 12.7, 3.27, 1.55, 127, 42.2),
    "L 100x65x10":    (10, 6.5, 1, 1, 0.5, 15.6, 3.36, 1.63, 154, 51),
    "L 100x75x7":     (10, 7.5, 0.7, 1.1, 0.55, 11.9, 3.06, 1.83, 118, 56.9),
    "L 100x75x8":     (10, 7.5, 0.8, 1.1, 0.55, 13.5, 3.1, 1.87, 133, 64.1),
    "L 100x75x10":    (10, 7.5, 1, 1.1, 0.55, 16.6, 3.19, 1.95, 162, 77.6),
    "L 100x75x12":    (10, 7.5, 1.2, 1.1, 0.55, 19.7, 3.27, 2.03, 189, 90.2),
    "L 120x80x8":     (12, 8, 0.8, 1.1, 0.55, 15.5, 3.83, 1.87, 226, 80.8),
    "L 120x80x10":    (12, 8, 1, 1.1, 0.55, 19.1, 3.92, 1.95, 276, 98.1),
    "L 120x80x12":    (12, 8, 1.2, 1.1, 0.55, 22.7, 4, 2.03, 323, 114),
    "L 120x80x14":    (12, 8, 1.4, 1.1, 0.55, 26.2, 4.08, 2.1, 368, 130),
    "L 125x75x8":     (12.5, 7.5, 0.8, 1.1, 0.55, 15.5, 4.14, 1.68, 247, 67.6),
    "L 125x75x10":    (12.5, 7.5, 1, 1.1, 0.55, 19.1, 4.23, 1.76, 302, 82.1),
    "L 125x75x12":    (12.5, 7.5, 1.2, 1.1, 0.55, 22.7, 4.31, 1.84, 354, 95.5),
    "L 150x75x10":    (15, 7.5, 1, 1.1, 0.55, 21.6, 5.32, 1.61, 501, 85.8),
    "L 150x75x12":    (15, 7.5, 1.2, 1.1, 0.55, 25.7, 5.41, 1.69, 589, 99.9),
    "L 150x75x15":    (15, 7.5, 1.5, 1.1, 0.55, 31.6, 5.53, 1.81, 713, 120),
    "L 150x90x10":    (15, 9, 1, 1.2, 0.6, 23.2, 5, 2.04, 533, 146),
    "L 150x90x12":    (15, 9, 1.2, 1.2, 0.6, 27.5, 5.08, 2.12, 627, 171),
    "L 150x90x15":    (15, 9, 1.5, 1.2, 0.6, 33.9, 5.21, 2.23, 761, 205),
    "L 150x100x10":   (15, 10, 1, 1.2, 0.6, 24.2, 4.8, 2.34, 552, 198),
    "L 150x100x12":   (15, 10, 1.2, 1.2, 0.6, 28.7, 4.89, 2.42, 650, 232),
    "L 150x100x14":   (15, 10, 1.4, 1.2, 0.6, 33.2, 4.97, 2.5, 743, 264),
    "L 150x100x15":   (15, 10, 1.5, 1.2, 0.6, 35.4, 5.01, 2.54, 789, 280),
    "L 200x100x10":   (20, 10, 1, 1.4, 0.7, 29.2, 6.93, 2.01, 1220, 210),
    "L 200x100x12":   (20, 10, 1.2, 1.4, 0.7, 34.8, 7.03, 2.1, 1440, 247),
    "L 200x100x15":   (20, 10, 1.5, 1.4, 0.7, 43, 7.16, 2.22, 1758, 299),
    "L 200x150x12":   (20, 15, 1.2, 1.5, 0.75, 40.8, 6.08, 3.61, 1652, 803),
    "L 200x150x15":   (20, 15, 1.5, 1.5, 0.75, 50.5, 6.21, 3.73, 2022, 979),
    "L 200x150x18":   (20, 15, 1.8, 1.5, 0.75, 60, 6.33, 3.85, 2376, 1146),
}

ANGLES = {**EQUAL_ANGLES, **UNEQUAL_ANGLES}


def _thin_walled_torsion(parts: list[tuple[float, float]]) -> float:
    """``J = sum(b t^3) / 3`` over the plate elements ``(b, t)``."""
    return sum(b * t ** 3 for b, t in parts) / 3


# ---------------------------------------------------------------- angle
@dataclass
class Angle:
    """A single angle. All dimensions in **cm**.

    The heel is at the origin, leg ``b`` runs up the y axis and leg ``d`` along
    the x axis. Every property is computed from the outline, root radius ``r1``
    and toe radii ``r2`` included.

    >>> a = get_angle("L 100x100x10")
    >>> round(a.A, 2), round(a.ex, 2), round(a.Ix, 1)
    (19.15, 2.82, 176.7)
    >>> round(a.alpha, 1), round(a.rv, 2), round(a.ru, 2)
    (45.0, 1.95, 3.83)
    """

    name: str
    b: float             #: vertical leg (the long leg of a catalogue angle)
    d: float             #: horizontal leg
    t: float             #: thickness
    r1: float = 0.0      #: root radius
    r2: float = 0.0      #: toe radius
    A: float = field(init=False)      #: area, cm2
    ex: float = field(init=False)     #: centroid from the back of leg b
    ey: float = field(init=False)     #: centroid from the back of leg d
    Ix: float = field(init=False)     #: about the centroidal axis parallel to d
    Iy: float = field(init=False)     #: about the centroidal axis parallel to b
    Ixy: float = field(init=False)    #: product of inertia (negative here)
    Iu: float = field(init=False)     #: major principal axis
    Iv: float = field(init=False)     #: minor principal axis
    alpha: float = field(init=False)  #: angle from x to u, degrees
    J: float = field(init=False)      #: St Venant torsion constant, cm4
    Cw: float = field(init=False)     #: warping constant, cm6

    def __post_init__(self) -> None:
        b, d, t = self.b, self.d, self.t
        if not 0 < t < min(b, d):
            raise ValueError(f"{self.name}: need 0 < t < both legs")
        if self.r2 > t or t + self.r1 > min(b, d):
            raise ValueError(f"{self.name}: radii do not fit the legs")
        p = props(self.outline())
        self.A = p.A
        self.ex, self.ey = p.xc, p.yc
        self.Ix, self.Iy, self.Ixy = p.Ix, p.Iy, p.Ixy
        mean = (p.Ix + p.Iy) / 2
        rad = math.hypot((p.Ix - p.Iy) / 2, p.Ixy)
        self.Iu, self.Iv = mean + rad, mean - rad
        self.alpha = math.degrees(0.5 * math.atan2(-2 * p.Ixy, p.Ix - p.Iy))
        bm, dm = b - t / 2, d - t / 2          # leg lengths to the mid-line
        self.J = _thin_walled_torsion([(bm, t), (dm, t)])
        self.Cw = t ** 3 * (bm ** 3 + dm ** 3) / 36

    def outline(self) -> list[tuple[float, float]]:
        """Counter-clockwise outline, heel at the origin."""
        b, d, t, r1, r2 = self.b, self.d, self.t, self.r1, self.r2
        return ([(0.0, 0.0), (d, 0.0)]
                + arc(d - r2, t - r2, r2, 0, 90)          # toe of leg d
                + arc(t + r1, t + r1, r1, 270, 180)       # root
                + arc(t - r2, b - r2, r2, 0, 90)          # toe of leg b
                + [(0.0, b)])

    # ------------------------------------------------------------ derived
    @property
    def equal(self) -> bool:
        return self.b == self.d

    @property
    def rx(self) -> float:
        return math.sqrt(self.Ix / self.A)

    @property
    def ry(self) -> float:
        return math.sqrt(self.Iy / self.A)

    @property
    def ru(self) -> float:
        """Radius of gyration about the major principal axis."""
        return math.sqrt(self.Iu / self.A)

    @property
    def rv(self) -> float:
        """Radius of gyration about the minor principal axis - the least r."""
        return math.sqrt(self.Iv / self.A)

    @property
    def x0(self) -> float:
        """Shear centre relative to the centroid, x. The shear centre of an angle
        is at the intersection of the leg mid-lines, (t/2, t/2) from the heel."""
        return self.t / 2 - self.ex

    @property
    def y0(self) -> float:
        return self.t / 2 - self.ey

    @property
    def u0(self) -> float:
        """Shear centre relative to the centroid, along the major axis u."""
        a = math.radians(self.alpha)
        return self.x0 * math.cos(a) + self.y0 * math.sin(a)

    @property
    def v0(self) -> float:
        """Shear centre relative to the centroid, along the minor axis v."""
        a = math.radians(self.alpha)
        return -self.x0 * math.sin(a) + self.y0 * math.cos(a)

    @property
    def ro2(self) -> float:
        """Polar radius of gyration about the shear centre, squared, cm2::

            ro^2 = x0^2 + y0^2 + (Ix + Iy) / A
        """
        return self.x0 ** 2 + self.y0 ** 2 + (self.Ix + self.Iy) / self.A

    @property
    def mass(self) -> float:
        """Mass per metre, kg/m."""
        return self.A * DENSITY / 10

    def __str__(self) -> str:
        return (f"{self.name}: A={self.A:0.2f} cm2 ex={self.ex:0.2f} "
                f"ey={self.ey:0.2f} Iu={self.Iu:0.1f} Iv={self.Iv:0.1f} "
                f"rv={self.rv:0.2f} ({self.mass:0.2f} kg/m)")


def _angle_key(name: str) -> str:
    k = (name.upper().replace("*", "X").replace("×", "X")
         .replace(" ", "").lstrip("L"))
    for cand in ANGLES:
        if cand[2:].upper() == k:
            return cand
    raise KeyError(f"{name!r} is not in the angle catalogue; {len(ANGLES)} "
                   f"angles available, e.g. 'L 100x100x10', 'L 100x75x8'")


def get_angle(name: str) -> Angle:
    """Look a catalogue angle up by name.

    >>> get_angle("L100x75x8").name
    'L 100x75x8'
    >>> get_angle("l 60X60X6").equal
    True
    """
    key = _angle_key(name)
    row = ANGLES[key]
    if key in EQUAL_ANGLES:
        a, t, r1, r2 = row[:4]
        return Angle(key, b=a, d=a, t=t, r1=r1, r2=r2)
    a, b, t, r1, r2 = row[:5]
    return Angle(key, b=a, d=b, t=t, r1=r1, r2=r2)


def search_angles(A_min: float = 0.0, rv_min: float = 0.0,
                  equal: Optional[bool] = None) -> list[dict]:
    """Catalogue angles meeting the criteria, **lightest first**.

    Args:
        A_min: minimum area, cm2.
        rv_min: minimum least radius of gyration, cm.
        equal: True for equal angles only, False for unequal only, None for both.

    Returns:
        A list of dicts with ``name``, ``mass``, ``A`` and ``rv``.

    >>> search_angles(A_min=12, rv_min=1.9, equal=True)[0]["name"]
    'L 100x100x8'
    """
    out = []
    for name in ANGLES:
        a = get_angle(name)
        if equal is not None and a.equal != equal:
            continue
        if a.A >= A_min and a.rv >= rv_min:
            out.append({"name": name, "mass": a.mass, "A": a.A, "rv": a.rv})
    return sorted(out, key=lambda row: row["mass"])


def verify_angles(tol: float = 0.03) -> list[str]:
    """Re-derive A, centroid, Ix and Iy of every catalogue angle from its
    outline and report disagreements with the tabulated witnesses.

    >>> verify_angles()
    []
    """
    bad = []
    for name, row in ANGLES.items():
        a = get_angle(name)
        if name in EQUAL_ANGLES:
            A, e, Ix = row[4:]
            given = {"A": A, "e": e, "Ix": Ix}
            calc = {"A": a.A, "e": a.ex, "Ix": a.Ix}
        else:
            A, ea, eb, Ix, Iy = row[5:]
            given = {"A": A, "ea": ea, "eb": eb, "Ix": Ix, "Iy": Iy}
            calc = {"A": a.A, "ea": a.ey, "eb": a.ex, "Ix": a.Ix, "Iy": a.Iy}
        for label, g in given.items():
            err = abs(g - calc[label]) / g
            if err > tol:
                bad.append(f"{name:<13s} {label:<3s} table={g:<10.4g} "
                           f"geometry={calc[label]:<10.4g} ({err * 100:4.1f}% off)")
    return bad


# ---------------------------------------------------------------- tee
@dataclass
class Tee:
    """A T section, flange at the bottom, stem up the y axis (the axis of
    symmetry). All dimensions in **cm**.

    ``J`` and ``Cw`` default to the thin-walled expressions; :func:`tee_from`
    takes J as half the parent's tabulated value instead, which keeps the
    fillet contribution.

    >>> t = tee_from("IPE 400")
    >>> t.name, t.d, round(t.A, 2)
    ('1/2 IPE 400', 20.0, 42.23)
    >>> round(t.yc, 2), round(t.y0, 2)
    (4.52, -3.85)
    """

    name: str
    d: float                   #: overall depth, flange face to stem tip
    bf: float                  #: flange width
    tf: float                  #: flange thickness
    tw: float                  #: stem thickness
    r: float = 0.0             #: root radius; 0 for a welded tee
    rolled: bool = True        #: rolled, or cut from a rolled section
    J: Optional[float] = None
    Cw: Optional[float] = None
    A: float = field(init=False)    #: area, cm2
    yc: float = field(init=False)   #: centroid from the outer face of the flange
    Ix: float = field(init=False)   #: about the centroidal axis parallel to the flange
    Iy: float = field(init=False)   #: about the axis of symmetry

    def __post_init__(self) -> None:
        if not (0 < self.tf < self.d and 0 < self.tw < self.bf):
            raise ValueError(f"{self.name}: inconsistent plate dimensions")
        p: Props = props(self.outline())
        self.A, self.yc, self.Ix, self.Iy = p.A, p.yc, p.Ix, p.Iy
        hs = self.d - self.tf / 2              # stem length to flange mid-line
        if self.J is None:
            self.J = _thin_walled_torsion([(self.bf, self.tf), (hs, self.tw)])
        if self.Cw is None:
            self.Cw = (self.tf ** 3 * self.bf ** 3 / 144
                       + self.tw ** 3 * hs ** 3 / 36)

    def outline(self) -> list[tuple[float, float]]:
        """Counter-clockwise outline, flange underside centred on the origin."""
        bf, tf, tw, d, r = self.bf, self.tf, self.tw, self.d, self.r
        return ([(-bf / 2, 0.0), (bf / 2, 0.0), (bf / 2, tf)]
                + arc(tw / 2 + r, tf + r, r, 270, 180)
                + [(tw / 2, d), (-tw / 2, d)]
                + arc(-tw / 2 - r, tf + r, r, 0, -90)
                + [(-bf / 2, tf)])

    @property
    def rx(self) -> float:
        return math.sqrt(self.Ix / self.A)

    @property
    def ry(self) -> float:
        """About the axis of symmetry."""
        return math.sqrt(self.Iy / self.A)

    @property
    def x0(self) -> float:
        return 0.0

    @property
    def y0(self) -> float:
        """Shear centre relative to the centroid, on the axis of symmetry. The
        shear centre is at the flange mid-thickness, so y0 is negative."""
        return self.tf / 2 - self.yc

    @property
    def ro2(self) -> float:
        """Polar radius of gyration about the shear centre, squared, cm2."""
        return self.y0 ** 2 + (self.Ix + self.Iy) / self.A

    @property
    def mass(self) -> float:
        return self.A * DENSITY / 10

    def __str__(self) -> str:
        return (f"{self.name}: A={self.A:0.2f} cm2 yc={self.yc:0.2f} "
                f"Ix={self.Ix:0.1f} Iy={self.Iy:0.1f} y0={self.y0:0.2f} "
                f"({self.mass:0.2f} kg/m)")


def tee_from(parent: str | ISection) -> Tee:
    """A tee cut at mid-depth from an I section - ``1/2 IPE 400`` and so on.

    J is half the parent's (tabulated) value.
    """
    p = get_I(parent) if isinstance(parent, str) else parent
    return Tee(f"1/2 {p.name}", d=p.d / 2, bf=p.bf, tf=p.tf, tw=p.tw, r=p.r,
               rolled=p.rolled, J=p.J / 2)


# ---------------------------------------------------------------- double angle
@dataclass
class DoubleAngle:
    """Two identical angles back to back with their vertical legs separated by
    ``gap`` (the gusset thickness). y is the axis of symmetry; the outstanding
    legs are at the bottom.

    The pair is treated as one section; the shear centre is taken at the
    mid-thickness of the outstanding legs on the axis of symmetry. ``Cw`` is set
    to zero, which is conservative.

    >>> p = double_angle("L 100x100x10", gap=1.0)
    >>> round(p.A, 2), round(p.Iy, 1), round(p.y0, 2)
    (38.31, 776.2, -2.32)
    """

    angle: Angle
    gap: float = 0.0
    name: str = ""

    def __post_init__(self) -> None:
        if self.gap < 0:
            raise ValueError("gap must be >= 0")
        if not self.name:
            self.name = f"2{self.angle.name} gap {self.gap:g}"

    @property
    def A(self) -> float:
        return 2 * self.angle.A

    @property
    def yc(self) -> float:
        """Centroid from the underside of the outstanding legs."""
        return self.angle.ey

    @property
    def Ix(self) -> float:
        return 2 * self.angle.Ix

    @property
    def Iy(self) -> float:
        """About the axis of symmetry."""
        a = self.angle
        return 2 * (a.Iy + a.A * (self.gap / 2 + a.ex) ** 2)

    @property
    def J(self) -> float:
        return 2 * self.angle.J

    @property
    def Cw(self) -> float:
        return 0.0

    @property
    def rx(self) -> float:
        return math.sqrt(self.Ix / self.A)

    @property
    def ry(self) -> float:
        return math.sqrt(self.Iy / self.A)

    @property
    def x0(self) -> float:
        return 0.0

    @property
    def y0(self) -> float:
        return self.angle.t / 2 - self.angle.ey

    @property
    def ro2(self) -> float:
        return self.y0 ** 2 + (self.Ix + self.Iy) / self.A

    @property
    def mass(self) -> float:
        return self.A * DENSITY / 10

    def __str__(self) -> str:
        return (f"{self.name}: A={self.A:0.2f} cm2 Ix={self.Ix:0.1f} "
                f"Iy={self.Iy:0.1f} rx={self.rx:0.2f} ry={self.ry:0.2f} "
                f"({self.mass:0.2f} kg/m)")


def double_angle(name: str, gap: float = 0.0, legs: str = "long") -> DoubleAngle:
    """Two catalogue angles back to back.

    Args:
        name: the angle, e.g. ``"L 100x75x8"``.
        gap: clear distance between the back-to-back legs, cm - the gusset
            thickness.
        legs: which legs are back to back for an unequal angle, ``"long"`` or
            ``"short"``.
    """
    a = get_angle(name)
    if legs not in ("long", "short"):
        raise ValueError("legs must be 'long' or 'short'")
    tag = ""
    if not a.equal:
        tag = " LLBB" if legs == "long" else " SLBB"
        if legs == "short":
            a = Angle(a.name, b=a.d, d=a.b, t=a.t, r1=a.r1, r2=a.r2)
    return DoubleAngle(a, gap, f"2{a.name}{tag} gap {gap:g}")
