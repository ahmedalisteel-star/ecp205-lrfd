"""Cross-sections: the :class:`ISection` geometry object and a catalogue of the
European hot-rolled profiles (IPE, HEA, HEB) used in Egyptian practice.

Catalogue values are **nominal**. Before a section goes on a drawing, confirm it
against the mill certificate or the supplier's table - rolling tolerances and
national variants exist, and a profile that is not actually available is not a
design. :func:`verify` re-derives A, Ix and Zx from the plate geometry so a typo
in the table cannot pass silently.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from .core import DENSITY

__all__ = ["ISection", "IPE", "HEA", "HEB", "TABLES", "CATALOGUE",
           "get", "search", "weight", "verify"]


@dataclass
class ISection:
    """A doubly symmetric I or H section. All dimensions in **cm**.

    Section properties left as ``None`` are computed from the plate geometry,
    including the four root fillets (total area ``(4 - pi) r^2`` taken at the
    mid-depth of the web). Including the fillets brings a rolled profile's
    computed Ix and Zx to within about 1 % of the catalogue value; omitting them
    under-reads by roughly 5 %.

    >>> s = ISection("IPE 400", d=40, bf=18, tf=1.35, tw=0.86, r=2.1)
    >>> round(s.Sx)
    1160
    """

    name: str
    d: float                    #: overall depth
    bf: float                   #: flange width
    tf: float                   #: flange thickness
    tw: float                   #: web thickness
    r: float = 0.0              #: root radius; 0 for a welded section
    rolled: bool = True         #: rolled (True) or welded (False) - Table 2.12c
    A: Optional[float] = None
    Ix: Optional[float] = None
    Iy: Optional[float] = None
    Zx: Optional[float] = None  #: plastic section modulus about x
    Zy: Optional[float] = None
    J: Optional[float] = None
    Cw: Optional[float] = None

    def __post_init__(self) -> None:
        hw = self.d - 2 * self.tf
        Ar = (4 - math.pi) * self.r ** 2          # four root fillets
        if self.A is None:
            self.A = 2 * self.bf * self.tf + hw * self.tw + Ar
        if self.Ix is None:
            self.Ix = ((self.bf * self.d ** 3 - (self.bf - self.tw) * hw ** 3) / 12
                       + Ar * (hw / 2) ** 2)
        if self.Iy is None:
            self.Iy = (2 * self.tf * self.bf ** 3 + hw * self.tw ** 3) / 12
        if self.Zx is None:
            self.Zx = (self.bf * self.tf * (self.d - self.tf)
                       + self.tw * hw ** 2 / 4 + Ar * hw / 2)
        if self.Zy is None:
            self.Zy = self.tf * self.bf ** 2 / 2 + hw * self.tw ** 2 / 4
        if self.J is None:
            self.J = (2 * self.bf * self.tf ** 3
                      + (self.d - self.tf) * self.tw ** 3) / 3
        if self.Cw is None:
            self.Cw = self.Iy * (self.d - self.tf) ** 2 / 4

    # ------------------------------------------------------------ derived
    @property
    def hw(self) -> float:
        """Clear web depth between the flanges."""
        return self.d - 2 * self.tf

    @property
    def dw(self) -> float:
        """Web depth for classification - clear depth less the fillets (Table 2.12a)."""
        return self.hw - 2 * self.r

    @property
    def c(self) -> float:
        """Outstand flange width (Table 2.12c)."""
        return (self.bf - self.tw) / 2 - self.r

    @property
    def Sx(self) -> float:
        """Elastic section modulus about the major axis."""
        return self.Ix / (self.d / 2)

    @property
    def Sy(self) -> float:
        return self.Iy / (self.bf / 2)

    @property
    def rx(self) -> float:
        return math.sqrt(self.Ix / self.A)

    @property
    def ry(self) -> float:
        return math.sqrt(self.Iy / self.A)

    @property
    def rT(self) -> float:
        """Radius of gyration of the compression flange plus one sixth of the web
        area, about the minor axis (cl. 5.1.3.1)."""
        Af = self.bf * self.tf
        Aw6 = self.hw * self.tw / 6
        I = self.tf * self.bf ** 3 / 12 + Aw6 * self.tw ** 2 / 12
        return math.sqrt(I / (Af + Aw6))

    @property
    def Aw(self) -> float:
        """Shear area - **overall depth** times web thickness (cl. 5.2.1)."""
        return self.d * self.tw

    @property
    def Af(self) -> float:
        """Area of one flange."""
        return self.bf * self.tf

    @property
    def mass(self) -> float:
        """Mass per metre, kg/m."""
        return self.A * DENSITY / 10

    def __str__(self) -> str:
        return (f"{self.name}: d={self.d} bf={self.bf} tf={self.tf} tw={self.tw} "
                f"| A={self.A:0.1f} cm2 Ix={self.Ix:0.0f} Sx={self.Sx:0.0f} "
                f"Zx={self.Zx:0.0f} ry={self.ry:0.2f} rT={self.rT:0.2f} "
                f"({self.mass:0.1f} kg/m)")


#            d      bf    tw     tf     r     A      Ix       Iy      Zx      J
IPE = {
    "IPE 100": (10.0,  5.5, 0.41, 0.57, 0.7,  10.3,    171,   15.9,   39.4,  1.20),
    "IPE 120": (12.0,  6.4, 0.44, 0.63, 0.7,  13.2,    318,   27.7,   60.7,  1.74),
    "IPE 140": (14.0,  7.3, 0.47, 0.69, 0.7,  16.4,    541,   44.9,   88.3,  2.45),
    "IPE 160": (16.0,  8.2, 0.50, 0.74, 0.9,  20.1,    869,   68.3,    124,  3.60),
    "IPE 180": (18.0,  9.1, 0.53, 0.80, 0.9,  23.9,   1317,    101,    166,  4.79),
    "IPE 200": (20.0, 10.0, 0.56, 0.85, 1.2,  28.5,   1943,    142,    221,  6.98),
    "IPE 220": (22.0, 11.0, 0.59, 0.92, 1.2,  33.4,   2772,    205,    285,  9.07),
    "IPE 240": (24.0, 12.0, 0.62, 0.98, 1.5,  39.1,   3892,    284,    367,  12.9),
    "IPE 270": (27.0, 13.5, 0.66, 1.02, 1.5,  45.9,   5790,    420,    484,  15.9),
    "IPE 300": (30.0, 15.0, 0.71, 1.07, 1.5,  53.8,   8356,    604,    628,  20.1),
    "IPE 330": (33.0, 16.0, 0.75, 1.15, 1.8,  62.6,  11770,    788,    804,  28.2),
    "IPE 360": (36.0, 17.0, 0.80, 1.27, 1.8,  72.7,  16270,   1043,   1019,  37.3),
    "IPE 400": (40.0, 18.0, 0.86, 1.35, 2.1,  84.5,  23130,   1318,   1307,  51.1),
    "IPE 450": (45.0, 19.0, 0.94, 1.46, 2.1,  98.8,  33740,   1676,   1702,  66.9),
    "IPE 500": (50.0, 20.0, 1.02, 1.60, 2.1, 116.0,  48200,   2142,   2194,  89.3),
    "IPE 550": (55.0, 21.0, 1.11, 1.72, 2.4, 134.0,  67120,   2668,   2787, 123.0),
    "IPE 600": (60.0, 22.0, 1.20, 1.90, 2.4, 156.0,  92080,   3387,   3512, 165.0),
}

HEA = {
    "HEA 100": ( 9.6, 10.0, 0.50, 0.80, 1.2,  21.2,  349.2,  133.8,   83.0,  5.24),
    "HEA 120": (11.4, 12.0, 0.50, 0.80, 1.2,  25.3,  606.2,  230.9,  119.5,  5.99),
    "HEA 140": (13.3, 14.0, 0.55, 0.85, 1.2,  31.4,   1033,  389.3,  173.5,  8.13),
    "HEA 160": (15.2, 16.0, 0.60, 0.90, 1.5,  38.8,   1673,  615.6,  245.1,  12.2),
    "HEA 180": (17.1, 18.0, 0.60, 0.95, 1.5,  45.3,   2510,  924.6,  324.9,  14.8),
    "HEA 200": (19.0, 20.0, 0.65, 1.00, 1.8,  53.8,   3692,   1336,  429.5,  21.0),
    "HEA 220": (21.0, 22.0, 0.70, 1.10, 1.8,  64.3,   5410,   1955,  568.5,  28.5),
    "HEA 240": (23.0, 24.0, 0.75, 1.20, 2.1,  76.8,   7763,   2769,  744.6,  41.6),
    "HEA 260": (25.0, 26.0, 0.75, 1.25, 2.4,  86.8,  10450,   3668,  919.8,  52.4),
    "HEA 280": (27.0, 28.0, 0.80, 1.30, 2.4,  97.3,  13670,   4763,   1112,  62.1),
    "HEA 300": (29.0, 30.0, 0.85, 1.40, 2.7, 112.5,  18260,   6310,   1383,  85.2),
    "HEA 320": (31.0, 30.0, 0.90, 1.55, 2.7, 124.4,  22930,   6985,   1628, 108.0),
    "HEA 340": (33.0, 30.0, 0.95, 1.65, 2.7, 133.5,  27690,   7436,   1850, 127.0),
    "HEA 360": (35.0, 30.0, 1.00, 1.75, 2.7, 142.8,  33090,   7887,   2088, 149.0),
    "HEA 400": (39.0, 30.0, 1.10, 1.90, 2.7, 159.0,  45070,   8564,   2562, 189.0),
    "HEA 450": (44.0, 30.0, 1.15, 2.10, 2.7, 178.0,  63720,   9465,   3216, 244.0),
    "HEA 500": (49.0, 30.0, 1.20, 2.30, 2.7, 197.5,  86970,  10370,   3949, 309.0),
    "HEA 600": (59.0, 30.0, 1.30, 2.50, 2.7, 226.5, 141200,  11270,   5350, 398.0),
}

HEB = {
    "HEB 100": (10.0, 10.0, 0.60, 1.00, 1.2,  26.0,  449.5,  167.3,  104.2,  9.25),
    "HEB 120": (12.0, 12.0, 0.65, 1.10, 1.2,  34.0,    864,  317.5,  165.2,  13.8),
    "HEB 140": (14.0, 14.0, 0.70, 1.20, 1.2,  43.0,   1509,  549.7,  245.4,  20.1),
    "HEB 160": (16.0, 16.0, 0.80, 1.30, 1.5,  54.3,   2492,  889.2,  354.0,  31.2),
    "HEB 180": (18.0, 18.0, 0.85, 1.40, 1.5,  65.3,   3831,   1363,  481.4,  42.2),
    "HEB 200": (20.0, 20.0, 0.90, 1.50, 1.8,  78.1,   5696,   2003,  642.5,  59.3),
    "HEB 220": (22.0, 22.0, 0.95, 1.60, 1.8,  91.0,   8091,   2843,  827.0,  76.6),
    "HEB 240": (24.0, 24.0, 1.00, 1.70, 2.1, 106.0,  11260,   3923,   1053, 102.7),
    "HEB 260": (26.0, 26.0, 1.00, 1.75, 2.4, 118.4,  14920,   5135,   1283, 123.8),
    "HEB 280": (28.0, 28.0, 1.05, 1.80, 2.4, 131.4,  19270,   6595,   1534, 143.7),
    "HEB 300": (30.0, 30.0, 1.10, 1.90, 2.7, 149.1,  25170,   8563,   1869, 185.0),
    "HEB 320": (32.0, 30.0, 1.15, 2.05, 2.7, 161.3,  30820,   9239,   2149, 225.0),
    "HEB 340": (34.0, 30.0, 1.20, 2.15, 2.7, 170.9,  36660,   9690,   2408, 257.0),
    "HEB 360": (36.0, 30.0, 1.25, 2.25, 2.7, 180.6,  43190,  10140,   2683, 292.0),
    "HEB 400": (40.0, 30.0, 1.35, 2.40, 2.7, 197.8,  57680,  10820,   3232, 356.0),
    "HEB 450": (45.0, 30.0, 1.40, 2.60, 2.7, 218.0,  79890,  11720,   3982, 440.0),
    "HEB 500": (50.0, 30.0, 1.45, 2.80, 2.7, 238.6, 107200,  12620,   4815, 538.0),
    "HEB 600": (60.0, 30.0, 1.55, 3.00, 2.7, 270.0, 171000,  13530,   6425, 667.0),
}

TABLES = {"IPE": IPE, "HEA": HEA, "HEB": HEB}
CATALOGUE = {**IPE, **HEA, **HEB}


def _key(name: str) -> str:
    k = name.upper().replace("-", " ").replace("_", " ").strip()
    if k in CATALOGUE:
        return k
    squashed = k.replace(" ", "")
    for cand in CATALOGUE:
        if cand.replace(" ", "") == squashed:
            return cand
    raise KeyError(f"{name!r} is not in the catalogue; "
                   f"{len(CATALOGUE)} profiles available, e.g. 'IPE 400', 'HEB 300'")


def get(name: str) -> ISection:
    """Look a catalogue profile up by name.

    >>> get("IPE 400").Zx
    1307
    >>> get("heb300").name
    'HEB 300'
    """
    d, bf, tw, tf, r, A, Ix, Iy, Zx, J = CATALOGUE[_key(name)]
    return ISection(_key(name), d=d, bf=bf, tw=tw, tf=tf, r=r, rolled=True,
                    A=A, Ix=Ix, Iy=Iy, Zx=Zx, J=J)


def weight(name: str) -> float:
    """Mass per metre, kg/m.

    >>> round(weight("IPE 400"), 1)
    66.3
    """
    return CATALOGUE[_key(name)][5] * DENSITY / 10


def search(family: Optional[str] = None, Zx_min: float = 0.0, Sx_min: float = 0.0,
           Ix_min: float = 0.0, d_max: float = 1e9,
           A_max: float = 1e9) -> list[dict]:
    """Catalogue profiles meeting the criteria, **lightest first**.

    Args:
        family: ``"IPE"``, ``"HEA"``, ``"HEB"``, or None for all.
        Zx_min: minimum plastic modulus, cm3.
        Sx_min: minimum elastic modulus, cm3.
        Ix_min: minimum second moment of area, cm4 - the deflection criterion.
        d_max: maximum overall depth, cm.
        A_max: maximum area, cm2.

    Returns:
        A list of dicts with ``name``, ``mass``, ``d``, ``A``, ``Ix``, ``Sx``, ``Zx``.

    >>> search(family="IPE", Zx_min=1500)[0]["name"]
    'IPE 450'
    """
    table = TABLES[family.upper()] if family else CATALOGUE
    out = []
    for name, row in table.items():
        d, A, Ix, Zx = row[0], row[5], row[6], row[8]
        Sx = Ix / (d / 2)
        if (Zx >= Zx_min and Sx >= Sx_min and Ix >= Ix_min
                and d <= d_max and A <= A_max):
            out.append({"name": name, "mass": A * DENSITY / 10, "d": d,
                        "A": A, "Ix": Ix, "Sx": Sx, "Zx": Zx})
    return sorted(out, key=lambda row: row["mass"])


def verify(tol: float = 0.03) -> list[str]:
    """Re-derive A, Ix and Zx from the plate geometry and report disagreements.

    Agreement is normally within about 1 %, so anything reported here is a
    suspect table entry rather than rolling tolerance.

    >>> verify()
    []
    """
    bad = []
    for name, row in CATALOGUE.items():
        d, bf, tw, tf, r, A, Ix, Zx = (row[0], row[1], row[2], row[3],
                                       row[4], row[5], row[6], row[8])
        hw = d - 2 * tf
        Ar = (4 - math.pi) * r ** 2
        calc = {
            "A": 2 * bf * tf + hw * tw + Ar,
            "Ix": (bf * d ** 3 - (bf - tw) * hw ** 3) / 12 + Ar * (hw / 2) ** 2,
            "Zx": bf * tf * (d - tf) + tw * hw ** 2 / 4 + Ar * hw / 2,
        }
        for label, given in (("A", A), ("Ix", Ix), ("Zx", Zx)):
            err = abs(given - calc[label]) / given
            if err > tol:
                bad.append(f"{name:<9s} {label:<3s} table={given:<10.4g} "
                           f"geometry={calc[label]:<10.4g} ({err * 100:4.1f}% off)")
    return bad
