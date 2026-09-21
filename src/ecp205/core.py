"""Constants, resistance factors and the result container shared by every module.

Units throughout the package: **ton, centimetre, t/cm2**. Moments in t.cm.
This is the unit system the Egyptian code is written in, and every tabulated
constant in it (127/sqrt(Fy), 80*ry/sqrt(Fy), 1380*Af/(d*Lb), 20700, 9500, ...)
is only valid there.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "E", "G", "NU", "DENSITY", "ALPHA_T", "PHI", "LAMBDA_MAX", "Check",
    "t_to_kN", "tcm2_to_MPa", "tcm_to_kNm",
]

# ---------------------------------------------------------------- cl. 1.3.2
E = 2100.0        #: Young's modulus, t/cm2
G = 810.0         #: shear modulus, t/cm2
NU = 0.3          #: Poisson's ratio
DENSITY = 7.85    #: mass density, t/m3
ALPHA_T = 1.2e-5  #: coefficient of thermal expansion, 1/degC

#: Resistance factors, by limit state, with the clause each comes from.
PHI = {
    "tension_yield": 0.85,      # cl. 3.1.1a
    "tension_fracture": 0.70,   # cl. 3.1.1b
    "compression": 0.80,        # cl. 4.2.1
    "flexure": 0.85,            # cl. 5.1.2
    "shear": 0.85,              # cl. 5.2.2
    "bolt_shear": 0.60,         # cl. 8.5.2
    "bolt_bearing": 0.70,       # cl. 8.5.3
    "bolt_tension": 0.70,       # cl. 8.5.4
    "slip": 1.00,               # cl. 8.6.3.2 (0.85 if the slot is parallel to force)
    "weld": 0.70,               # cl. 9.6.4.2
    "groove_weld": 0.85,        # cl. 9.5.3.3
    "rupture": 0.70,            # cl. 8.9
    "bearing_concrete": 0.60,   # cl. 8.11
    "bearing_surface": 0.70,    # cl. 8.12
    "web_yielding": 0.95,       # cl. 10.3
    "web_crippling": 0.70,      # cl. 10.4
    "flange_bending": 0.85,     # cl. 10.2
    "sidesway": 0.80,           # cl. 10.5
}

#: Table 2.3 - maximum slenderness ratio KL/r.
LAMBDA_MAX = {"compression": 180, "bracing": 200, "tension": 300}


# ---------------------------------------------------------------- conversions
def t_to_kN(x: float) -> float:
    """tons -> kN."""
    return x * 9.80665


def tcm2_to_MPa(x: float) -> float:
    """t/cm2 -> MPa."""
    return x * 98.0665


def tcm_to_kNm(x: float) -> float:
    """t.cm -> kN.m."""
    return x * 0.0980665


def _fmt(x) -> str:
    return f"{x:0.4g}" if isinstance(x, float) else str(x)


@dataclass
class Check:
    """The result of one limit-state check.

    Carries the demand, the design capacity, the clause it came from, the
    governing failure mode and every intermediate value - which is what a
    calculation sheet needs in order to be reviewable.

    >>> c = Check("Flexure", "5.1.3", demand=2400, capacity=2666)
    >>> round(c.ratio, 3)
    0.9
    >>> c.ok
    True
    """

    name: str
    clause: str
    demand: float = 0.0
    capacity: float = 0.0
    values: dict = field(default_factory=dict)
    governing: str = ""
    note: str = ""

    @property
    def ratio(self) -> float:
        """Utilisation, demand / capacity. ``inf`` if the capacity is zero."""
        return float("inf") if self.capacity == 0 else abs(self.demand) / self.capacity

    @property
    def ok(self) -> bool:
        """True when the utilisation is at most 1.0."""
        return self.ratio <= 1.0

    @property
    def status(self) -> str:
        return "OK" if self.ok else "FAIL"

    def as_dict(self) -> dict:
        """Flat dict, for writing a results table to CSV / DataFrame."""
        return {
            "check": self.name, "clause": self.clause,
            "demand": self.demand, "capacity": self.capacity,
            "ratio": self.ratio, "status": self.status,
            "governing": self.governing, **self.values,
        }

    def __str__(self) -> str:
        head = (f"{self.name:<28s} [{self.clause}]  "
                f"Ru={_fmt(self.demand)}  Rd={_fmt(self.capacity)}  "
                f"ratio={self.ratio:0.3f}  {self.status}")
        lines = [head]
        if self.governing:
            lines.append(f"    governing: {self.governing}")
        if self.values:
            lines.append("    " + "  ".join(f"{k}={_fmt(v)}"
                                            for k, v in self.values.items()))
        if self.note:
            lines.append(f"    note: {self.note}")
        return "\n".join(lines)
