"""Steel and bolt materials - ECP 205-LRFD chapter 1 and clause 8.1."""

from __future__ import annotations

__all__ = [
    "STEEL_GRADES", "BOLT_GRADES", "BOLT_LOW_SHEAR", "BOLT_AREAS",
    "FRICTION", "SLIP_FS", "steel", "bolt",
]

#: Table 1.1 (cl. 1.3.3), Egyptian Standard 260/2004.
#: grade -> ((Fy, Fu) for t <= 40 mm, (Fy, Fu) for 40 < t <= 100 mm), t/cm2
STEEL_GRADES = {
    "St 37": ((2.40, 3.70), (2.15, 3.40)),
    "St 44": ((2.80, 4.40), (2.55, 4.10)),
    "St 52": ((3.60, 5.20), (3.35, 4.90)),
}

#: Table 8.1 - bolt grade -> (Fyb, Fub), t/cm2
BOLT_GRADES = {
    "4.6": (2.4, 4.0), "4.8": (3.2, 4.0),
    "5.6": (3.0, 5.0), "5.8": (4.0, 5.0),
    "6.8": (4.8, 6.0),
    "8.8": (6.4, 8.0), "10.9": (9.0, 10.0),
}

#: Grades whose shear strength uses 0.5*Fub instead of 0.6*Fub (cl. 8.5.2b, eq. 8.3).
BOLT_LOW_SHEAR = frozenset({"4.8", "5.8", "6.8", "10.9"})

#: Nominal diameter (mm) -> (shank area A, tensile stress area As), cm2.
BOLT_AREAS = {
    12: (1.13, 0.84), 16: (2.01, 1.57), 20: (3.14, 2.45), 22: (3.80, 3.03),
    24: (4.52, 3.53), 27: (5.73, 4.59), 30: (7.06, 5.61), 36: (10.18, 8.17),
}

#: Friction coefficient by surface class (cl. 8.6.2.2).
#: A - blasted, unpainted or metallized; B - blasted + alkali-zinc silicate 50-80 um;
#: C - wire brushed or flame cleaned.
FRICTION = {"A": 0.50, "B": 0.40, "C": 0.30}

#: Factor of safety against slip (cl. 8.6.3.2).
#: Case I  - primary stresses (dead, live, dynamic, centrifugal).
#: Case II - case I plus wind/earthquake, braking, shock, temperature, settlement.
SLIP_FS = {"I": 1.25, "II": 1.05, "crane_I": 1.60, "crane_II": 1.35}


def steel(grade: str = "St 37", t_mm: float = 10.0) -> tuple[float, float]:
    """Yield and ultimate strength for a grade and element thickness (Table 1.1).

    Args:
        grade: ``"St 37"``, ``"St 44"`` or ``"St 52"``.
        t_mm: thickness of the element being checked, in **millimetres**.
            Fy and Fu drop above 40 mm.

    Returns:
        ``(Fy, Fu)`` in t/cm2.

    Raises:
        ValueError: unknown grade, or a thickness beyond the 100 mm the table covers.

    >>> steel("St 37", 15)
    (2.4, 3.7)
    >>> steel("St 52", 60)
    (3.35, 4.9)
    """
    if grade not in STEEL_GRADES:
        raise ValueError(f"unknown grade {grade!r}; Table 1.1 has {list(STEEL_GRADES)}")
    if t_mm > 100:
        raise ValueError(
            "Table 1.1 stops at t = 100 mm - specify Fy and Fu explicitly, "
            "and see cl. 1.3.1.3 for heavy sections over 50 mm in tension")
    thin, thick = STEEL_GRADES[grade]
    return thin if t_mm <= 40 else thick


def bolt(grade: str) -> tuple[float, float]:
    """``(Fyb, Fub)`` in t/cm2 for a bolt grade (Table 8.1).

    >>> bolt("10.9")
    (9.0, 10.0)
    """
    if grade not in BOLT_GRADES:
        raise ValueError(f"unknown bolt grade {grade!r}; Table 8.1 has "
                         f"{list(BOLT_GRADES)}. Grades below 4.6 or above 10.9 "
                         f"are not permitted (cl. 8.1.1).")
    return BOLT_GRADES[grade]
