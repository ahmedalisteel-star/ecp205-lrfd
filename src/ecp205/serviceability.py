"""Serviceability - ECP 205-LRFD chapter 14.

Every check here uses **unfactored (service)** loads.
"""

from __future__ import annotations

from .core import Check, E

__all__ = ["DEFLECTION_LIMITS", "DRIFT_LIMITS", "deflection_check",
           "drift_check", "udl_deflection", "required_Ix"]

#: Table 14.1 - (delta2_max, delta_max) as span/x. L is the span, or **twice**
#: the projecting length of a cantilever.
DEFLECTION_LIMITS = {
    "brittle_finish": (300, 250),           # beams and trusses carrying plaster
    "floor_supporting_columns": (500, 400),
    "other_beams": (200, 160),
    "cantilever": (180, 140),
    "crane_track": (800, 650),
    "appearance": (None, 250),              # where delta_max impairs appearance
}

#: Table 14.2 - horizontal deflection (drift) limits, height/x.
DRIFT_LIMITS = {
    "single_storey": 300,      # top of columns, single-storey, not portal frames
    "storey": 300,             # each storey of a multi-storey building
    "building_total": 500,     # top of a multi-storey building
    "portal_no_crane": 140,
    "portal_with_crane": 140,  # and per the crane manufacturer
}


def deflection_check(delta: float, span: float, member: str = "other_beams",
                     which: str = "delta_max") -> Check:
    """Vertical deflection against Table 14.1 (cl. 14.2.1.2).

    Args:
        delta: the computed **service** deflection, cm.
        span: span in cm. For a cantilever pass **twice** the projecting length.
        member: a key of :data:`DEFLECTION_LIMITS`.
        which: ``"delta2_max"`` - variable load plus time-dependent deformation
            under the permanent load; or ``"delta_max"`` - the total sagging in
            the final state, after any camber (eq. 14.1).

    >>> c = deflection_check(1.8, 600, "other_beams")
    >>> round(c.capacity, 3), c.ok
    (3.75, True)
    """
    if member not in DEFLECTION_LIMITS:
        raise ValueError(f"member must be one of {list(DEFLECTION_LIMITS)}")
    if which not in ("delta2_max", "delta_max"):
        raise ValueError("which must be 'delta2_max' or 'delta_max'")
    d2, dmax = DEFLECTION_LIMITS[member]
    div = d2 if which == "delta2_max" else dmax
    if div is None:
        raise ValueError(f"Table 14.1 gives no {which} limit for {member!r}")
    return Check(f"Deflection ({member})", "14.2.1.2 Table 14.1", delta,
                 span / div, {"span": span, "limit": f"L/{div}"},
                 f"{which}, service loads")


def drift_check(delta: float, height: float, case: str = "storey") -> Check:
    """Horizontal deflection against Table 14.2 (cl. 14.2.3), unfactored loads.

    >>> round(drift_check(1.2, 400, "storey").capacity, 4)
    1.3333
    """
    if case not in DRIFT_LIMITS:
        raise ValueError(f"case must be one of {list(DRIFT_LIMITS)}")
    div = DRIFT_LIMITS[case]
    return Check(f"Drift ({case})", "14.2.3 Table 14.2", delta, height / div,
                 {"height": height, "limit": f"H/{div}"}, "unfactored loads")


def udl_deflection(w: float, L: float, Ix: float) -> float:
    """Mid-span deflection of a simply supported beam under a UDL, cm::

        delta = 5 w L^4 / (384 E Ix)

    Args:
        w: **service** load per unit length, t/cm.
        L: span, cm.
        Ix: second moment of area, cm4.

    >>> round(udl_deflection(0.024, 800, 23130), 3)
    2.635
    """
    return 5 * w * L ** 4 / (384 * E * Ix)


def required_Ix(w: float, L: float, span_over: float) -> float:
    """Second moment of area needed to hold a UDL deflection to ``L/span_over``.

    >>> round(required_Ix(0.024, 800, 250), 0)
    19048.0
    """
    return 5 * w * L ** 4 * span_over / (384 * E * L)
