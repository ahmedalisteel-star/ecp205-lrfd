"""Connections - bolts (ch. 8), welds (ch. 9) and local effects of concentrated
forces (ch. 10)."""

from __future__ import annotations

import math

from .core import PHI, Check, E
from .materials import BOLT_AREAS, BOLT_GRADES, BOLT_LOW_SHEAR, FRICTION, SLIP_FS

__all__ = [
    # bolts
    "bolt_shear", "bolt_bearing", "bolt_tension", "bolt_combined",
    "bolt_pretension", "installation_torque", "slip_resistance",
    "prying_force", "end_plate_thickness", "end_plate_exact",
    "long_joint_factor", "packing_factor", "base_plate_bearing",
    # welds
    "min_fillet_size", "fillet_weld", "fillet_weld_combined", "groove_weld",
    # concentrated forces
    "flange_local_bending", "web_local_yielding", "web_crippling",
]


# ============================================================ ch. 8 bolts
def _areas(d_mm: int) -> tuple[float, float]:
    if d_mm not in BOLT_AREAS:
        raise ValueError(f"M{d_mm} is not tabulated; available: "
                         f"{sorted(BOLT_AREAS)}")
    return BOLT_AREAS[d_mm]


def bolt_shear(d_mm: int, grade: str, n_planes: int = 1,
               threads_in_plane: bool = True, Ru: float = 0.0) -> Check:
    """Design shear strength of one bolt, cl. 8.5.2, phi_v = 0.60.

    ``Rnv = phi_v (0.6 Fub) As n`` for grades 4.6, 5.6 and 8.8 (eq. 8.2);
    ``0.5 Fub`` for grades 4.8, 5.8, 6.8 and 10.9 (eq. 8.3).

    Args:
        d_mm: nominal diameter, mm.
        grade: bolt grade, Table 8.1.
        n_planes: number of shear planes.
        threads_in_plane: True uses the tensile stress area As; False uses the
            gross shank area A.

    >>> round(bolt_shear(20, "8.8").capacity, 3)
    7.056
    """
    A, As = _areas(d_mm)
    if grade not in BOLT_GRADES:
        raise ValueError(f"unknown bolt grade {grade!r}")
    Fub = BOLT_GRADES[grade][1]
    area = As if threads_in_plane else A
    k = 0.5 if grade in BOLT_LOW_SHEAR else 0.6
    cap = PHI["bolt_shear"] * k * Fub * area * n_planes
    return Check(f"Bolt shear M{d_mm}-{grade}", "8.5.2", Ru, cap,
                 {"Fub": Fub, "area": area, "planes": n_planes, "k": k},
                 f"eq. {'8.3' if k == 0.5 else '8.2'} ({k} Fub)")


def bolt_bearing(d_mm: int, sum_t: float, Fu: float, e1: float, Ru: float = 0.0,
                 long_slot_perp: bool = False) -> Check:
    """Design bearing strength of one bolt, cl. 8.5.3 eq. 8.4, phi_br = 0.70.

    ``Rbr = phi_br d (sum t) (alpha Fu)``, with alpha from Table 8.2 and capped
    by ``alpha = 0.8 e1/d <= 2.40``.

    Args:
        d_mm: nominal diameter, mm.
        sum_t: smallest sum of plate thicknesses in one direction, cm.
        Fu: ultimate strength of the connected **plates**, t/cm2.
        e1: end distance along the line of force, cm.
        long_slot_perp: long slots perpendicular to the force reduce alpha 20 %.

    >>> round(bolt_bearing(20, 1.2, 3.7, e1=6.0).values["alpha"], 2)
    2.4
    """
    d = d_mm / 10
    if e1 < 1.5 * d:
        raise ValueError(
            f"end distance e1 = {e1:g} cm is below the cl. 8.4.2 minimum of "
            f"1.5d = {1.5 * d:g} cm")
    if e1 >= 3 * d:
        a = 2.4
    elif e1 >= 2.5 * d:
        a = 2.0
    elif e1 >= 2.0 * d:
        a = 1.6
    else:
        a = 1.2
    a = min(a, 0.8 * e1 / d)
    if long_slot_perp:
        a *= 0.8
    cap = PHI["bolt_bearing"] * d * sum_t * a * Fu
    return Check(f"Bolt bearing M{d_mm}", "8.5.3", Ru, cap,
                 {"alpha": a, "d": d, "sum_t": sum_t, "Fu": Fu, "e1": e1},
                 "eq. 8.4, alpha from Table 8.2")


def bolt_tension(d_mm: int, grade: str, Ru: float = 0.0,
                 pretensioned: bool = False, prying: float = 0.0) -> Check:
    """Design tension strength of one bolt, phi_t = 0.70.

    ``0.66 Fub As`` for snug-tight grades 4.6-6.8 (cl. 8.5.4 eq. 8.5);
    ``0.80 Fub As`` for pretensioned grades 8.8 and 10.9 (cl. 8.6.4 eq. 8.10).

    Args:
        prying: prying force P per bolt from :func:`prying_force`, added to
            the demand (cl. 8.6.4 eq. 8.9).

    Raises:
        ValueError: grades 8.8 and 10.9 in tension must be pretensioned
            (cl. 8.5.4b).
    """
    _, As = _areas(d_mm)
    if grade not in BOLT_GRADES:
        raise ValueError(f"unknown bolt grade {grade!r}")
    if not pretensioned and grade in ("8.8", "10.9"):
        raise ValueError(
            f"cl. 8.5.4b: grade {grade} bolts in tension must be pretensioned "
            f"- pass pretensioned=True and detail to cl. 8.6.8")
    Fub = BOLT_GRADES[grade][1]
    k, eq, clause = ((0.80, "eq. 8.10", "8.6.4") if pretensioned
                     else (0.66, "eq. 8.5", "8.5.4"))
    cap = PHI["bolt_tension"] * k * Fub * As
    return Check(f"Bolt tension M{d_mm}-{grade}", clause, Ru + prying, cap,
                 {"Fub": Fub, "As": As, "k": k, "prying P": prying}, eq)


def bolt_combined(Rut: float, phiRnt: float, Ruv: float, phiRnv: float) -> Check:
    """Combined shear and tension, cl. 8.5.5 eq. 8.6 - circular interaction::

        (Rut/phiRnt)^2 + (Ruv/phiRnv)^2 <= 1.0

    >>> round(bolt_combined(3, 5, 4, 7).ratio, 4)
    0.8286
    """
    ratio = math.sqrt((Rut / phiRnt) ** 2 + (Ruv / phiRnv) ** 2)
    return Check("Bolt shear+tension", "8.5.5", ratio, 1.0,
                 {"Rut/phiRnt": Rut / phiRnt, "Ruv/phiRnv": Ruv / phiRnv},
                 "eq. 8.6")


def bolt_pretension(d_mm: int, grade: str = "10.9") -> float:
    """Pretension force ``T = 0.7 Fyb As``, tons - cl. 8.6.2.1 eq. 8.7.

    Reproduces Table 8.3: M20 grade 10.9 gives 15.44 t against 15.43 tabulated.

    >>> round(bolt_pretension(20, "10.9"), 2)
    15.44
    """
    Fyb = BOLT_GRADES[grade][0]
    return 0.7 * Fyb * _areas(d_mm)[1]


def installation_torque(d_mm: int, grade: str = "10.9") -> float:
    """Installation torque ``Ma = 0.2 d T`` - cl. 8.6.8 eq. 8.17.

    Returned in **t.cm**, consistent with the rest of the package. Table 8.3
    prints the same quantity in kg.m, which is 10x this value: M20 grade 10.9
    gives 6.17 t.cm here against 62 kg.m tabulated.

    >>> round(installation_torque(20, "10.9"), 2)
    6.17
    """
    return 0.2 * (d_mm / 10) * bolt_pretension(d_mm, grade)


def slip_resistance(d_mm: int, grade: str = "10.9", surface: str = "B",
                    case: str = "I", n_surfaces: int = 1,
                    slot_parallel: bool = False, Qb: float = 0.0) -> Check:
    """Slip resistance per bolt, cl. 8.6.3.2 eq. 8.8::

        Ps = phi mu T / F.S.

    This is a **serviceability** check at unfactored load, with a factor of
    safety rather than a phi factor. The same bolts must also satisfy cl. 8.5.2
    shear and cl. 8.5.3 bearing at factored load (cl. 8.6.3.3).

    Args:
        surface: friction class A (0.50), B (0.40) or C (0.30), cl. 8.6.2.2.
        case: ``"I"`` primary loads, ``"II"`` primary plus wind/earthquake;
            ``"crane_I"`` / ``"crane_II"`` for cranes and crane girders.
        n_surfaces: number of friction surfaces.
        slot_parallel: phi drops from 1.0 to 0.85 if the slot is parallel to
            the line of force.
        Qb: the **service** shear per bolt, tons.

    Reproduces Table 8.3: M20 grade 10.9, class B, case I gives 4.94 t against
    4.93 tabulated.

    >>> round(slip_resistance(20, "10.9", "B", "I").capacity, 2)
    4.94
    """
    if surface not in FRICTION:
        raise ValueError(f"surface class must be one of {list(FRICTION)}")
    if case not in SLIP_FS:
        raise ValueError(f"loading case must be one of {list(SLIP_FS)}")
    if grade not in ("8.8", "10.9"):
        raise ValueError("cl. 8.6.1: only grades 8.8 and 10.9 may be used in "
                         "slip-critical connections")
    T = bolt_pretension(d_mm, grade)
    mu, fs = FRICTION[surface], SLIP_FS[case]
    phi = 0.85 if slot_parallel else 1.0
    cap = phi * mu * T / fs * n_surfaces
    return Check(f"Slip resistance M{d_mm}-{grade}", "8.6.3.2", Qb, cap,
                 {"T": T, "mu": mu, "F.S.": fs, "phi": phi,
                  "surfaces": n_surfaces},
                 "eq. 8.8 - SERVICEABILITY limit state, service loads")


def prying_force(b: float, a: float, tp: float, T_ext: float, w: float = 0.0,
                 As: float = 0.0, refined: bool = False) -> float:
    """Prying force P per bolt, cl. 8.8.3. Lengths cm, forces tons.

    Simple, eq. 8.26::

        P = [ 5b/(8a) - tp^3/100 ] T_ext

    Refined, eq. 8.27::

        P = [ (1/2 - k) / ((3a/4b)(a/4b + 1) + k) ] T_ext ,
        k = w tp^4 / (30 a b^2 As)

    Args:
        b: inner bolt dimension with respect to the tee-stub web.
        a: outer overhang of the tee-stub flange.
        tp: flange or end plate thickness.
        T_ext: factored external tension per bolt.
        w: plate breadth tributary to one bolt - eq. 8.27 only.
        As: bolt tensile stress area - eq. 8.27 only.

    Returns:
        P, clipped at zero - a negative result means no prying occurs.

    >>> round(prying_force(5, 4, 2.0, 12), 3)
    8.415
    """
    if a <= 0 or b <= 0:
        raise ValueError("a and b must be positive")
    if refined:
        if not (w and As):
            raise ValueError("eq. 8.27 needs w (plate breadth) and As")
        k = w * tp ** 4 / (30 * a * b ** 2 * As)
        P = (0.5 - k) / ((3 * a / (4 * b)) * (a / (4 * b) + 1) + k) * T_ext
    else:
        P = (5 * b / (8 * a) - tp ** 3 / 100) * T_ext
    return max(P, 0.0)


def end_plate_thickness(Mu: float, b: float, s: float, tb: float, db: float,
                        w: float, Fy: float, bolts_per_row: int = 2) -> float:
    """Approximate end plate / tee-stub thickness, cl. 8.8.4 eq. 8.28 / 8.29, cm.

    ``tp = k sqrt( Mu (2b + 2s + tb) / (db w Fy) )``, k = 0.50 for two bolts per
    row (four around the tension flange) or 0.353 for four per row (eight).

    Follow with :func:`end_plate_exact` for the eq. 8.30 / 8.31 thickness.
    """
    if bolts_per_row not in (2, 4):
        raise ValueError("bolts_per_row must be 2 (eq. 8.28) or 4 (eq. 8.29)")
    k = 0.50 if bolts_per_row == 2 else 0.353
    return k * math.sqrt(Mu * (2 * b + 2 * s + tb) / (db * w * Fy))


def end_plate_exact(P: float, a: float, b: float, T_extM: float, w: float,
                    Fy: float) -> dict:
    """Exact end plate thickness, cl. 8.8.4 eq. 8.30 / 8.31::

        M1 = P a ;  M2 = P a - T_ext,b,M b ;  tp = sqrt( 4 max(M1,M2) / (w Fy) )
    """
    M1 = P * a
    M2 = P * a - T_extM * b
    return {"M1": M1, "M2": M2, "tp": math.sqrt(4 * max(M1, M2) / (w * Fy)),
            "clause": "8.8.4 eq. 8.30/8.31"}


def long_joint_factor(Lj: float, d_mm: float) -> float:
    """Long-joint reduction on shear and bearing, cl. 8.10.2 eq. 8.37::

        beta_L = 1 - (Lj - 15d)/(200 d) ,  0.75 <= beta_L <= 1.0

    Applies when Lj > 15d; Lj must not exceed 65d. Not applicable where the
    force transfer is uniform along the joint.

    >>> long_joint_factor(25, 20)          # Lj <= 15d, no reduction
    1.0
    >>> round(long_joint_factor(60, 20), 4)
    0.925
    """
    d = d_mm / 10
    if Lj > 65 * d:
        raise ValueError(f"cl. 8.10.2: Lj must not exceed 65d = {65 * d:g} cm")
    if Lj <= 15 * d:
        return 1.0
    return max(min(1 - (Lj - 15 * d) / (200 * d), 1.0), 0.75)


def packing_factor(tp: float, d_mm: float) -> float:
    """Packing reduction on shear, cl. 8.10.4 eq. 8.39 - ``9d/(8d + 3tp) <= 1``,
    applicable when the total packing thickness exceeds d/3."""
    d = d_mm / 10
    return min(9 * d / (8 * d + 3 * tp), 1.0) if tp > d / 3 else 1.0


def base_plate_bearing(A1: float, fcu: float, A2: float = 0.0,
                       Pu: float = 0.0) -> Check:
    """Bearing on concrete under a column base, cl. 8.11, phi_c = 0.60.

    ``Pp = 0.85 f'c A1`` on the full area (eq. 8.41), or
    ``Pp = 0.85 f'c A1 sqrt(A2/A1) <= 0.85 f'c (2 A1)`` on a partial area
    (eq. 8.42).

    Args:
        A1: area of the steel plate bearing on the support, cm2.
        fcu: specified 28-day compressive strength, t/cm2.
        A2: geometrically similar concentric area of the support, cm2;
            0 or equal to A1 for bearing on the full area.
    """
    if A2 and A2 > A1:
        Pp = min(0.85 * fcu * A1 * math.sqrt(A2 / A1), 0.85 * fcu * 2 * A1)
        gov = "partial area, eq. 8.42"
    else:
        Pp = 0.85 * fcu * A1
        gov = "full area, eq. 8.41"
    return Check("Base plate bearing", "8.11", Pu,
                 PHI["bearing_concrete"] * Pp,
                 {"A1": A1, "A2": A2, "fcu": fcu, "Pp": Pp}, gov)


# ============================================================ ch. 9 welds
def min_fillet_size(t_thicker_mm: float) -> float:
    """Minimum fillet leg in mm from the thicker part joined, Table 9.6.

    >>> min_fillet_size(10), min_fillet_size(25)
    (5.0, 8.0)
    """
    if t_thicker_mm <= 6:
        return 3.0
    if t_thicker_mm <= 12:
        return 5.0
    if t_thicker_mm <= 18:
        return 6.0
    return 8.0


def fillet_weld(s: float, Fu: float, length: float = 1.0,
                Ru: float = 0.0) -> Check:
    """Design strength of a fillet weld, cl. 9.6.4.2 eq. 9.7, phi = 0.70::

        Rnw = phi s (0.4 Fu)   per unit length

    Args:
        s: **leg** size, cm.
        Fu: ultimate strength of the **base** metal, t/cm2 (the smaller Fu for
            a hybrid joint).
        length: effective length, cm - the overall length less 2s for the end
            craters (cl. 9.6.4.2).

    >>> round(fillet_weld(0.6, 3.7, length=1.0).capacity, 4)
    0.6216
    """
    q = PHI["weld"] * s * 0.4 * Fu
    return Check("Fillet weld", "9.6.4.2", Ru, q * length,
                 {"s": s, "Fu": Fu, "length": length, "q per cm": q}, "eq. 9.7")


def fillet_weld_combined(R_perp: float, R_par: float, R_shear_perp: float,
                         s: float, Fu: float) -> Check:
    """Combined normal and shear on a fillet weld, cl. 9.6.4.2 eq. 9.8::

        Rw,eff = sqrt( R_perp^2 + 3 (R_par^2 + R_shear_perp^2) )
               <= 0.77 s (0.4 Fu)

    The 10 % increase permitted for this combination is already in the 0.77.
    All the R terms are forces per unit length of weld.
    """
    eff = math.sqrt(R_perp ** 2 + 3 * (R_par ** 2 + R_shear_perp ** 2))
    return Check("Fillet weld combined", "9.6.4.2", eff, 0.77 * s * 0.4 * Fu,
                 {"R_perp": R_perp, "R_par": R_par,
                  "R_shear_perp": R_shear_perp, "s": s, "Fu": Fu}, "eq. 9.8")


def groove_weld(t: float, Fy: float, shear: bool = False, Ru: float = 0.0,
                length: float = 1.0) -> Check:
    """Groove (butt) weld, cl. 9.5.3.3, phi = 0.85::

        tension / compression   Rnw = 0.85 t Fy        (eq. 9.1)
        shear                   Rnw = 0.85 t (0.6 Fy)  (eq. 9.2)

    Args:
        t: base metal thickness for a complete penetration weld, or the actual
            penetrated depth tg for a partial penetration weld (eq. 9.3 / 9.4).
    """
    cap = 0.85 * t * (0.6 * Fy if shear else Fy) * length
    return Check("Groove weld", "9.5.3.3", Ru, cap,
                 {"t": t, "Fy": Fy, "shear": shear},
                 "eq. 9.2" if shear else "eq. 9.1")


# ============================================================ ch. 10
def flange_local_bending(tf: float, Fyf: float, Ru: float = 0.0,
                         near_end: bool = False) -> Check:
    """Flange local bending under a tensile concentrated force, cl. 10.2
    eq. 10.1 - ``Rn = 6.25 tf^2 Fyf``, phi = 0.85.

    Args:
        near_end: True halves Rn when the force acts within 10 tf of the member
            end.

    Need not be checked when the force is applied over a length below 0.15 bf.
    """
    Rn = 6.25 * tf ** 2 * Fyf
    if near_end:
        Rn *= 0.5
    return Check("Flange local bending", "10.2", Ru,
                 PHI["flange_bending"] * Rn,
                 {"tf": tf, "Fyf": Fyf, "near_end": near_end}, "eq. 10.1")


def web_local_yielding(k: float, N: float, Fyw: float, tw: float,
                       Ru: float = 0.0, near_end: bool = False) -> Check:
    """Web local yielding, cl. 10.3, phi = 0.95.

    ``Rn = (5k + N) Fyw tw`` (eq. 10.2), or ``(2.5k + N) Fyw tw`` within d of
    the member end (eq. 10.3).

    Args:
        k: flange outer face to the web toe of the fillet, cm (tf + r).
        N: bearing length, cm - not less than k for an end reaction.
    """
    Rn = ((2.5 if near_end else 5.0) * k + N) * Fyw * tw
    return Check("Web local yielding", "10.3", Ru,
                 PHI["web_yielding"] * Rn,
                 {"k": k, "N": N, "tw": tw, "Fyw": Fyw},
                 "eq. 10.3" if near_end else "eq. 10.2")


def web_crippling(tw: float, tf: float, d: float, N: float, Fyw: float,
                  Ru: float = 0.0, near_end: bool = False) -> Check:
    """Web crippling under a compressive concentrated force, cl. 10.4, phi = 0.70.

    eq. 10.4 at d/2 or more from the end; eq. 10.5a / 10.5b nearer than d/2,
    split on N/d = 0.2.

    A short bearing near a beam end is the usual failure - the end-zone
    coefficient is half the interior one.
    """
    base = math.sqrt(E * Fyw * tf / tw)
    if not near_end:
        Rn = 0.3627 * tw ** 2 * (1 + 3 * (N / d) * (tw / tf) ** 1.5) * base
        eq = "eq. 10.4"
    elif N / d <= 0.2:
        Rn = 0.1813 * tw ** 2 * (1 + 3 * (N / d) * (tw / tf) ** 1.5) * base
        eq = "eq. 10.5a"
    else:
        Rn = 0.1813 * tw ** 2 * (1 + (4 * N / d - 0.2) * (tw / tf) ** 1.5) * base
        eq = "eq. 10.5b"
    return Check("Web crippling", "10.4", Ru, PHI["web_crippling"] * Rn,
                 {"tw": tw, "tf": tf, "N": N, "N/d": N / d, "Fyw": Fyw}, eq)
