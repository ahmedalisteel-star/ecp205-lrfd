"""Member strength - tension (ch. 3), compression (ch. 4), flexure and shear
(ch. 5), and the beam-column interaction (ch. 2.2.2 and 7)."""

from __future__ import annotations

import math
from typing import Optional

from .classification import classify, flange_limits, web_limits
from .core import LAMBDA_MAX, PHI, Check, E
from .sections import ISection

__all__ = [
    # tension
    "tension", "shear_lag_U", "block_shear",
    # compression
    "compression", "effective_width_unstiffened", "effective_width_stiffened",
    "Q_factor",
    # flexure and shear
    "Cb_ends", "Cb_quarter", "flexure_I", "shear_web",
    # beam-column
    "Pe", "Cm_endmoments", "B1", "B2_drift", "B2_euler", "beam_column",
]


# ============================================================ ch. 3 tension
def shear_lag_U(x_bar: float, L_conn: float) -> float:
    """Shear lag reduction coefficient, cl. 2.1.3 eq. 2.2: ``U = 1 - x/l <= 0.9``.

    Args:
        x_bar: connection eccentricity - distance from the connection plane to
            the centroid of the part resisting the force, cm.
        L_conn: length of the connection along the line of force, cm.

    >>> round(shear_lag_U(2.5, 20.0), 3)
    0.875
    >>> shear_lag_U(0.1, 100.0)
    0.9
    """
    if L_conn <= 0:
        raise ValueError("connection length must be positive")
    return min(1 - x_bar / L_conn, 0.9)


def tension(Ag: float, Ae: float, Fy: float, Fu: float, Pu: float = 0.0) -> Check:
    """Design tensile strength, cl. 3.1.1 - the lower of gross yielding
    (eq. 3.1, phi = 0.85) and net fracture (eq. 3.2, phi = 0.70).

    Args:
        Ag: gross area, cm2.
        Ae: effective net area from cl. 2.1.3 (shear lag), cm2.
        Fy, Fu: t/cm2.
        Pu: required tensile strength, tons.

    >>> c = tension(84.5, 75.0, 2.4, 3.7, Pu=120)
    >>> round(c.capacity, 1), c.ok
    (172.4, True)
    """
    py = PHI["tension_yield"] * Fy * Ag
    pf = PHI["tension_fracture"] * Fu * Ae
    cap, gov = ((py, "yielding of the gross section, eq. 3.1") if py <= pf
                else (pf, "fracture of the net section, eq. 3.2"))
    return Check("Tension", "3.1.1", Pu, cap,
                 {"Ag": Ag, "Ae": Ae, "Fy": Fy, "Fu": Fu,
                  "phi*Fy*Ag": py, "phi*Fu*Ae": pf}, gov)


def block_shear(Agv: float, Anv: float, Agt: float, Ant: float,
                Fy: float, Fu: float, Ru: float = 0.0) -> Check:
    """Block shear rupture, cl. 8.9.3 eq. 8.35 / 8.36, phi = 0.70.

    Where rupture on the net section governs one segment, yielding on the gross
    section is used on the perpendicular segment.

    Args:
        Agv, Anv: gross and net area subject to shear, cm2.
        Agt, Ant: gross and net area subject to tension, cm2.
    """
    phi = PHI["rupture"]
    if Fu * Ant >= 0.6 * Fu * Anv:
        cap = phi * (0.6 * Fy * Agv + Fu * Ant)
        gov = "shear yield + tension rupture, eq. 8.35"
    else:
        cap = phi * (0.6 * Fu * Anv + Fy * Agt)
        gov = "shear rupture + tension yield, eq. 8.36"
    return Check("Block shear", "8.9.3", Ru, cap,
                 {"Agv": Agv, "Anv": Anv, "Agt": Agt, "Ant": Ant}, gov)


# ============================================================ ch. 4 compression
def compression(Ag: float, r: float, KL: float, Fy: float, Q: float = 1.0,
                Pu: float = 0.0, member: str = "compression") -> Check:
    """Design compressive strength for flexural buckling, cl. 4.2, phi_c = 0.80.

    ``lambda_c <= 1.1``: ``Fcr = Fy Q (1 - 0.384 Q lambda_c^2)`` (eq. 4.2 / 4.6)
    ``lambda_c >  1.1``: ``Fcr = 0.648 Fy / lambda_c^2``          (eq. 4.3 / 4.7)

    Args:
        Ag: gross area, cm2 - always the actual gross section, even when Q < 1.
        r: radius of gyration about the buckling axis, cm.
        KL: effective buckling length, cm (K from cl. 2.2.1.3).
        Fy: t/cm2.
        Q: Ae/Ag reduction for slender elements (cl. 4.2.2); 1.0 otherwise.
        Pu: required compressive strength, tons.
        member: which Table 2.3 slenderness limit to report against -
            ``"compression"`` (180), ``"bracing"`` (200).

    >>> c = compression(84.5, 3.95, 400, 2.4, Pu=60)
    >>> round(c.capacity, 1)
    88.3
    """
    if r <= 0 or KL <= 0:
        raise ValueError("KL and r must be positive")
    lam = KL / r
    Fe = math.pi ** 2 * E / lam ** 2                       # eq. 4.5
    lc = math.sqrt(Fy / Fe)                                # eq. 4.4
    if lc * math.sqrt(Q) <= 1.1:
        Fcr = Fy * Q * (1 - 0.384 * Q * lc ** 2)
        gov = f"inelastic buckling, eq. {'4.2' if Q == 1 else '4.6'}"
    else:
        Fcr = 0.648 * Fy / lc ** 2
        gov = f"elastic buckling, eq. {'4.3' if Q == 1 else '4.7'}"
    chk = Check("Compression", "4.2", Pu, PHI["compression"] * Ag * Fcr,
                {"KL/r": lam, "Fe": Fe, "lambda_c": lc, "Q": Q,
                 "Fcr": Fcr, "Ag": Ag}, gov)
    limit = LAMBDA_MAX[member]
    if lam > limit:
        chk.note = (f"KL/r = {lam:0.0f} exceeds the Table 2.3 limit of {limit} "
                    f"for a {member} member")
    return chk


def effective_width_unstiffened(b: float, t: float, Fy: float) -> float:
    """cl. 4.2.2.1 eq. 4.8 - outstand flanges, angle legs, T stems."""
    ratio = b / t
    be = 0.78 * t * math.sqrt(E / Fy) * (1 - (0.13 / ratio) * math.sqrt(E / Fy))
    return min(max(be, 0.0), b)


def effective_width_stiffened(b: float, t: float, Fy: float) -> float:
    """cl. 4.2.2.2 eq. 4.9 - webs of I sections, flanges of box sections."""
    ratio = b / t
    be = 1.92 * t * math.sqrt(E / Fy) * (1 - (0.385 / ratio) * math.sqrt(E / Fy))
    return min(max(be, 0.0), b)


def Q_factor(sec: ISection, Fy: float) -> float:
    """``Q = Ae/Ag`` for an I section in uniform compression (cl. 4.2.2).

    Returns 1.0 when no element exceeds its lambda_r.
    """
    _, lrf = flange_limits(Fy, sec.rolled)
    _, lrw = web_limits(Fy, alpha=1.0, psi=1.0)
    be_f = (sec.c if sec.c / sec.tf <= lrf
            else effective_width_unstiffened(sec.c, sec.tf, Fy))
    be_w = (sec.dw if sec.dw / sec.tw <= lrw
            else effective_width_stiffened(sec.dw, sec.tw, Fy))
    Ae = sec.A - 4 * (sec.c - be_f) * sec.tf - (sec.dw - be_w) * sec.tw
    return min(Ae / sec.A, 1.0)


# ============================================================ ch. 5 flexure
def Cb_ends(M1: float, M2: float) -> float:
    """Bending coefficient from the end moments, cl. 5.1.1 eq. 5.1::

        Cb = 1.75 + 1.05 (M1/M2) + 0.3 (M1/M2)^2 <= 2.3

    M1/M2 is **positive for reverse curvature**. Valid only for a straight-line
    moment diagram within the unbraced length; if the moment anywhere inside the
    segment exceeds both end values, use :func:`Cb_quarter` or Cb = 1.0.

    For a cantilever or overhang with a free unbraced end, Cb = 1.0 (cl. 5.1.1).

    >>> Cb_ends(0, 100)
    1.75
    >>> round(Cb_ends(-50, 100), 4)
    1.3
    """
    ratio = 0.0 if M2 == 0 else M1 / M2
    return min(1.75 + 1.05 * ratio + 0.3 * ratio ** 2, 2.3)


def Cb_quarter(Mmax: float, Ma: float, Mb: float, Mc: float) -> float:
    """Bending coefficient from the quarter-point moments, cl. 5.1.1 eq. 5.2::

        Cb = 12.5 Mmax / (2.5 Mmax + 3 Ma + 4 Mb + 3 Mc)

    Ma, Mb, Mc at the quarter, mid and three-quarter points. Absolute values.

    >>> round(Cb_quarter(100, 75, 100, 75), 3)
    1.136
    """
    den = 2.5 * abs(Mmax) + 3 * abs(Ma) + 4 * abs(Mb) + 3 * abs(Mc)
    return 12.5 * abs(Mmax) / den if den else 1.0


def flexure_I(sec: ISection, Fy: float, Lb: float, Cb: float = 1.0,
              Mu: float = 0.0, Fyw: Optional[float] = None) -> Check:
    """Design flexural strength of an I section about the major axis, cl. 5.1.3.

    Covers compact sections (cl. 5.1.3.1) and non-compact sections
    (cl. 5.1.3.2) across all three unbraced-length regimes. phi_b = 0.85.

    Slender sections (class 3) return an eq. 5.19 estimate using the **gross**
    Sx and carry a note - they need effective section properties (cl. 2.3.1.3)
    and, if the web is slender, chapter 6.

    Args:
        sec: the section.
        Fy: flange yield strength, t/cm2.
        Lb: laterally unbraced length of the compression flange, cm.
        Cb: bending coefficient, :func:`Cb_ends` or :func:`Cb_quarter`.
        Mu: required flexural strength, t.cm.
        Fyw: web yield strength if hybrid; defaults to Fy.

    >>> from .sections import get
    >>> c = flexure_I(get("IPE 400"), 2.4, Lb=150)
    >>> round(c.capacity), c.governing
    (2666, 'full plastic moment, eq. 5.3 (Lb <= Lp)')
    """
    if Lb < 0:
        raise ValueError("Lb must be non-negative")
    Fyw = Fy if Fyw is None else Fyw
    cl = classify(sec, Fy, "bending")
    Sx = sec.Sx
    My = Fy * Sx
    Mp = min(Fy * sec.Zx, 1.5 * My)                        # cl. 5.1.3.1

    FL = 0.75 * Fy if sec.rolled else 0.60 * min(Fy, Fyw)
    Mr = FL * Sx                                           # eq. 5.7
    Lp = 80 * sec.ry / math.sqrt(Fy)                       # eq. 5.4
    X = (0.104 * sec.rT * sec.d / sec.Af) ** 2             # eq. 5.8
    Lr = (1380 * sec.Af / (sec.d * FL)) * math.sqrt(
        0.5 * (1 + math.sqrt(1 + (2 * X * FL) ** 2)))      # eq. 5.7

    def Mcr(Lb_: float) -> float:                          # eq. 5.12
        t1 = 1380 * sec.Af / (sec.d * Lb_)
        t2 = 20700 / (Lb_ / sec.rT) ** 2
        return min(Sx * math.sqrt(t1 ** 2 + t2 ** 2), Mp)

    vals = {"class": cl["class"], "Mp": Mp, "My": My, "Mr": Mr,
            "Lp": Lp, "Lr": Lr, "Lb": Lb, "Cb": Cb, "FL": FL}
    note = ""

    if cl["class"] == 1:
        if Lb <= Lp:
            Mn, gov = Mp, "full plastic moment, eq. 5.3 (Lb <= Lp)"
        elif Lb <= Lr:
            Mn = min((Mp - (Mp - Mr) * (Lb - Lp) / (Lr - Lp)) * Cb, Mp)
            gov = "inelastic lateral-torsional buckling, eq. 5.6"
        else:
            Mn = min(Cb * Mcr(Lb), Mp)
            gov = "elastic lateral-torsional buckling, eq. 5.11/5.12"

    elif cl["class"] == 2:
        cands = [Mp - (Mp - Mr) * (p["lambda"] - p["lp"]) / (p["lr"] - p["lp"])
                 for p in (cl["flange"], cl["web"]) if p["class"] >= 2]
        Mn_loc = min(min(cands), Mp) if cands else Mp      # eq. 5.16
        Lp_dash = Lp + (Lr - Lp) * (Mp - Mn_loc) / (Mp - Mr)   # eq. 5.17
        vals["L'p"], vals["Mn_local"] = Lp_dash, Mn_loc
        if Lb <= Lp_dash:
            Mn, gov = Mn_loc, "flange/web local buckling, eq. 5.16"
        elif Lb <= Lr:
            Mltb = min((Mp - (Mp - Mr) * (Lb - Lp) / (Lr - Lp)) * Cb, Mp)
            Mn = min(Mltb, Mn_loc)
            gov = ("inelastic lateral-torsional buckling, eq. 5.6"
                   if Mltb < Mn_loc else "flange/web local buckling, eq. 5.16")
        else:
            Mn = min(Cb * Mcr(Lb), Mp)
            gov = "elastic lateral-torsional buckling, eq. 5.11/5.12"

    else:
        Mn = FL * Sx
        gov = "slender section, eq. 5.19 with GROSS Sx - provisional"
        note = ("class 3: recompute with effective section properties "
                "(cl. 2.3.1.3), and check the web as a plate girder (ch. 6) "
                "if h/tw exceeds the eq. 6.1 limit")

    vals["Mn"] = Mn
    chk = Check("Flexure Mx", "5.1.3", Mu, PHI["flexure"] * Mn, vals, gov)
    chk.note = note
    return chk


def shear_web(sec: ISection, Fyw: float, Vu: float = 0.0,
              h: Optional[float] = None) -> Check:
    """Design shear strength of an unstiffened web, cl. 5.2.2 eq. 5.22-5.24.

    phi_v = 0.85. ``Aw = d * tw`` - the **overall** depth (cl. 5.2.1).

    Args:
        h: clear web depth for the h/tw ratio; defaults to ``sec.dw``
            (clear distance between flanges less the fillets).

    Raises:
        ValueError: h/tw > 260, where transverse stiffeners are mandatory and
            the member must be designed as a plate girder (cl. 6.3).

    >>> from .sections import get
    >>> round(shear_web(get("IPE 400"), 2.4, Vu=30).capacity, 1)
    42.1
    """
    h = sec.dw if h is None else h
    ratio = h / sec.tw
    s = math.sqrt(Fyw)
    if ratio <= 112 / s:
        Vn, gov = 0.6 * Fyw * sec.Aw, "web shear yielding, eq. 5.22"
    elif ratio <= 139 / s:
        Vn = 0.6 * Fyw * sec.Aw * (112 / s) / ratio
        gov = "inelastic web shear buckling, eq. 5.23"
    elif ratio <= 260:
        Vn = 9500 / ratio ** 2 * sec.Aw
        gov = "elastic web shear buckling, eq. 5.24"
    else:
        raise ValueError(
            f"h/tw = {ratio:0.0f} > 260: cl. 5.2.2 does not apply. Transverse "
            f"stiffeners are required - design as a plate girder to cl. 6.3.")
    return Check("Shear Vy", "5.2.2", Vu, PHI["shear"] * Vn,
                 {"h/tw": ratio, "112/sqrt(Fy)": 112 / s,
                  "139/sqrt(Fy)": 139 / s, "Aw": sec.Aw, "Vn": Vn}, gov)


# ============================================================ beam-columns
def Pe(I: float, KL: float) -> float:
    """Euler load ``pi^2 E I / (KL)^2``, tons - eq. 2.10 (K <= 1, braced) and
    eq. 2.14 (K >= 1, unbraced)."""
    return math.pi ** 2 * E * I / KL ** 2


def Cm_endmoments(M1: float, M2: float) -> float:
    """cl. 2.2.2 eq. 2.11 - ``Cm = 0.6 - 0.4 (M1/M2)`` for a member without
    transverse loading. M1/M2 positive for reverse curvature.

    With transverse loading use Cm = 0.85 (restrained ends) or 1.0 (simply
    supported), or Table 2.11 for a refined value.

    >>> Cm_endmoments(0, 100)
    0.6
    """
    ratio = 0.0 if M2 == 0 else M1 / M2
    return 0.6 - 0.4 * ratio


def B1(Cm: float, Pu: float, Pe1: float) -> float:
    """P-delta magnifier, cl. 2.2.2 eq. 2.9 - ``B1 = Cm/(1 - Pu/Pe1) >= 1.0``.

    A rigorous second-order analysis is recommended once B1 exceeds about 1.2.

    >>> round(B1(0.6, 50, 500), 4)
    1.0
    """
    if Pe1 <= 0:
        raise ValueError("Pe1 must be positive")
    if Pu >= Pe1:
        raise ValueError(
            f"Pu = {Pu:g} >= Pe1 = {Pe1:g}: the member is unstable in the "
            f"braced mode - increase the section or reduce KL")
    return max(Cm / (1 - Pu / Pe1), 1.0)


def B2_drift(sum_Pu: float, sum_H: float, delta_oh: float, L: float) -> float:
    """P-Delta magnifier from the storey drift, cl. 2.2.2 eq. 2.12::

        B2 = 1 / (1 - sum_Pu * delta_oh / (sum_H * L)) >= 1.0

    Args:
        sum_Pu: total factored gravity load above the level, tons.
        sum_H: sum of the storey horizontal forces producing delta_oh, tons.
        delta_oh: lateral inter-storey drift, cm.
        L: storey height, cm.
    """
    denom = 1 - sum_Pu * delta_oh / (sum_H * L)
    if denom <= 0:
        raise ValueError("B2 is negative or infinite - the storey is unstable")
    return max(1 / denom, 1.0)


def B2_euler(sum_Pu: float, sum_Pe2: float) -> float:
    """P-Delta magnifier from the Euler loads, cl. 2.2.2 eq. 2.13.

    ``sum_Pe2`` excludes leaning columns - they do not resist sidesway.
    """
    if sum_Pu >= sum_Pe2:
        raise ValueError("sum_Pu >= sum_Pe2: the storey is unstable")
    return max(1 / (1 - sum_Pu / sum_Pe2), 1.0)


def beam_column(Pu: float, phiPn: float, Mux: float, phiMnx: float,
                Muy: float = 0.0, phiMny: float = 1.0) -> Check:
    """Combined bending and axial force, cl. 7.1 eq. 7.1a / 7.1b.

    ``Mu`` must already include the second-order effects - either from a
    second-order analysis or via :func:`B1` and :func:`B2_drift` (cl. 2.2.2).

    Works for compression (phiPn from cl. 4.2, phi_c = 0.80) and for tension
    (phiPn from cl. 3.1, phi_t = 0.85) - cl. 7.3 uses the same equations.

    >>> c = beam_column(50, 200, 1000, 2500)
    >>> round(c.ratio, 4)
    0.6056
    """
    if phiPn <= 0 or phiMnx <= 0:
        raise ValueError("capacities must be positive")
    p = Pu / phiPn
    mx = Mux / phiMnx
    my = Muy / phiMny if Muy else 0.0
    if p >= 0.20:
        ratio = p + (8 / 9) * (mx + my)
        gov = "eq. 7.1a  (Pu/phiPn >= 0.20)"
    else:
        ratio = p / 2 + (mx + my)
        gov = "eq. 7.1b  (Pu/phiPn < 0.20)"
    return Check("Beam-column interaction", "7.1", ratio, 1.0,
                 {"Pu/phiPn": p, "Mux/phiMnx": mx, "Muy/phiMny": my}, gov)
