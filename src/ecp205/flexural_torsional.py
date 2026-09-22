"""Flexural-torsional buckling of compression members - ECP 205 cl. 4.3 - and the
modified slenderness of built-up members, cl. 4.4.1.3.

- :func:`compression_tee` - tees and double angles, cl. 4.3.1, eq. 4.11-4.14,
  with the Table 4.1 exemption for tees and eq. 4.25/4.26 for double angles.
- :func:`compression_single_angle` - single angles, cl. 4.3.2.
- :func:`compression_ft` - the general method of cl. 4.3.3, eq. 4.15-4.24, for
  doubly symmetric, singly symmetric and unsymmetric sections.

phi_c = 0.80 throughout (cl. 4.2.1). Units: ton, cm, t/cm2.

Cl. 4.3.1 and cl. 4.3.3 define the polar radius differently, and the difference
is deliberate: eq. 4.13/4.14 use ``ro^2 = yo^2 + (Ix + Iy)/A`` because xo = 0 for
the singly symmetric sections cl. 4.3.1 covers, while eq. 4.23/4.24 carry
``xo^2 + yo^2``. The two are kept as separate expressions here.
"""

from __future__ import annotations

import math
from typing import Optional

from .classification import angle_limit, tee_limit
from .core import LAMBDA_MAX, PHI, Check, E, G
from .members import compression
from .open_sections import Angle, DoubleAngle, Tee
from .sections import ISection

__all__ = ["modified_slenderness", "table_4_1_exempt", "compression_tee",
           "compression_single_angle", "compression_ft", "fe_unsymmetric"]


# ------------------------------------------------------------ cl. 4.4.1.3
def modified_slenderness(KL_r: float, lz: float, ri: float,
                         connectors: str) -> float:
    """Modified slenderness of an opened built-up member, cl. 4.4.1.3.

    ``connectors="lacing"`` (lacing with end batten plates), eq. 4.25::

        (KL/r)m = sqrt[ (KL/r)^2 + (lz/ri)^2 ]

    ``connectors="battens"`` (batten plates only), eq. 4.26::

        (KL/r)m = sqrt[ (KL/r)^2 + (1.25 lz/ri)^2 ]

    Args:
        KL_r: slenderness of the built-up member as a whole.
        lz: unsupported length of each separate part - the connector spacing, cm.
        ri: least radius of gyration of one part, cm.
        connectors: ``"lacing"`` or ``"battens"``.

    >>> round(modified_slenderness(80, 100, 1.95, "lacing"), 2)
    95.03
    >>> round(modified_slenderness(80, 100, 1.95, "battens"), 2)
    102.51
    """
    if connectors not in ("lacing", "battens"):
        raise ValueError("connectors must be 'lacing' (eq. 4.25) or "
                         "'battens' (eq. 4.26)")
    if lz < 0 or ri <= 0:
        raise ValueError("lz must be >= 0 and ri > 0")
    k = 1.0 if connectors == "lacing" else 1.25
    return math.sqrt(KL_r ** 2 + (k * lz / ri) ** 2)


# ------------------------------------------------------------ cl. 4.3.1
def table_4_1_exempt(tee: Tee) -> bool:
    """Table 4.1 - a T section need not be checked for flexural-torsional
    buckling when **both** ratios exceed the limits (strict inequalities)::

                          full flange width / depth    flange t / stem t
        built-up T                > 0.5                     > 1.25
        rolled T                  > 0.5                     > 1.10

    A tee cut from a rolled section counts as rolled.

    >>> from .open_sections import tee_from
    >>> table_4_1_exempt(tee_from("HEB 200"))     # bf/d = 2.0, tf/tw = 1.67
    True
    >>> table_4_1_exempt(tee_from("IPE 400"))     # bf/d = 0.9, tf/tw = 1.57
    True
    >>> table_4_1_exempt(Tee("WT", d=30, bf=20, tf=1.2, tw=1.0, rolled=False))
    False
    """
    t_limit = 1.10 if tee.rolled else 1.25
    return tee.bf / tee.d > 0.5 and tee.tf / tee.tw > t_limit


def _fcr_flexural(lam: float, Fy: float) -> tuple[float, float, str]:
    """Cl. 4.2.1 eq. 4.2-4.5 on a slenderness KL/r: returns (Fcr, lambda_c, eq)."""
    lc = lam / math.pi * math.sqrt(Fy / E)                  # eq. 4.4 / 4.5
    if lc <= 1.1:
        return Fy * (1 - 0.384 * lc ** 2), lc, "eq. 4.2"
    return 0.648 * Fy / lc ** 2, lc, "eq. 4.3"


def _slender_element(sec: Tee | DoubleAngle, Fy: float) -> Optional[str]:
    """The first element beyond its Table 2.12d non-compact limit, or None.

    Tee stem: d/tw against 30/sqrt(Fy). Angles of a double angle: b/t against
    23/sqrt(Fy), or (b+d)/2t against 17/sqrt(Fy) for unequal angles.
    """
    if isinstance(sec, Tee):
        lam, lim = sec.d / sec.tw, tee_limit(Fy)
        return None if lam <= lim else (
            f"stem d/tw = {lam:0.1f} > 30/sqrt(Fy) = {lim:0.1f}")
    a = sec.angle
    if a.equal:
        lam, lim, what = a.b / a.t, angle_limit(Fy), "b/t"
    else:
        lam, lim, what = (a.b + a.d) / (2 * a.t), angle_limit(Fy, True), "(b+d)/2t"
    return None if lam <= lim else f"angle {what} = {lam:0.2f} > {lim:0.2f}"


def compression_tee(sec: Tee | DoubleAngle, KLx: float, KLy: float,
                    Fy: float, Pu: float = 0.0, Q: float = 1.0,
                    lz: Optional[float] = None, connectors: Optional[str] = None,
                    member: str = "compression") -> Check:
    """Design compressive strength of a tee or double angle, cl. 4.3.1,
    phi_c = 0.80. y is the axis of symmetry.

    The lower of

    - flexural buckling about x, cl. 4.2.1 on KLx/rx, and
    - flexural-torsional buckling, eq. 4.11-4.14::

          Pn    = Ag Fcrft                                                (4.11)
          Fcrft = (Fcry + Fcrz)/(2H) [1 - sqrt(1 - 4 Fcry Fcrz H/(Fcry+Fcrz)^2)]
                                                                          (4.12)
          Fcrz  = G J / (A ro^2)                                          (4.13)
          ro^2  = yo^2 + (Ix + Iy)/A
          H     = 1 - yo^2/ro^2                                           (4.14)

      with Fcry from cl. 4.2.1 on the slenderness about the axis of symmetry.

    A tee meeting Table 4.1 is exempt from the flexural-torsional check; its
    strength is then the lower of flexural buckling about x and about y.

    For a **double angle**, cl. 4.4.1.3 applies to cl. 4.3 as well as 4.2:
    buckling about the axis of symmetry shears the connectors, so KLy/ry is
    replaced by (KL/r)m from eq. 4.25/4.26 before Fcry is computed. ``lz`` and
    ``connectors`` are therefore required for a double angle. Flexure about x
    moves both angles together and is not modified.

    Args:
        sec: a :class:`~ecp205.open_sections.Tee` or
            :class:`~ecp205.open_sections.DoubleAngle`.
        KLx, KLy: effective lengths about x and about the axis of symmetry y, cm.
        Fy: t/cm2.
        Pu: required compressive strength, tons.
        Q: must be 1.0. Eq. 4.12 takes Fcry from cl. 4.2.1, which has no Q.
            The section is classified against Table 2.12d and a slender one is
            refused; use :func:`compression_ft` (cl. 4.3.3, which carries Q).
        lz: connector spacing along a double angle, cm (cl. 4.4.1.3).
        connectors: ``"lacing"`` (eq. 4.25) or ``"battens"`` (eq. 4.26).
        member: Table 2.3 limit to report against.
    """
    if Q != 1.0:
        raise ValueError(
            "cl. 4.3.1 takes Fcry from cl. 4.2.1, which has no Q - for a section "
            "with slender elements use compression_ft() (cl. 4.3.3, eq. 4.15/4.16)")
    if KLx <= 0 or KLy <= 0:
        raise ValueError("KLx and KLy must be positive")
    is_pair = isinstance(sec, DoubleAngle)
    if not isinstance(sec, (Tee, DoubleAngle)):
        raise TypeError("compression_tee takes a Tee or a DoubleAngle")
    slender = _slender_element(sec, Fy)
    if slender:
        raise ValueError(
            f"{sec.name}: {slender} (Table 2.12d) - a slender section. Cl. 4.3.1 "
            "takes Fcry from cl. 4.2.1, which has no Q; use compression_ft() "
            "(cl. 4.3.3, eq. 4.15/4.16) with Q = Ae/Ag")

    lam_x = KLx / sec.rx
    lam_y = KLy / sec.ry
    values: dict = {"KLx/rx": lam_x, "KLy/ry": lam_y}
    if is_pair:
        if lz is None or connectors is None:
            raise ValueError(
                "cl. 4.4.1.3: a double angle needs the connector spacing lz and "
                "the connector type ('lacing' or 'battens') - its slenderness "
                "about the axis of symmetry is (KL/r)m, eq. 4.25/4.26")
        ri = sec.angle.rv
        lam_y = modified_slenderness(lam_y, lz, ri, connectors)
        values.update({"lz": lz, "ri": ri, "(KL/r)m": lam_y})

    Fcrx, lcx, eqx = _fcr_flexural(lam_x, Fy)
    Fcry, lcy, eqy = _fcr_flexural(lam_y, Fy)
    values.update({"Fcrx": Fcrx, "Fcry": Fcry})
    note = ""

    exempt = (not is_pair) and table_4_1_exempt(sec)
    if exempt:
        Fcr_t, gov_t = Fcry, f"flexural buckling about y, {eqy}"
        note = "Table 4.1: flexural-torsional buckling need not be checked"
    else:
        ro2 = sec.y0 ** 2 + (sec.Ix + sec.Iy) / sec.A          # 4.3.1, no xo
        Fcrz = G * sec.J / (sec.A * ro2)                       # eq. 4.13
        H = 1 - sec.y0 ** 2 / ro2                              # eq. 4.14
        s = Fcry + Fcrz
        Fcrft = s / (2 * H) * (1 - math.sqrt(1 - 4 * Fcry * Fcrz * H / s ** 2))
        values.update({"yo": sec.y0, "ro^2": ro2, "Fcrz": Fcrz, "H": H,
                       "Fcrft": Fcrft})
        Fcr_t, gov_t = Fcrft, f"flexural-torsional buckling, eq. 4.12 (Fcry {eqy})"

    if Fcrx <= Fcr_t:
        Fcr, gov = Fcrx, f"flexural buckling about x, {eqx}"
    else:
        Fcr, gov = Fcr_t, gov_t
    values.update({"Fcr": Fcr, "Ag": sec.A})
    chk = Check(f"Compression {sec.name}", "4.3.1", Pu,
                PHI["compression"] * sec.A * Fcr, values, gov, note)
    lam_max, limit = max(lam_x, lam_y), LAMBDA_MAX[member]
    if lam_max > limit:
        chk.note = "; ".join(filter(None, [chk.note, (
            f"KL/r = {lam_max:0.0f} exceeds the Table 2.3 limit of {limit} "
            f"for a {member} member")]))
    return chk


# ------------------------------------------------------------ cl. 4.3.2
def compression_single_angle(angle: Angle, KL: float, Fy: float, Pu: float = 0.0,
                             Q: float = 1.0, gusset: Optional[str] = None,
                             member: str = "compression") -> Check:
    """Design compressive strength of a single angle, cl. 4.3.2, phi_c = 0.80.

    Flexural buckling on the least radius of gyration r_v: eq. 4.1-4.3 for a
    non-compact angle, eq. 4.6-4.9 with Q for a slender one. Cl. 4.3.2 does not
    name the axis; r_v follows from cl. 4.2, where r is the radius of gyration
    about the axis of buckling, and for a single angle that is the minor
    principal axis. Cl. 4.3.3 (:func:`compression_ft`) may be used instead to
    consider flexural-torsional buckling directly.

    Single angles **connected to gusset plates** carry a moment from the load
    eccentricity. The code gives two alternatives, and they are alternatives -
    not both:

    - ``gusset="chapter7"``: design as a beam-column to Chapter 7. The capacity
      returned is the axial strength phi_c Pn only, for use in eq. 7.1.
    - ``gusset="reduce"``: reduce the design strength by 40 %.

    Args:
        angle: the angle.
        KL: effective length, cm.
        Q: Ae/Ag for a slender angle (cl. 4.2.2), 1.0 otherwise.
        gusset: None for a concentrically loaded angle, else ``"reduce"`` or
            ``"chapter7"``.

    >>> from .open_sections import get_angle
    >>> a = get_angle("L 100x100x10")
    >>> c0 = compression_single_angle(a, 200, 2.4)
    >>> c1 = compression_single_angle(a, 200, 2.4, gusset="reduce")
    >>> round(c1.capacity / c0.capacity, 2)
    0.6
    """
    if gusset not in (None, "reduce", "chapter7"):
        raise ValueError("gusset must be None, 'reduce' or 'chapter7'")
    base = compression(angle.A, angle.rv, KL, Fy, Q=Q, Pu=Pu, member=member)
    values = {"r": "r_v", **base.values}
    cap, note = base.capacity, base.note
    if gusset == "reduce":
        cap *= 0.6
        values["gusset reduction"] = 0.6
        note = "; ".join(filter(None, [
            "gusset-connected: design strength reduced by 40 % (cl. 4.3.2)", note]))
    elif gusset == "chapter7":
        note = "; ".join(filter(None, [
            "gusset-connected: axial strength only - check the eccentricity "
            "moment as a beam-column, Chapter 7 (cl. 4.3.2)", note]))
    return Check(f"Compression {angle.name}", "4.3.2", Pu, cap, values,
                 base.governing, note)


# ------------------------------------------------------------ cl. 4.3.3
def fe_unsymmetric(Fex: float, Fey: float, Fez: float, xo: float, yo: float,
                   ro2: float) -> float:
    """Lowest root of the cubic, cl. 4.3.3 eq. 4.20::

        (Fe-Fex)(Fe-Fey)(Fe-Fez) - Fe^2 (Fe-Fey)(xo/ro)^2
                                 - Fe^2 (Fe-Fex)(yo/ro)^2 = 0

    The lowest root lies in (0, min(Fex, Fey, Fez)]; it is bracketed by a scan
    from zero and then bisected.

    >>> round(fe_unsymmetric(10, 20, 30, 0, 0, 1), 6)      # uncoupled: min
    10.0
    """
    ax, ay = xo ** 2 / ro2, yo ** 2 / ro2

    def f(Fe: float) -> float:
        return ((Fe - Fex) * (Fe - Fey) * (Fe - Fez)
                - Fe ** 2 * (Fe - Fey) * ax - Fe ** 2 * (Fe - Fex) * ay)

    # f(0) = -Fex Fey Fez < 0 and the lowest root lies in (0, min]. Scan for
    # the first point where f stops being negative, then bisect - this finds
    # the lowest root even when two of Fex, Fey, Fez coincide and f touches
    # zero again at the upper bound.
    top = min(Fex, Fey, Fez)
    n = 2000
    lo, hi = 0.0, top
    for i in range(1, n + 1):
        x = top * i / n
        if f(x) >= 0:
            lo, hi = top * (i - 1) / n, x
            break
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _principal(sec) -> tuple[str, float, float, float, float, float, float, float]:
    """(name, A, Ix, Iy, xo, yo, J, Cw) about the principal axes."""
    if isinstance(sec, Angle):
        return sec.name, sec.A, sec.Iu, sec.Iv, sec.u0, sec.v0, sec.J, sec.Cw
    if isinstance(sec, Tee):
        return sec.name, sec.A, sec.Ix, sec.Iy, 0.0, sec.y0, sec.J, sec.Cw
    if isinstance(sec, ISection):
        return sec.name, sec.A, sec.Ix, sec.Iy, 0.0, 0.0, sec.J, sec.Cw
    raise TypeError(f"unsupported section type {type(sec).__name__}")


def compression_ft(sec, KLx: float, KLy: float, KLz: float, Fy: float,
                   Pu: float = 0.0, Q: float = 1.0,
                   member: str = "compression") -> Check:
    """Design compressive strength for torsional and flexural-torsional buckling,
    the general method of cl. 4.3.3, phi_c = 0.80.

    ::

        lambda_e sqrt(Q) <= 1.1 :  Fcr = Fy Q (1 - 0.384 Q lambda_e^2)    (4.15)
        lambda_e sqrt(Q) >  1.1 :  Fcr = 0.648 Fy / lambda_e^2            (4.16)
        lambda_e = sqrt(Fy / Fe)                                          (4.17)

        Fex = pi^2 E / (KxL/rx)^2                                         (4.21)
        Fey = pi^2 E / (KyL/ry)^2                                         (4.22)
        Fez = [pi^2 E Cw/(KzL)^2 + G J] / (A ro^2)                        (4.23)
        ro^2 = xo^2 + yo^2 + (Ix + Iy)/A
        H   = 1 - (xo^2 + yo^2)/ro^2                                      (4.24)

    Fe is the elastic buckling stress of the governing mode:

    - doubly symmetric (xo = yo = 0): torsional, eq. 4.18
      ``Fe = [pi^2 E Cw/(KzL)^2 + G J] / (Ix + Iy)``, against flexure about
      x and y;
    - singly symmetric about y (xo = 0): flexural-torsional, eq. 4.19
      ``Fe = (Fey+Fez)/(2H) [1 - sqrt(1 - 4 Fey Fez H/(Fey+Fez)^2)]``, against
      flexure about x;
    - unsymmetric: the lowest root of eq. 4.20, which already includes the
      flexural modes.

    Axes are the **principal** axes. For a single angle x is the major axis u and
    y the minor axis v, so ``KLx`` and ``KLy`` are the effective lengths about u
    and v. For a tee y is the axis of symmetry. A section
    symmetric about x instead (an equal angle, symmetric about u) has its axes
    relabelled internally so that y is the axis of symmetry, as eq. 4.19 is
    written. Results are labelled with the section's own axis names (u, v for
    an angle), and ``values["symmetry axis"]`` says which one is the axis of
    symmetry.

    Args:
        sec: ISection, Tee or Angle.
        KLx, KLy: flexural effective lengths about the principal axes, cm.
        KLz: effective length for torsional buckling, Kz L, cm.
        Q: Ae/Ag, cl. 4.2.2; 1.0 when no element exceeds lambda_r.

    Raises:
        TypeError: for a double angle. Cl. 4.4.1.3 modifies its slenderness in
            cl. 4.3 as well; use :func:`compression_tee`.
    """
    if isinstance(sec, DoubleAngle):
        raise TypeError(
            "cl. 4.4.1.3 applies to cl. 4.3: a double angle's slenderness about "
            "the axis of symmetry is (KL/r)m - use compression_tee()")
    if min(KLx, KLy, KLz) <= 0:
        raise ValueError("KLx, KLy and KLz must be positive")
    name, A, Ix, Iy, xo, yo, J, Cw = _principal(sec)
    tol = 1e-9 * math.sqrt((Ix + Iy) / A)
    xo = 0.0 if abs(xo) < tol else xo
    yo = 0.0 if abs(yo) < tol else yo
    n1, n2 = ("u", "v") if isinstance(sec, Angle) else ("x", "y")
    swapped = xo != 0 and yo == 0
    if swapped:
        n1, n2 = n2, n1
        # symmetric about x (an equal angle: the major axis u). Eq. 4.19 is
        # written with y as the axis of symmetry, so relabel the axes.
        Ix, Iy, xo, yo, KLx, KLy = Iy, Ix, yo, xo, KLy, KLx

    rx, ry = math.sqrt(Ix / A), math.sqrt(Iy / A)
    Fex = math.pi ** 2 * E / (KLx / rx) ** 2                   # eq. 4.21
    Fey = math.pi ** 2 * E / (KLy / ry) ** 2                   # eq. 4.22
    ro2 = xo ** 2 + yo ** 2 + (Ix + Iy) / A                    # 4.3.3, with xo
    tors = math.pi ** 2 * E * Cw / KLz ** 2 + G * J
    Fez = tors / (A * ro2)                                     # eq. 4.23
    H = 1 - (xo ** 2 + yo ** 2) / ro2                          # eq. 4.24
    values: dict = {"symmetry axis": n2 if (xo == 0) != (yo == 0) else "-",
                    f"KL{n1}/r{n1}": KLx / rx, f"KL{n2}/r{n2}": KLy / ry,
                    f"{n1}o": xo, f"{n2}o": yo, "ro^2": ro2,
                    f"Fe{n1}": Fex, f"Fe{n2}": Fey, "Fez": Fez, "H": H}

    if xo == 0 and yo == 0:
        Fet = tors / (Ix + Iy)                                 # eq. 4.18
        modes = [(Fex, f"flexural about {n1}, eq. 4.21"),
                 (Fey, f"flexural about {n2}, eq. 4.22"),
                 (Fet, "torsional, eq. 4.18")]
        values["Fe torsional"] = Fet
    elif xo == 0:
        s = Fey + Fez
        Feft = s / (2 * H) * (1 - math.sqrt(1 - 4 * Fey * Fez * H / s ** 2))
        modes = [(Fex, f"flexural about {n1}, eq. 4.21"),
                 (Feft, "flexural-torsional, eq. 4.19")]
        values["Fe flexural-torsional"] = Feft
    else:
        modes = [(fe_unsymmetric(Fex, Fey, Fez, xo, yo, ro2),
                  "flexural-torsional, eq. 4.20")]
    Fe, mode = min(modes)

    le = math.sqrt(Fy / Fe)                                    # eq. 4.17
    if le * math.sqrt(Q) <= 1.1:
        Fcr, eq = Fy * Q * (1 - 0.384 * Q * le ** 2), "eq. 4.15"
    else:
        Fcr, eq = 0.648 * Fy / le ** 2, "eq. 4.16"
    values.update({"Fe": Fe, "lambda_e": le, "Q": Q, "Fcr": Fcr, "Ag": A})
    chk = Check(f"Compression {name}", "4.3.3", Pu,
                PHI["compression"] * A * Fcr, values, f"{mode}; {eq}")
    lam_max, limit = max(KLx / rx, KLy / ry), LAMBDA_MAX[member]
    if lam_max > limit:
        chk.note = (f"KL/r = {lam_max:0.0f} exceeds the Table 2.3 limit of "
                    f"{limit} for a {member} member")
    return chk
