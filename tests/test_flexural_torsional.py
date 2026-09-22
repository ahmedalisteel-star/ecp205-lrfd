"""Flexural-torsional buckling, cl. 4.3, and modified slenderness, cl. 4.4.1.3.

ECP 205 prints no worked example for cl. 4.3, so these tests pin the
implementation three ways:

- worked examples whose every intermediate value is written out below, from the
  printed equations, so they can be checked by hand against the code;
- identities between the code's own equations - eq. 4.20 with xo = 0 must
  reproduce eq. 4.19, and cl. 4.3.3 on a doubly symmetric I must reproduce
  cl. 4.2 when flexure governs;
- the exact wording of the provisions: Table 4.1 strict inequalities, the 40 %
  gusset reduction as an alternative to Chapter 7, phi_c = 0.80.
"""

import math

import pytest

from ecp205 import (
    PHI,
    E,
    G,
    ISection,
    Tee,
    compression,
    compression_ft,
    compression_single_angle,
    compression_tee,
    double_angle,
    fe_unsymmetric,
    get,
    get_angle,
    modified_slenderness,
    table_4_1_exempt,
    tee_from,
)

FY = 2.4   # St 37


# ------------------------------------------------------------ cl. 4.4.1.3
def test_eq_4_25_lacing():
    # sqrt(80^2 + (100/2)^2) = sqrt(6400 + 2500)
    assert modified_slenderness(80, 100, 2.0, "lacing") == pytest.approx(
        math.sqrt(8900))


def test_eq_4_26_battens():
    # sqrt(80^2 + (1.25 * 100/2)^2) = sqrt(6400 + 3906.25)
    assert modified_slenderness(80, 100, 2.0, "battens") == pytest.approx(
        math.sqrt(10306.25))


def test_modified_slenderness_zero_spacing_is_unmodified():
    assert modified_slenderness(80, 0, 2.0, "battens") == pytest.approx(80)


def test_modified_slenderness_requires_connector_type():
    with pytest.raises(ValueError, match="lacing"):
        modified_slenderness(80, 100, 2.0, "stitch plates")


# ------------------------------------------------------------ Table 4.1
def _tee(bf_over_d, tf_over_tw, rolled):
    d, tw = 20.0, 1.0
    return Tee("T", d=d, bf=bf_over_d * d, tf=tf_over_tw * tw, tw=tw,
               rolled=rolled)


@pytest.mark.parametrize("bf_d,tf_tw,rolled,exempt", [
    (0.60, 1.20, True, True),
    (0.50, 1.20, True, False),     # > 0.5 is strict
    (0.60, 1.10, True, False),     # > 1.10 is strict (rolled)
    (0.60, 1.11, True, True),
    (0.60, 1.20, False, False),    # built-up needs > 1.25
    (0.60, 1.25, False, False),    # strict
    (0.60, 1.26, False, True),
    (0.40, 2.00, True, False),     # both ratios must pass
])
def test_table_4_1(bf_d, tf_tw, rolled, exempt):
    assert table_4_1_exempt(_tee(bf_d, tf_tw, rolled)) is exempt


def test_rolled_cut_tees_are_rolled():
    assert tee_from("IPE 300").rolled


# ------------------------------------------------------------ cl. 4.3.1
def test_worked_example_welded_tee_ft_governs():
    """Welded tee, bf = 20, tf = 1.2, d = 30, tw = 1.6 cm, sharp corners, St 37,
    KLx = KLy = 300 cm.

    Stem d/tw = 18.75 <= 30/sqrt(2.4) = 19.36, not slender (Table 2.12d).
    Built-up with tf/tw = 0.75, not > 1.25, so Table 4.1 does not exempt it.

        A    = 20(1.2) + 28.8(1.6)                      = 70.08 cm2
        yc   = (24(0.6) + 46.08(15.6)) / 70.08          = 10.4630 cm
        yo   = 0.6 - 10.4630                            = -9.8630 cm
        Ix   = 6738.61 cm4
        Iy   = 1.2(20^3)/12 + 28.8(1.6^3)/12            = 809.830 cm4
        J    = (20(1.2^3) + 29.4(1.6^3)) / 3            = 51.661 cm4
        KLy/ry = 300 / sqrt(809.83/70.08)               = 88.251
        lc   = 88.251/pi sqrt(2.4/2100)                 = 0.9497  (<= 1.1)
        Fcry = 2.4 (1 - 0.384 (0.9497^2))               = 1.5689  eq. 4.2
        ro^2 = 9.8630^2 + (6738.61 + 809.83)/70.08      = 204.991 cm2
        Fcrz = 810 (51.661) / (70.08 (204.991))         = 2.9128  eq. 4.13
        H    = 1 - 9.8630^2 / 204.991                   = 0.5254  eq. 4.14
        Fcrft                                           = 1.1840  eq. 4.12
        Fcrx on KLx/rx = 30.49                          = 2.3001  (does not govern)
        phi Pn = 0.80 (70.08) (1.1840)                  = 66.38 t eq. 4.11
    """
    t = Tee("WT", d=30, bf=20, tf=1.2, tw=1.6, rolled=False)
    c = compression_tee(t, KLx=300, KLy=300, Fy=FY, Pu=40)
    v = c.values
    assert t.A == pytest.approx(70.08)
    assert t.y0 == pytest.approx(-9.8630, abs=1e-4)
    assert t.Iy == pytest.approx(809.830, abs=1e-3)
    assert t.J == pytest.approx(51.661, abs=1e-3)
    assert v["KLy/ry"] == pytest.approx(88.251, abs=1e-3)
    assert v["Fcry"] == pytest.approx(1.5689, abs=1e-4)
    assert v["ro^2"] == pytest.approx(204.991, abs=1e-3)
    assert v["Fcrz"] == pytest.approx(2.9128, abs=1e-4)
    assert v["H"] == pytest.approx(0.5254, abs=1e-4)
    assert v["Fcrft"] == pytest.approx(1.1840, abs=1e-4)
    assert v["Fcrx"] == pytest.approx(2.3001, abs=1e-4)
    assert c.capacity == pytest.approx(66.38, abs=0.01)
    assert "eq. 4.12" in c.governing
    assert c.clause == "4.3.1"


def test_worked_example_double_angle_battened():
    """2L 100x100x10, 10 mm gusset, battens at 100 cm, St 37, KL = 300 cm.

        KLy/ry  = 300 / 4.5013                     = 66.647
        (KL/r)m = sqrt(66.647^2 + (1.25 (100)/1.9523)^2) = 92.420  eq. 4.26
        Fcry on (KL/r)m                            = 1.4885  eq. 4.2
        Fcrft                                      = 1.4373  eq. 4.12
        KLx/rx  = 300 / 3.0370 = 98.78  ->  Fcrx   = 1.3587  eq. 4.2 (governs)
        phi Pn  = 0.80 (38.309) (1.3587)           = 41.64 t
    """
    p = double_angle("L 100x100x10", gap=1.0)
    c = compression_tee(p, 300, 300, FY, lz=100, connectors="battens")
    v = c.values
    assert v["KLy/ry"] == pytest.approx(66.647, abs=1e-3)
    assert v["(KL/r)m"] == pytest.approx(92.420, abs=1e-3)
    assert v["Fcry"] == pytest.approx(1.4885, abs=1e-4)
    assert v["Fcrft"] == pytest.approx(1.4373, abs=1e-4)
    assert v["Fcrx"] == pytest.approx(1.3587, abs=1e-4)
    assert c.capacity == pytest.approx(41.64, abs=0.01)
    assert "about x" in c.governing


def test_double_angle_modified_slenderness_lowers_fcry():
    """Cl. 4.4.1.3 applies to cl. 4.3: wider connector spacing, lower Fcry."""
    p = double_angle("L 80x80x8", gap=1.0)
    close = compression_tee(p, 250, 250, FY, lz=40, connectors="lacing")
    wide = compression_tee(p, 250, 250, FY, lz=120, connectors="lacing")
    batt = compression_tee(p, 250, 250, FY, lz=120, connectors="battens")
    assert wide.values["Fcry"] < close.values["Fcry"]
    assert batt.values["Fcry"] < wide.values["Fcry"]
    assert batt.values["Fcrft"] < wide.values["Fcrft"]


def test_double_angle_requires_connectors():
    p = double_angle("L 80x80x8", gap=1.0)
    with pytest.raises(ValueError, match="4.4.1.3"):
        compression_tee(p, 250, 250, FY)


def test_fcrz_has_no_xo_and_no_cw():
    """Eq. 4.13/4.14: ro^2 = yo^2 + (Ix+Iy)/A and Fcrz = GJ/(A ro^2) - neither
    Cw nor xo appears."""
    t = Tee("WT", d=30, bf=20, tf=1.2, tw=1.6, rolled=False, Cw=1e9)
    v = compression_tee(t, 300, 300, FY).values
    ro2 = t.y0 ** 2 + (t.Ix + t.Iy) / t.A
    assert v["ro^2"] == pytest.approx(ro2)
    assert v["Fcrz"] == pytest.approx(G * t.J / (t.A * ro2))


def test_fcrft_never_exceeds_fcry_or_fcrz():
    for name in ("L 60x60x6", "L 100x100x10", "L 150x100x12"):
        v = compression_tee(double_angle(name, 1.0), 200, 200, FY, lz=60,
                            connectors="battens").values
        assert v["Fcrft"] <= min(v["Fcry"], v["Fcrz"]) + 1e-12


def test_exempt_tee_uses_flexural_buckling_both_axes():
    t = tee_from("HEB 200")
    assert table_4_1_exempt(t)
    c = compression_tee(t, 300, 300, FY)
    assert "Fcrft" not in c.values
    assert "Table 4.1" in c.note
    lower = min(compression(t.A, t.rx, 300, FY).capacity,
                compression(t.A, t.ry, 300, FY).capacity)
    assert c.capacity == pytest.approx(lower)


def test_slender_tee_stem_refused():
    """1/2 IPE 400: stem d/tw = 20/0.86 = 23.3 > 30/sqrt(2.4) = 19.4. Cl. 4.3.1
    has no Q, so the check is refused and cl. 4.3.3 is named instead."""
    with pytest.raises(ValueError, match="Table 2.12d.*compression_ft"):
        compression_tee(tee_from("IPE 400"), 300, 300, FY)


def test_slender_double_angle_refused():
    """L 100x100x6: b/t = 16.7 > 23/sqrt(2.4) = 14.8."""
    with pytest.raises(ValueError, match="b/t"):
        compression_tee(double_angle("L 100x100x6", 1.0), 200, 200, FY,
                        lz=60, connectors="battens")


def test_slender_tee_via_cl_4_3_3_with_q():
    c = compression_ft(tee_from("IPE 400"), 300, 300, 300, FY, Q=0.9)
    assert c.clause == "4.3.3" and c.values["Q"] == 0.9


def test_cl_4_3_3_refuses_double_angles():
    with pytest.raises(TypeError, match="4.4.1.3"):
        compression_ft(double_angle("L 80x80x8", 1.0), 200, 200, 200, FY)


def test_tee_rejects_q_below_one():
    with pytest.raises(ValueError, match="4.3.3"):
        compression_tee(tee_from("IPE 400"), 300, 300, FY, Q=0.9)


def test_phi_is_ecp_not_aisc():
    """phi_c = 0.80 (cl. 4.2.1), not AISC's 0.90."""
    t = Tee("WT", d=30, bf=20, tf=1.2, tw=1.6, rolled=False)
    c = compression_tee(t, 300, 300, FY)
    assert PHI["compression"] == 0.80
    assert c.capacity == pytest.approx(0.80 * t.A * c.values["Fcr"])


def test_slenderness_limit_note():
    p = double_angle("L 50x50x5", gap=1.0)
    c = compression_tee(p, 400, 400, FY, lz=100, connectors="battens")
    assert "Table 2.3" in c.note


# ------------------------------------------------------------ cl. 4.3.2
def test_single_angle_uses_rv():
    a = get_angle("L 100x75x8")
    c = compression_single_angle(a, 200, FY)
    ref = compression(a.A, a.rv, 200, FY)
    assert c.capacity == pytest.approx(ref.capacity)
    assert c.values["KL/r"] == pytest.approx(200 / a.rv)
    assert c.clause == "4.3.2"


def test_single_angle_slender_uses_q():
    a = get_angle("L 100x100x6")
    assert (compression_single_angle(a, 200, FY, Q=0.9).capacity
            < compression_single_angle(a, 200, FY).capacity)


def test_gusset_reduction_is_40_percent():
    a = get_angle("L 80x80x8")
    base = compression_single_angle(a, 150, FY)
    red = compression_single_angle(a, 150, FY, gusset="reduce")
    assert red.capacity == pytest.approx(0.6 * base.capacity)
    assert "40 %" in red.note


def test_gusset_chapter7_is_an_alternative_not_an_addition():
    """Chapter 7 route: the axial strength is not reduced; the note sends the
    eccentricity moment to the beam-column check."""
    a = get_angle("L 80x80x8")
    base = compression_single_angle(a, 150, FY)
    ch7 = compression_single_angle(a, 150, FY, gusset="chapter7")
    assert ch7.capacity == pytest.approx(base.capacity)
    assert "Chapter 7" in ch7.note


def test_gusset_option_validated():
    with pytest.raises(ValueError, match="gusset"):
        compression_single_angle(get_angle("L 80x80x8"), 150, FY, gusset="both")


# ------------------------------------------------------------ cl. 4.3.3
def test_eq_4_20_with_xo_zero_reproduces_eq_4_19():
    """With xo = 0 the cubic factors as (Fe - Fex)(eq. 4.19's quadratic), so
    its lowest root is min(Fex, eq. 4.19)."""
    Fex, Fey, Fez, yo, ro2 = 7.9, 3.6, 7.6, -3.85, 64.7
    H = 1 - yo ** 2 / ro2
    s = Fey + Fez
    fe19 = s / (2 * H) * (1 - math.sqrt(1 - 4 * Fey * Fez * H / s ** 2))
    assert fe_unsymmetric(Fex, Fey, Fez, 0.0, yo, ro2) == pytest.approx(
        min(Fex, fe19), rel=1e-10)


@pytest.mark.parametrize("args", [
    (6.2, 1.33, 7.17, -2.58, -1.66, 23.97),
    (3.0, 2.0, 1.5, 1.0, 2.0, 12.0),
    (10.0, 10.0, 10.0, 1.0, 1.0, 5.0),
])
def test_eq_4_20_root_is_a_root_and_lowest(args):
    Fex, Fey, Fez, xo, yo, ro2 = args
    Fe = fe_unsymmetric(*args)

    def cubic(f):
        return ((f - Fex) * (f - Fey) * (f - Fez)
                - f ** 2 * (f - Fey) * xo ** 2 / ro2
                - f ** 2 * (f - Fex) * yo ** 2 / ro2)
    assert cubic(Fe) == pytest.approx(0, abs=1e-9 * Fex * Fey * Fez)
    assert 0 < Fe <= min(Fex, Fey, Fez)
    assert all(cubic(Fe * k) * cubic(1e-12) > 0 for k in (0.2, 0.5, 0.9))


def test_doubly_symmetric_matches_cl_4_2_when_flexure_governs():
    b = get("HEB 200")
    c = compression_ft(b, 400, 400, 400, FY)
    assert c.capacity == pytest.approx(compression(b.A, b.ry, 400, FY).capacity)
    assert "flexural about y" in c.governing


def test_eq_4_18_torsional():
    """A short cruciform-like section: tiny Cw and J, stocky flexurally, so the
    torsional mode of eq. 4.18 governs."""
    s = ISection("X", d=30, bf=30, tf=1.0, tw=1.0, J=1.0, Cw=1.0)
    c = compression_ft(s, 50, 50, 500, FY)
    Fe = (math.pi ** 2 * E * 1.0 / 500 ** 2 + G * 1.0) / (s.Ix + s.Iy)
    assert c.values["Fe torsional"] == pytest.approx(Fe)
    assert "eq. 4.18" in c.governing


def test_eq_4_24_carries_xo():
    a = get_angle("L 150x90x10")
    v = compression_ft(a, 200, 200, 200, FY).values
    assert v["uo"] != 0 and v["vo"] != 0
    assert v["ro^2"] == pytest.approx(a.u0 ** 2 + a.v0 ** 2
                                      + (a.Iu + a.Iv) / a.A)
    assert v["H"] == pytest.approx(1 - (a.u0 ** 2 + a.v0 ** 2) / v["ro^2"])
    assert v["Fez"] == pytest.approx(
        (math.pi ** 2 * E * a.Cw / 200 ** 2 + G * a.J) / (a.A * v["ro^2"]))


def test_unequal_angle_ft_below_pure_flexure():
    a = get_angle("L 100x75x8")
    c = compression_ft(a, 200, 200, 200, FY)
    assert "eq. 4.20" in c.governing
    assert c.values["Fe"] < c.values["Fev"]


def test_equal_angle_symmetry_axis_is_u():
    a = get_angle("L 100x100x10")
    c = compression_ft(a, 200, 200, 200, FY)
    assert c.values["symmetry axis"] == "u"
    assert "Fe flexural-torsional" in c.values
    # here flexure about v governs, so 4.3.3 agrees with 4.3.2
    assert c.capacity == pytest.approx(
        compression_single_angle(a, 200, FY).capacity)


def test_cl_4_3_3_q_enters_eq_4_15():
    a = get_angle("L 100x100x6")
    assert (compression_ft(a, 100, 100, 100, FY, Q=0.85).capacity
            < compression_ft(a, 100, 100, 100, FY).capacity)
