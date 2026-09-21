"""Tests that pin the package to values printed in ECP 205-LRFD itself.

These are the tests that matter: they check the implementation against the
code's own tabulated numbers, not against the implementation's own output.
"""

import math

import pytest

from ecp205 import (
    PHI,
    E,
    G,
    bolt_pretension,
    classify,
    flange_limits,
    flexure_I,
    get,
    installation_torque,
    load_combinations,
    slip_resistance,
    steel,
    tube_limits,
    verify,
    web_limits,
)


# ------------------------------------------------------------ ch. 1 materials
def test_elastic_constants():
    """cl. 1.3.2."""
    assert E == 2100.0
    assert G == 810.0


@pytest.mark.parametrize("grade,t,expected", [
    ("St 37", 10, (2.40, 3.70)),
    ("St 37", 40, (2.40, 3.70)),      # boundary: 40 mm is still the thin band
    ("St 37", 41, (2.15, 3.40)),
    ("St 44", 20, (2.80, 4.40)),
    ("St 44", 60, (2.55, 4.10)),
    ("St 52", 20, (3.60, 5.20)),
    ("St 52", 60, (3.35, 4.90)),
])
def test_table_1_1(grade, t, expected):
    """Table 1.1, cl. 1.3.3 - every cell."""
    assert steel(grade, t) == expected


def test_steel_rejects_unknown_grade():
    with pytest.raises(ValueError, match="unknown grade"):
        steel("S355", 10)


def test_steel_rejects_thickness_beyond_table():
    with pytest.raises(ValueError, match="100 mm"):
        steel("St 37", 120)


def test_resistance_factors():
    """The phi values quoted through chapters 3 to 10."""
    assert PHI["tension_yield"] == 0.85
    assert PHI["tension_fracture"] == 0.70
    assert PHI["compression"] == 0.80       # NOT the AISC 0.90
    assert PHI["flexure"] == 0.85
    assert PHI["shear"] == 0.85
    assert PHI["web_yielding"] == 0.95
    assert PHI["web_crippling"] == 0.70
    assert PHI["weld"] == 0.70


# ------------------------------------------------------------ ch. 1 loads
def test_load_combinations_eq_1_1_to_1_6():
    """cl. 1.4.1 - the wind factor is 1.3, not the ASCE 7 value."""
    c = load_combinations(D=10, L=15, Lr=3, W=8, EQ=5)
    assert c["1.1  1.4D"] == pytest.approx(14.0)
    assert c["1.2  1.2D+1.6L+0.5Lr"] == pytest.approx(12 + 24 + 1.5)
    assert c["1.4a 1.2D+1.3W+0.5L+0.5Lr"] == pytest.approx(12 + 10.4 + 7.5 + 1.5)
    assert c["1.6a 0.9D+1.3W"] == pytest.approx(9 + 10.4)
    assert c["1.6c 0.9D+1.0EQ"] == pytest.approx(9 + 5)


def test_wind_and_eq_reversal_present():
    """Combinations 1.3 to 1.6 must be run both ways."""
    c = load_combinations(D=10, W=8, EQ=5)
    assert c["1.6b 0.9D-1.3W"] == pytest.approx(9 - 10.4)
    assert c["1.6d 0.9D-1.0EQ"] == pytest.approx(9 - 5)


def test_heavy_live_load_factor():
    """cl. 1.4.1 - L factor becomes 1.0 where L > 500 kg/m2."""
    normal = load_combinations(D=10, L=20)["1.5a 1.2D+1.0EQ+0.5L"]
    heavy = load_combinations(D=10, L=20, heavy_live=True)["1.5a 1.2D+1.0EQ+0.5L"]
    assert normal == pytest.approx(12 + 10)
    assert heavy == pytest.approx(12 + 20)


# ------------------------------------------------------------ Table 2.12
def test_table_2_12a_web_limits():
    """Table 2.12a, Fy in t/cm2."""
    Fy = 2.4
    s = math.sqrt(Fy)
    lp, lr = web_limits(Fy, alpha=0.5, psi=-1.0)          # pure bending
    assert lp == pytest.approx(127 / s)
    assert lr == pytest.approx(222 / s)
    lp, lr = web_limits(Fy, alpha=1.0, psi=1.0)           # pure compression
    assert lp == pytest.approx(58 / s)
    assert lr == pytest.approx(64 / s)


def test_table_2_12a_combined_bending_compression():
    """The alpha and psi branches of Table 2.12a."""
    Fy, s = 2.4, math.sqrt(2.4)
    assert web_limits(Fy, alpha=0.75)[0] == pytest.approx((699 / s) / (13 * 0.75 - 1))
    assert web_limits(Fy, alpha=0.3)[0] == pytest.approx((63.6 / 0.3) / s)
    assert web_limits(Fy, psi=0.0)[1] == pytest.approx((222 / s) / 2)
    assert web_limits(Fy, psi=-2.0)[1] == pytest.approx(
        111 * 3 * math.sqrt(2) / s)


def test_table_2_12c_flange_limits():
    """Table 2.12c - rolled and welded outstand flanges differ."""
    Fy, s = 2.4, math.sqrt(2.4)
    assert flange_limits(Fy, rolled=True) == pytest.approx((16.9 / s, 33 / s))
    assert flange_limits(Fy, rolled=False) == pytest.approx((15.3 / s, 28 / s))


def test_table_2_12d_tube_limits_divide_by_Fy_not_sqrt():
    """Table 2.12d - the tube limits use Fy, not sqrt(Fy)."""
    assert tube_limits(2.4) == pytest.approx((165 / 2.4, 211 / 2.4))


# ------------------------------------------------------------ Table 8.3
@pytest.mark.parametrize("d_mm,T_tabulated", [
    (12, 5.29), (16, 9.89), (20, 15.43), (22, 19.08),
    (24, 22.23), (27, 28.91), (30, 35.34), (36, 51.47),
])
def test_table_8_3_pretension(d_mm, T_tabulated):
    """Table 8.3 pretension force T = 0.7 Fyb As for grade 10.9 (eq. 8.7).

    This is the check that the whole ton/cm unit system is right: the values
    are read straight off the code's own table.
    """
    assert bolt_pretension(d_mm, "10.9") == pytest.approx(T_tabulated, rel=0.01)


@pytest.mark.parametrize("d_mm,Ps_tabulated", [
    (16, 3.16), (20, 4.93), (24, 7.11), (30, 11.30),
])
def test_table_8_3_slip_resistance(d_mm, Ps_tabulated):
    """Table 8.3 slip resistance, St 37/44 (mu = 0.4), ordinary steelwork,
    case I (F.S. = 1.25), one friction surface - eq. 8.8."""
    ps = slip_resistance(d_mm, "10.9", surface="B", case="I").capacity
    assert ps == pytest.approx(Ps_tabulated, rel=0.01)


@pytest.mark.parametrize("d_mm,Ma_kgm", [
    (16, 31), (20, 62), (24, 107), (30, 213), (36, 372),
])
def test_table_8_3_installation_torque(d_mm, Ma_kgm):
    """Table 8.3 torque Ma = 0.2 d T (eq. 8.17).

    The package returns t.cm; Table 8.3 prints kg.m, which is 10x larger.
    """
    # Table 8.3 rounds Ma to whole kg.m, and its own rounding of T carries
    # through - M36 is 370.6 against 372 tabulated, 0.4% out.
    assert installation_torque(d_mm, "10.9") * 10 == pytest.approx(
        Ma_kgm, rel=0.01, abs=1.0)


# ------------------------------------------------------------ ch. 5 constants
def test_Lp_is_80_ry_over_sqrt_Fy():
    """eq. 5.4 - dimensional, not the AISC 1.76 ry sqrt(E/Fy)."""
    sec = get("IPE 400")
    Fy = 2.4
    chk = flexure_I(sec, Fy, Lb=100)
    assert chk.values["Lp"] == pytest.approx(80 * sec.ry / math.sqrt(Fy))


def test_FL_rolled_vs_welded():
    """cl. 5.1.3.1 - FL = 0.75 Fy rolled, 0.60 Fy welded."""
    from ecp205 import ISection
    rolled = get("IPE 400")
    welded = ISection("plate girder", d=40, bf=18, tf=1.35, tw=0.86, r=0.0,
                      rolled=False)
    assert flexure_I(rolled, 2.4, Lb=100).values["FL"] == pytest.approx(1.8)
    assert flexure_I(welded, 2.4, Lb=100).values["FL"] == pytest.approx(1.44)


def test_shear_area_is_overall_depth():
    """cl. 5.2.1 - Aw = d * tw, using the OVERALL depth."""
    sec = get("IPE 400")
    assert sec.Aw == pytest.approx(sec.d * sec.tw)


# ------------------------------------------------------------ catalogue
def test_catalogue_is_self_consistent():
    """Every profile's A, Ix and Zx must agree with its own plate geometry."""
    assert verify(tol=0.03) == []


def test_known_section_properties():
    """Spot-check against the standard profile tables."""
    ipe400 = get("IPE 400")
    assert ipe400.Sx == pytest.approx(1156, rel=0.01)
    assert ipe400.ry == pytest.approx(3.95, rel=0.01)
    assert ipe400.mass == pytest.approx(66.3, rel=0.01)
    assert classify(ipe400, 2.4)["class"] == 1
