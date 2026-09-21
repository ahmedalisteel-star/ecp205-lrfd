"""Behavioural tests on the member strength functions."""

import math

import pytest

from ecp205 import (
    B1,
    Cb_ends,
    Cb_quarter,
    Cm_endmoments,
    ISection,
    beam_column,
    block_shear,
    compression,
    deflection_check,
    drift_check,
    flexure_I,
    get,
    shear_lag_U,
    shear_web,
    tension,
    udl_deflection,
)


# ------------------------------------------------------------ tension
def test_tension_takes_the_lower_limit_state():
    c = tension(Ag=84.5, Ae=75.0, Fy=2.4, Fu=3.7)
    assert c.capacity == pytest.approx(0.85 * 2.4 * 84.5)
    assert "gross" in c.governing


def test_tension_fracture_governs_with_a_small_net_area():
    c = tension(Ag=84.5, Ae=45.0, Fy=2.4, Fu=3.7)
    assert c.capacity == pytest.approx(0.70 * 3.7 * 45.0)
    assert "fracture" in c.governing


def test_shear_lag_U_capped_at_0_9():
    """eq. 2.2 - U = 1 - x/l but never above 0.9."""
    assert shear_lag_U(0.01, 100) == 0.9
    assert shear_lag_U(2.5, 20) == pytest.approx(0.875)


def test_block_shear_picks_the_right_branch():
    """cl. 8.9.3 - eq. 8.35 when tension rupture dominates, else eq. 8.36."""
    # eq. 8.35 needs Fu*Ant >= 0.6*Fu*Anv, i.e. Ant >= 0.6*Anv
    a = block_shear(Agv=30, Anv=20, Agt=16, Ant=15, Fy=2.4, Fu=3.7)
    assert "8.35" in a.governing
    b = block_shear(Agv=60, Anv=50, Agt=4, Ant=2, Fy=2.4, Fu=3.7)
    assert "8.36" in b.governing


# ------------------------------------------------------------ compression
def test_compression_inelastic_branch():
    """eq. 4.2 for lambda_c <= 1.1."""
    c = compression(Ag=84.5, r=3.95, KL=400, Fy=2.4)
    lam_c = c.values["lambda_c"]
    assert lam_c <= 1.1
    assert c.values["Fcr"] == pytest.approx(2.4 * (1 - 0.384 * lam_c ** 2))


def test_compression_elastic_branch():
    """eq. 4.3 for lambda_c > 1.1."""
    c = compression(Ag=84.5, r=3.95, KL=900, Fy=2.4)
    lam_c = c.values["lambda_c"]
    assert lam_c > 1.1
    assert c.values["Fcr"] == pytest.approx(0.648 * 2.4 / lam_c ** 2)


def test_compression_branches_are_continuous_at_lambda_c_1_1():
    """The two curves must not jump at the changeover."""
    Fy = 2.4
    r, Ag = 3.95, 84.5
    # KL that makes lambda_c exactly 1.1
    KL = math.pi * r * math.sqrt(2100 / Fy) * 1.1
    below = compression(Ag, r, KL * 0.999, Fy).values["Fcr"]
    above = compression(Ag, r, KL * 1.001, Fy).values["Fcr"]
    assert below == pytest.approx(above, rel=0.005)


def test_compression_flags_slenderness_over_table_2_3():
    c = compression(Ag=20, r=2.0, KL=400, Fy=2.4)   # KL/r = 200
    assert "Table 2.3" in c.note


def test_compression_rejects_nonsense():
    with pytest.raises(ValueError):
        compression(Ag=50, r=0, KL=300, Fy=2.4)


# ------------------------------------------------------------ Cb
def test_Cb_uniform_moment():
    """M1 = 0 gives the base value 1.75."""
    assert Cb_ends(0, 100) == pytest.approx(1.75)


def test_Cb_capped_at_2_3():
    assert Cb_ends(100, 100) == 2.3


def test_Cb_sign_convention():
    """M1/M2 positive for reverse curvature raises Cb; single curvature lowers it."""
    assert Cb_ends(50, 100) > Cb_ends(0, 100) > Cb_ends(-50, 100)


def test_Cb_quarter_simply_supported_udl():
    """The classic value for a UDL on a simple span is about 1.14."""
    assert Cb_quarter(1.0, 0.75, 1.0, 0.75) == pytest.approx(1.136, rel=0.01)


# ------------------------------------------------------------ flexure
def test_flexure_three_regimes():
    sec, Fy = get("IPE 400"), 2.4
    short = flexure_I(sec, Fy, Lb=100)
    mid = flexure_I(sec, Fy, Lb=400)
    long = flexure_I(sec, Fy, Lb=900)
    assert "5.3" in short.governing
    assert "5.6" in mid.governing
    assert "5.11" in long.governing
    assert short.capacity > mid.capacity > long.capacity


def test_flexure_capped_at_Mp():
    """Neither eq. 5.6 nor eq. 5.11 may exceed Mp."""
    sec, Fy = get("IPE 400"), 2.4
    c = flexure_I(sec, Fy, Lb=250, Cb=2.3)
    assert c.values["Mn"] <= c.values["Mp"] * (1 + 1e-9)


def test_flexure_Mp_limited_to_1_5_My():
    """cl. 5.1.3.1 - Mp = Fy Z but not more than 1.5 My."""
    stocky = ISection("stub", d=10, bf=20, tf=3.0, tw=2.0, r=0.0)
    c = flexure_I(stocky, 2.4, Lb=10)
    assert c.values["Mp"] <= 1.5 * c.values["My"] * (1 + 1e-9)


def test_flexure_rejects_negative_Lb():
    with pytest.raises(ValueError):
        flexure_I(get("IPE 400"), 2.4, Lb=-1)


# ------------------------------------------------------------ shear
def test_shear_three_regimes():
    Fy = 2.4
    stocky = ISection("stocky", d=40, bf=18, tf=1.35, tw=1.5, r=2.1)
    assert "5.22" in shear_web(stocky, Fy).governing
    slender = ISection("slender", d=150, bf=30, tf=2.0, tw=0.8, r=0.0)
    assert "5.24" in shear_web(slender, Fy).governing


def test_shear_rejects_h_over_tw_above_260():
    """cl. 5.2.2 stops at 260 - beyond that it is a plate girder."""
    girder = ISection("deep", d=300, bf=40, tf=2.0, tw=0.8, r=0.0)
    with pytest.raises(ValueError, match="plate girder"):
        shear_web(girder, 2.4)


# ------------------------------------------------------------ beam-column
def test_beam_column_picks_eq_7_1a_above_0_2():
    c = beam_column(Pu=100, phiPn=200, Mux=500, phiMnx=2500)
    assert "7.1a" in c.governing
    assert c.ratio == pytest.approx(0.5 + (8 / 9) * 0.2)


def test_beam_column_picks_eq_7_1b_below_0_2():
    c = beam_column(Pu=20, phiPn=200, Mux=500, phiMnx=2500)
    assert "7.1b" in c.governing
    assert c.ratio == pytest.approx(0.05 + 0.2)


def test_beam_column_curves_meet_on_the_interaction_boundary():
    """eq. 7.1a and 7.1b meet at (Pu/phiPn, Mu/phiMn) = (0.2, 0.9).

    They are two different curves, so they only coincide on the boundary
    itself - not for an arbitrary moment ratio.
    """
    phiPn, phiMnx = 200.0, 2500.0
    a = beam_column(0.2 * phiPn * 1.0001, phiPn, 0.9 * phiMnx, phiMnx).ratio
    b = beam_column(0.2 * phiPn * 0.9999, phiPn, 0.9 * phiMnx, phiMnx).ratio
    assert a == pytest.approx(1.0, rel=1e-3)
    assert b == pytest.approx(1.0, rel=1e-3)


def test_B1_never_below_one():
    assert B1(0.4, 10, 1000) == 1.0


def test_B1_rejects_unstable_member():
    with pytest.raises(ValueError, match="unstable"):
        B1(1.0, 500, 400)


def test_Cm_end_moments():
    assert Cm_endmoments(0, 100) == pytest.approx(0.6)
    assert Cm_endmoments(-100, 100) == pytest.approx(1.0)


# ------------------------------------------------------------ serviceability
def test_deflection_limits_table_14_1():
    assert deflection_check(1.0, 600, "other_beams", "delta_max").capacity == \
        pytest.approx(600 / 160)
    assert deflection_check(1.0, 600, "brittle_finish", "delta2_max").capacity == \
        pytest.approx(600 / 300)
    assert deflection_check(1.0, 600, "crane_track", "delta_max").capacity == \
        pytest.approx(600 / 650)


def test_appearance_row_has_no_delta2_limit():
    with pytest.raises(ValueError):
        deflection_check(1.0, 600, "appearance", "delta2_max")


def test_drift_limits_table_14_2():
    assert drift_check(1.0, 400, "storey").capacity == pytest.approx(400 / 300)
    assert drift_check(1.0, 900, "portal_no_crane").capacity == \
        pytest.approx(900 / 140)
    assert drift_check(1.0, 3000, "building_total").capacity == \
        pytest.approx(3000 / 500)


def test_udl_deflection():
    """5 w L^4 / (384 E I)."""
    assert udl_deflection(0.024, 800, 23130) == pytest.approx(
        5 * 0.024 * 800 ** 4 / (384 * 2100 * 23130))
