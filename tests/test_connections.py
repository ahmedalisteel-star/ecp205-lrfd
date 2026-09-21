"""Tests on the bolt, weld and concentrated-force provisions."""

import pytest

from ecp205 import (
    base_plate_bearing,
    bolt_bearing,
    bolt_combined,
    bolt_shear,
    bolt_tension,
    fillet_weld,
    fillet_weld_combined,
    flange_local_bending,
    groove_weld,
    long_joint_factor,
    min_fillet_size,
    packing_factor,
    prying_force,
    slip_resistance,
    web_crippling,
    web_local_yielding,
)


# ------------------------------------------------------------ bolts
def test_bolt_shear_grade_dependent_coefficient():
    """cl. 8.5.2 - 0.6 Fub for 8.8, but 0.5 Fub for 10.9 (eq. 8.3)."""
    assert bolt_shear(20, "8.8").values["k"] == 0.6
    assert bolt_shear(20, "10.9").values["k"] == 0.5
    assert bolt_shear(20, "4.8").values["k"] == 0.5


def test_bolt_shear_uses_stress_area_or_shank():
    threaded = bolt_shear(20, "8.8", threads_in_plane=True).capacity
    excluded = bolt_shear(20, "8.8", threads_in_plane=False).capacity
    assert excluded > threaded


def test_bolt_shear_scales_with_planes():
    single = bolt_shear(20, "8.8", n_planes=1).capacity
    double = bolt_shear(20, "8.8", n_planes=2).capacity
    assert double == pytest.approx(2 * single)


@pytest.mark.parametrize("e1_over_d,alpha", [
    (3.0, 2.4), (2.7, 2.0), (2.2, 1.6), (1.6, 1.2),
])
def test_bolt_bearing_table_8_2(e1_over_d, alpha):
    """Table 8.2 - alpha steps with the end distance."""
    d = 2.0   # M20
    c = bolt_bearing(20, sum_t=1.0, Fu=3.7, e1=e1_over_d * d)
    assert c.values["alpha"] == pytest.approx(min(alpha, 0.8 * e1_over_d))


def test_bolt_bearing_rejects_end_distance_below_minimum():
    """cl. 8.4.2 - e1 must be at least 1.5d."""
    with pytest.raises(ValueError, match="1.5d"):
        bolt_bearing(20, sum_t=1.0, Fu=3.7, e1=2.0)


def test_bolt_bearing_long_slot_reduction():
    plain = bolt_bearing(20, 1.0, 3.7, e1=6.0).capacity
    slot = bolt_bearing(20, 1.0, 3.7, e1=6.0, long_slot_perp=True).capacity
    assert slot == pytest.approx(0.8 * plain)


def test_high_grade_bolts_in_tension_must_be_pretensioned():
    """cl. 8.5.4b."""
    with pytest.raises(ValueError, match="pretensioned"):
        bolt_tension(20, "8.8", pretensioned=False)
    ok = bolt_tension(20, "8.8", pretensioned=True)
    assert ok.values["k"] == 0.80


def test_bolt_tension_snug_tight_uses_0_66_Fub():
    """cl. 8.5.4 eq. 8.5."""
    assert bolt_tension(20, "5.6").values["k"] == 0.66


def test_bolt_tension_adds_prying_to_demand():
    c = bolt_tension(20, "10.9", Ru=5.0, pretensioned=True, prying=2.0)
    assert c.demand == pytest.approx(7.0)


def test_bolt_combined_is_circular():
    """eq. 8.6 - the interaction is a circle, not a straight line."""
    c = bolt_combined(Rut=3, phiRnt=5, Ruv=4, phiRnv=5)
    assert c.ratio == pytest.approx(1.0)


def test_slip_only_for_high_grades():
    """cl. 8.6.1."""
    with pytest.raises(ValueError, match="8.8 and 10.9"):
        slip_resistance(20, "5.6")


def test_slip_surface_classes():
    """cl. 8.6.2.2 - class A 0.50, B 0.40, C 0.30."""
    a = slip_resistance(20, "10.9", surface="A").capacity
    b = slip_resistance(20, "10.9", surface="B").capacity
    c = slip_resistance(20, "10.9", surface="C").capacity
    assert a / b == pytest.approx(0.5 / 0.4)
    assert b / c == pytest.approx(0.4 / 0.3)


def test_slip_case_II_is_less_conservative():
    """F.S. drops from 1.25 to 1.05 for case II."""
    i = slip_resistance(20, "10.9", case="I").capacity
    ii = slip_resistance(20, "10.9", case="II").capacity
    assert ii > i


def test_slip_crane_cases_are_more_conservative():
    ordinary = slip_resistance(20, "10.9", case="I").capacity
    crane = slip_resistance(20, "10.9", case="crane_I").capacity
    assert crane < ordinary


def test_prying_force_eq_8_26():
    assert prying_force(b=5, a=4, tp=2.0, T_ext=12) == pytest.approx(8.415, rel=1e-3)


def test_prying_force_never_negative():
    """A thick plate gives no prying, not a negative force."""
    assert prying_force(b=2, a=10, tp=6.0, T_ext=10) == 0.0


def test_refined_prying_needs_w_and_As():
    with pytest.raises(ValueError, match="eq. 8.27"):
        prying_force(5, 4, 2.0, 12, refined=True)


def test_long_joint_factor_bounds():
    """eq. 8.37 - unity below 15d, floored at 0.75, rejected beyond 65d."""
    assert long_joint_factor(25, 20) == 1.0          # Lj <= 15d = 30 cm
    assert long_joint_factor(60, 20) == pytest.approx(0.925)
    assert long_joint_factor(100, 20) == pytest.approx(0.825)
    with pytest.raises(ValueError, match="65d"):
        long_joint_factor(200, 20)


def test_packing_factor_only_beyond_d_over_3():
    assert packing_factor(tp=0.5, d_mm=20) == 1.0
    assert packing_factor(tp=2.0, d_mm=20) < 1.0


def test_base_plate_bearing_partial_area_capped():
    """eq. 8.42 - never more than twice the full-area value."""
    full = base_plate_bearing(A1=400, fcu=0.25).capacity
    huge = base_plate_bearing(A1=400, fcu=0.25, A2=40000).capacity
    assert huge == pytest.approx(2 * full)


# ------------------------------------------------------------ welds
def test_fillet_weld_eq_9_7():
    """phi s (0.4 Fu) per unit length, St 37 -> 1.036 t/cm per cm of leg."""
    per_cm = fillet_weld(s=1.0, Fu=3.7, length=1.0).capacity
    assert per_cm == pytest.approx(0.7 * 0.4 * 3.7)
    assert per_cm == pytest.approx(1.036, rel=1e-3)


def test_fillet_weld_scales_with_length_and_leg():
    a = fillet_weld(0.6, 3.7, length=10).capacity
    b = fillet_weld(0.6, 3.7, length=20).capacity
    assert b == pytest.approx(2 * a)


def test_fillet_weld_combined_gets_the_ten_percent():
    """eq. 9.8 - the capacity is 0.77 s (0.4 Fu), not 0.70."""
    c = fillet_weld_combined(0.0, 0.0, 0.0, s=0.6, Fu=3.7)
    assert c.capacity == pytest.approx(0.77 * 0.6 * 0.4 * 3.7)


def test_fillet_weld_combined_von_mises_form():
    c = fillet_weld_combined(R_perp=1.0, R_par=1.0, R_shear_perp=0.0,
                             s=0.6, Fu=3.7)
    assert c.demand == pytest.approx((1 + 3) ** 0.5)


def test_groove_weld_shear_is_point_six():
    t, Fy = 1.0, 2.4
    assert groove_weld(t, Fy).capacity == pytest.approx(0.85 * t * Fy)
    assert groove_weld(t, Fy, shear=True).capacity == pytest.approx(
        0.85 * t * 0.6 * Fy)


@pytest.mark.parametrize("t_mm,s_mm", [(5, 3.0), (10, 5.0), (15, 6.0), (25, 8.0)])
def test_min_fillet_size_table_9_6(t_mm, s_mm):
    assert min_fillet_size(t_mm) == s_mm


# ------------------------------------------------------------ ch. 10
def test_flange_local_bending_halved_near_the_end():
    interior = flange_local_bending(tf=1.35, Fyf=2.4).capacity
    end = flange_local_bending(tf=1.35, Fyf=2.4, near_end=True).capacity
    assert end == pytest.approx(0.5 * interior)


def test_web_local_yielding_5k_vs_2_5k():
    """cl. 10.3 - the end-zone coefficient is 2.5k, not 5k."""
    interior = web_local_yielding(k=3.5, N=10, Fyw=2.4, tw=0.86).capacity
    end = web_local_yielding(k=3.5, N=10, Fyw=2.4, tw=0.86,
                             near_end=True).capacity
    assert interior > end


def test_web_crippling_end_zone_is_half():
    """0.3627 interior (eq. 10.4) vs 0.1813 at the end (eq. 10.5a).

    Only comparable while N/d <= 0.2, where both use the same 3(N/d) bracket.
    Above that the end uses eq. 10.5b with (4N/d - 0.2) instead, and the two
    are no longer a simple factor of two apart.
    """
    kw = dict(tw=0.86, tf=1.35, d=40, Fyw=2.4)
    interior = web_crippling(N=6, **kw).capacity            # N/d = 0.15
    end = web_crippling(N=6, near_end=True, **kw).capacity
    assert end == pytest.approx(0.5 * interior, rel=1e-3)


def test_web_crippling_end_brackets_differ_above_N_over_d_0_2():
    """Above N/d = 0.2 the end equation is 10.5b, not a halved 10.4."""
    kw = dict(tw=0.86, tf=1.35, d=40, Fyw=2.4)
    interior = web_crippling(N=10, **kw)                    # N/d = 0.25
    end = web_crippling(N=10, near_end=True, **kw)
    assert end.governing == "eq. 10.5b"
    assert end.capacity != pytest.approx(0.5 * interior.capacity, rel=1e-3)


def test_web_crippling_branch_at_N_over_d_0_2():
    """eq. 10.5a and 10.5b must meet at N/d = 0.2."""
    kw = dict(tw=0.86, tf=1.35, d=40, Fyw=2.4, near_end=True)
    a = web_crippling(N=0.2 * 40 * 0.999, **kw).capacity
    b = web_crippling(N=0.2 * 40 * 1.001, **kw).capacity
    assert a == pytest.approx(b, rel=0.01)


def test_web_crippling_improves_with_bearing_length():
    kw = dict(tw=0.86, tf=1.35, d=40, Fyw=2.4, near_end=True)
    assert (web_crippling(N=20, **kw).capacity
            > web_crippling(N=10, **kw).capacity)
