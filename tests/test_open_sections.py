"""Angles, double angles and tees (ecp205.open_sections).

Geometry is tested against independent derivations - closed-form rectangles,
the EN 10056-1 area formula, principal-axis invariants and a rotated re-analysis
- and the catalogue against its own tabulated witnesses. None of these compare
the implementation with its own output.
"""

import math

import pytest

from ecp205 import (
    ANGLES,
    EQUAL_ANGLES,
    UNEQUAL_ANGLES,
    Angle,
    DoubleAngle,
    ISection,
    Tee,
    double_angle,
    get,
    get_angle,
    search_angles,
    tee_from,
    verify,
)
from ecp205._polygon import arc, props
from ecp205.open_sections import verify_angles


def _rect(x0, y0, w, h):
    """(A, xc, yc, Ix_own, Iy_own) of a rectangle."""
    return w * h, x0 + w / 2, y0 + h / 2, w * h ** 3 / 12, h * w ** 3 / 12


def _combine(parts):
    """Closed-form composite of rectangles: A, xc, yc, Ix, Iy, Ixy."""
    A = sum(p[0] for p in parts)
    xc = sum(p[0] * p[1] for p in parts) / A
    yc = sum(p[0] * p[2] for p in parts) / A
    Ix = sum(p[3] + p[0] * (p[2] - yc) ** 2 for p in parts)
    Iy = sum(p[4] + p[0] * (p[1] - xc) ** 2 for p in parts)
    Ixy = sum(p[0] * (p[1] - xc) * (p[2] - yc) for p in parts)
    return A, xc, yc, Ix, Iy, Ixy


# ------------------------------------------------------------ polygon engine
def test_polygon_rectangle_exact():
    p = props([(1, 2), (5, 2), (5, 8), (1, 8)])
    assert p.A == pytest.approx(24)
    assert (p.xc, p.yc) == pytest.approx((3, 5))
    assert p.Ix == pytest.approx(4 * 6 ** 3 / 12)
    assert p.Iy == pytest.approx(6 * 4 ** 3 / 12)
    assert p.Ixy == pytest.approx(0, abs=1e-9)


def test_polygon_circle_converges():
    """A 128-chord circle is within 0.1 % of pi r^2 and pi r^4 / 4."""
    pts = arc(0, 0, 3, 0, 360)[:-1]
    p = props(pts)
    assert p.A == pytest.approx(math.pi * 9, rel=1e-3)
    assert p.Ix == pytest.approx(math.pi * 3 ** 4 / 4, rel=1e-3)


def test_polygon_rejects_clockwise():
    with pytest.raises(ValueError, match="counter-clockwise"):
        props([(0, 0), (0, 1), (1, 1), (1, 0)])


# ------------------------------------------------------------ single angles
@pytest.mark.parametrize("b,d,t", [(10, 10, 1), (10, 7.5, 0.8), (15, 9, 1.2)])
def test_sharp_angle_matches_closed_form(b, d, t):
    """With no radii the outline is two rectangles - an exact, independent check
    of A, the centroid, Ix, Iy and the (negative) product of inertia."""
    a = Angle("sharp", b=b, d=d, t=t)
    A, xc, yc, Ix, Iy, Ixy = _combine([_rect(0, 0, d, t), _rect(0, t, t, b - t)])
    assert a.A == pytest.approx(A)
    assert (a.ex, a.ey) == pytest.approx((xc, yc))
    assert (a.Ix, a.Iy, a.Ixy) == pytest.approx((Ix, Iy, Ixy))
    assert a.Ixy < 0


@pytest.mark.parametrize("name", list(ANGLES))
def test_angle_area_formula(name):
    """EN 10056-1 area: A = t(b + d - t) + (1 - pi/4)(r1^2 - 2 r2^2)."""
    a = get_angle(name)
    A = (a.t * (a.b + a.d - a.t)
         + (1 - math.pi / 4) * (a.r1 ** 2 - 2 * a.r2 ** 2))
    assert a.A == pytest.approx(A, rel=2e-4)


@pytest.mark.parametrize("name", ["L 100x100x10", "L 100x75x8", "L 200x100x15",
                                  "L 65x50x5"])
def test_principal_axes_by_rotation(name):
    """Rotate the outline by -alpha and analyse it again: the product of inertia
    must vanish and Ix, Iy must become Iu, Iv."""
    a = get_angle(name)
    c, s = math.cos(math.radians(-a.alpha)), math.sin(math.radians(-a.alpha))
    p = props([(x * c - y * s, x * s + y * c) for x, y in a.outline()])
    assert p.Ixy == pytest.approx(0, abs=1e-6 * a.Iu)
    assert p.Ix == pytest.approx(a.Iu)
    assert p.Iy == pytest.approx(a.Iv)


@pytest.mark.parametrize("name", list(ANGLES))
def test_principal_invariants(name):
    a = get_angle(name)
    assert a.Iu + a.Iv == pytest.approx(a.Ix + a.Iy)
    assert a.Iu * a.Iv == pytest.approx(a.Ix * a.Iy - a.Ixy ** 2)
    assert a.Iu >= a.Ix >= a.Iv and a.Iu >= a.Iy >= a.Iv


@pytest.mark.parametrize("name", list(EQUAL_ANGLES))
def test_equal_angle_symmetry(name):
    """The major axis of an equal angle is its axis of symmetry: alpha = 45 deg
    and the shear centre lies on it (v0 = 0)."""
    a = get_angle(name)
    assert a.alpha == pytest.approx(45)
    assert a.v0 == pytest.approx(0, abs=1e-9)
    assert a.ex == pytest.approx(a.ey)


def test_unequal_angle_orientation():
    """Long leg vertical: the centroid sits higher than it sits out, Ix > Iy, and
    the major axis u lies between x and 45 deg - nearer the short leg."""
    a = get_angle("L 150x75x10")
    assert a.ey > a.ex and a.Ix > a.Iy
    assert 0 < a.alpha < 45


def test_angle_shear_centre_and_polar_radius():
    a = get_angle("L 100x100x10")
    assert (a.x0, a.y0) == pytest.approx((0.5 - a.ex, 0.5 - a.ey))
    assert math.hypot(a.u0, a.v0) == pytest.approx(math.hypot(a.x0, a.y0))
    assert a.ro2 == pytest.approx(a.x0 ** 2 + a.y0 ** 2 + (a.Iu + a.Iv) / a.A)


def test_angle_torsion_constants():
    """Thin-walled: J = (b' + d') t^3 / 3 and Cw = t^3 (b'^3 + d'^3) / 36 on the
    mid-line leg lengths b' = b - t/2."""
    a = get_angle("L 100x100x10")
    assert a.J == pytest.approx(2 * 9.5 * 1 / 3)
    assert a.Cw == pytest.approx(2 * 9.5 ** 3 / 36)


def test_angle_rejects_bad_geometry():
    with pytest.raises(ValueError, match="0 < t"):
        Angle("bad", b=5, d=5, t=6)
    with pytest.raises(ValueError, match="radii"):
        Angle("bad", b=5, d=5, t=0.5, r1=5)


# ------------------------------------------------------------ catalogue
def test_catalogue_size():
    assert len(EQUAL_ANGLES) == 58
    assert len(UNEQUAL_ANGLES) == 42


def test_verify_covers_angles():
    """sections.verify() now includes the angles, and both are clean."""
    assert verify() == []
    assert verify_angles() == []


def test_catalogue_witnesses_within_1_5_percent():
    """The worst agreement in the shipped table is 1.3 %; hold it there so a
    new row cannot drift toward the 3 % production threshold unnoticed."""
    assert verify_angles(tol=0.015) == []


def test_verify_catches_a_typo(monkeypatch):
    row = list(EQUAL_ANGLES["L 100x100x10"])
    row[1] = 1.2                                   # t 10 mm -> 12 mm
    monkeypatch.setitem(EQUAL_ANGLES, "L 100x100x10", tuple(row))
    monkeypatch.setitem(ANGLES, "L 100x100x10", tuple(row))
    bad = verify_angles()
    assert any(line.startswith("L 100x100x10  A") for line in bad)


@pytest.mark.parametrize("alias", ["L100x100x10", "l 100X100X10", "100x100x10",
                                   "L 100*100*10"])
def test_get_angle_aliases(alias):
    assert get_angle(alias).name == "L 100x100x10"


def test_get_angle_unknown():
    with pytest.raises(KeyError, match="not in the angle catalogue"):
        get_angle("L 110x110x10")


def test_search_angles_lightest_first():
    rows = search_angles(A_min=10, rv_min=1.5)
    masses = [r["mass"] for r in rows]
    assert masses == sorted(masses)
    assert all(r["A"] >= 10 and r["rv"] >= 1.5 for r in rows)
    assert all(get_angle(r["name"]).equal
               for r in search_angles(A_min=10, equal=True))
    assert not any(get_angle(r["name"]).equal
                   for r in search_angles(A_min=10, equal=False))


# ------------------------------------------------------------ tees
def test_sharp_tee_matches_closed_form():
    t = Tee("T", d=20, bf=18, tf=1.35, tw=0.86)
    A, _, yc, Ix, Iy, _ = _combine([_rect(-9, 0, 18, 1.35),
                                    _rect(-0.43, 1.35, 0.86, 20 - 1.35)])
    assert t.A == pytest.approx(A)
    assert t.yc == pytest.approx(yc)
    assert (t.Ix, t.Iy) == pytest.approx((Ix, Iy))


@pytest.mark.parametrize("parent", ["IPE 400", "HEA 300", "HEB 200", "IPE 120"])
def test_tee_is_half_its_parent(parent):
    """Two tees make the parent: area equal to the ISection geometry exactly (both
    carry (4 - pi) r^2 of fillet), and Ix by the parallel axis theorem within
    the 1 % the ISection fillet approximation allows."""
    tab = get(parent)
    geo = ISection(parent, d=tab.d, bf=tab.bf, tf=tab.tf, tw=tab.tw, r=tab.r)
    t = tee_from(parent)
    assert 2 * t.A == pytest.approx(geo.A, rel=1e-4)
    Ix = 2 * (t.Ix + t.A * (tab.d / 2 - t.yc) ** 2)
    assert Ix == pytest.approx(geo.Ix, rel=0.01)
    assert 2 * t.Iy == pytest.approx(geo.Iy, rel=0.02)
    assert t.J == pytest.approx(tab.J / 2)
    assert t.rolled


def test_tee_shear_centre():
    """Shear centre at flange mid-thickness, below the centroid."""
    t = tee_from("IPE 400")
    assert t.y0 == pytest.approx(1.35 / 2 - t.yc)
    assert t.y0 < 0 and t.x0 == 0
    assert t.ro2 == pytest.approx(t.y0 ** 2 + (t.Ix + t.Iy) / t.A)


def test_welded_tee_default_torsion():
    t = Tee("WT", d=30, bf=20, tf=2, tw=1, rolled=False)
    assert t.J == pytest.approx((20 * 8 + 29 * 1) / 3)
    assert t.Cw == pytest.approx(8 * 8000 / 144 + 29 ** 3 / 36)


# ------------------------------------------------------------ double angles
def test_double_angle_parallel_axis():
    a = get_angle("L 80x80x8")
    p = double_angle("L 80x80x8", gap=1.2)
    assert p.A == pytest.approx(2 * a.A)
    assert p.Ix == pytest.approx(2 * a.Ix)
    assert p.Iy == pytest.approx(2 * (a.Iy + a.A * (0.6 + a.ex) ** 2))
    assert p.y0 == pytest.approx(0.4 - a.ey)
    assert p.Cw == 0.0 and p.J == pytest.approx(2 * a.J)


def test_double_angle_gap_only_raises_Iy():
    p0, p1 = double_angle("L 100x100x10"), double_angle("L 100x100x10", gap=1.5)
    assert p1.Ix == pytest.approx(p0.Ix)
    assert p1.Iy > p0.Iy


def test_double_angle_long_or_short_legs():
    ll = double_angle("L 150x90x10", gap=1.0, legs="long")
    sl = double_angle("L 150x90x10", gap=1.0, legs="short")
    assert "LLBB" in ll.name and "SLBB" in sl.name
    assert ll.A == pytest.approx(sl.A)
    assert ll.Ix > sl.Ix          # long legs vertical: deeper section
    assert sl.Iy > ll.Iy          # long legs outstanding: wider section


def test_double_angle_validation():
    with pytest.raises(ValueError, match="legs"):
        double_angle("L 100x75x8", legs="sideways")
    with pytest.raises(ValueError, match="gap"):
        DoubleAngle(get_angle("L 100x75x8"), gap=-1)
