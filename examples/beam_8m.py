"""Design an 8 m simply supported floor beam to ECP 205-LRFD.

Runs the full sequence: loads -> combinations -> trial section -> classify ->
every limit state -> serviceability. Shows why web crippling at a short end
bearing is the check that catches people out.

    python examples/beam_8m.py
"""

from ecp205 import (
    classify,
    deflection_check,
    flexure_I,
    get,
    governing,
    search,
    shear_web,
    steel,
    udl_deflection,
    web_crippling,
    web_local_yielding,
)

# ---------------------------------------------------------------- the problem
SPAN = 800.0          # cm
TRIB = 3.0            # m - tributary width
D_AREA = 0.500        # t/m2 - 120 slab + finishes
L_AREA = 0.300        # t/m2 - office live load, from ECP 201
BEARING = 10.0        # cm  - length of bearing at the support
GRADE = "St 37"

Fy, Fu = steel(GRADE, t_mm=15)
wD, wL = D_AREA * TRIB / 100, L_AREA * TRIB / 100     # t/cm

print(f"8 m floor beam, {GRADE}, {TRIB} m tributary")
print(f"  service:  wD = {wD * 100:.2f} t/m   wL = {wL * 100:.2f} t/m")

# ---------------------------------------------------------------- cl. 1.4.1
combo, wu = governing(D=wD, L=wL)
Mu, Vu = wu * SPAN ** 2 / 8, wu * SPAN / 2
print(f"  governing: {combo}  ->  {wu * 100:.2f} t/m")
print(f"  Mu = {Mu:.0f} t.cm ({Mu / 100:.2f} t.m)    Vu = {Vu:.2f} t\n")

# ---------------------------------------------------------------- trial
Zx_req = Mu / (0.85 * Fy)
print(f"first pass, laterally restrained: Zx >= {Zx_req:.0f} cm3")
for row in search(family="IPE", Zx_min=Zx_req)[:3]:
    print(f"   {row['name']:<9s} {row['mass']:5.1f} kg/m  Zx={row['Zx']:.0f}")

# ---------------------------------------------------------------- check
for name in ("IPE 400", "IPE 450"):
    sec = get(name)
    sw = sec.mass / 1000 / 100                        # t/cm
    wu2, ws2 = wu + 1.2 * sw, wD + wL + sw
    Mu2, Vu2 = wu2 * SPAN ** 2 / 8, wu2 * SPAN / 2

    cl = classify(sec, Fy)
    print(f"\n{'=' * 70}\n{sec}\n  class {cl['class']}  "
          f"(flange {cl['flange']['lambda']:.1f}<={cl['flange']['lp']:.1f}, "
          f"web {cl['web']['lambda']:.1f}<={cl['web']['lp']:.1f})")
    print(f"  with self weight: Mu={Mu2:.0f} t.cm  Vu={Vu2:.2f} t")

    checks = [
        ("flexure, braced by slab", flexure_I(sec, Fy, Lb=1.0, Mu=Mu2)),
        ("flexure, braced 1/3 pts", flexure_I(sec, Fy, Lb=SPAN / 3,
                                              Cb=1.14, Mu=Mu2)),
        ("flexure, unbraced", flexure_I(sec, Fy, Lb=SPAN, Cb=1.14, Mu=Mu2)),
        ("shear", shear_web(sec, Fy, Vu=Vu2)),
        ("web local yielding", web_local_yielding(
            k=sec.tf + sec.r, N=BEARING, Fyw=Fy, tw=sec.tw,
            Ru=Vu2, near_end=True)),
        ("web crippling", web_crippling(
            tw=sec.tw, tf=sec.tf, d=sec.d, N=BEARING, Fyw=Fy,
            Ru=Vu2, near_end=True)),
    ]
    for label, chk in checks:
        print(f"   {label:<26s} ratio {chk.ratio:5.2f}  {chk.status:<4s} "
              f"[{chk.clause}] {chk.governing}")

    dmax = udl_deflection(ws2, SPAN, sec.Ix)
    d2 = udl_deflection(wL, SPAN, sec.Ix)
    for row in ("other_beams", "brittle_finish"):
        cm = deflection_check(dmax, SPAN, row, "delta_max")
        c2 = deflection_check(d2, SPAN, row, "delta2_max")
        print(f"   deflection {row:<15s} "
              f"d2={d2:.2f}/{c2.capacity:.2f} r={c2.ratio:.2f} {c2.status:<4s} | "
              f"dmax={dmax:.2f}/{cm.capacity:.2f} r={cm.ratio:.2f} {cm.status}")

# ---------------------------------------------------------------- the fix
print(f"\n{'=' * 70}\nweb crippling vs bearing length, IPE 400 (cl. 10.4):")
sec = get("IPE 400")
Vu2 = (wu + 1.2 * sec.mass / 1000 / 100) * SPAN / 2
for N in (10, 15, 20, 25):
    c = web_crippling(tw=sec.tw, tf=sec.tf, d=sec.d, N=N, Fyw=Fy,
                      Ru=Vu2, near_end=True)
    print(f"   N = {N:2d} cm -> phiRn = {c.capacity:5.2f} t   "
          f"ratio {c.ratio:.2f}  {c.status}")

print("\nConclusion: IPE 450 with the compression flange braced by the slab and")
print("a bearing of at least 15 cm. Unbraced over the full 8 m no IPE works -")
print("an HEA 300 would be needed, 88 kg/m against 78.")
