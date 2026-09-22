# ecp205-lrfd

Steel member design to the **Egyptian Code of Practice for Steel Construction, Load
and Resistance Factor Design** — ECP 205-LRFD, 1st edition 2008, Ministerial Decree
359/2007, Housing and Building National Research Center.

Pure Python, no dependencies. Every function returns the demand, the design capacity,
the utilisation, **the clause it came from** and the governing limit state — so the
output can go straight into a calculation sheet a reviewer can follow.

```python
from ecp205 import steel, get, classify, flexure_I, shear_web

Fy, Fu = steel("St 37", t_mm=13.5)        # Table 1.1
beam = get("IPE 400")

print(flexure_I(beam, Fy, Lb=250, Cb=1.14, Mu=2400))
```

```
Flexure Mx                   [5.1.3]  Ru=2400  Rd=2666  ratio=0.900  OK
    governing: full plastic moment, eq. 5.3 (Lb <= Lp)
    class=1  Mp=3137  My=2776  Mr=2082  Lp=203.9  Lr=619.7  Lb=250  Cb=1.14  Mn=3137
```

## Why this exists

ECP 205 follows the AISC-LRFD structure closely, which makes it easy to reach for an
AISC library and assume the numbers carry over. They do not:

| | ECP 205-LRFD | AISC 360 |
|---|---|---|
| φ compression | **0.80** | 0.90 |
| φ flexure | **0.85** | 0.90 |
| φ web local yielding | **0.95** | 1.00 |
| wind load factor | **1.3** | 1.6 |
| Lp | **80 ry/√Fy** | 1.76 ry√(E/Fy) |
| shear area | **d × tw** (overall depth) | d × tw |
| units | **ton, cm, t/cm²** | kip, in, ksi |

Every constant in ECP 205 is dimensional. `127/√Fy`, `80ry/√Fy`, `1380Af/(d·Lb)`,
`20700`, `9500`, `6.25tf²` are only valid in **ton, centimetre, t/cm²**. This package
stays in those units throughout and converts only at the edges:

```python
from ecp205 import t_to_kN, tcm2_to_MPa
t_to_kN(15.44)      # 151.4 kN
tcm2_to_MPa(2.4)    # 235.4 MPa
```

## Install

```bash
pip install ecp205-lrfd
```

or from source:

```bash
git clone https://github.com/ahmedalisteel-star/ecp205-lrfd
cd ecp205-lrfd
pip install -e ".[dev]"
pytest
```

Python 3.9+. No runtime dependencies.

## What is covered

| Chapter | Provisions | Module |
|---|---|---|
| 1 | materials (Table 1.1), load combinations eq. 1.1–1.6, resistance factors | `materials`, `loads`, `core` |
| 2 | section classification Tables 2.12a–d, shear lag, B1/B2 second-order, slenderness limits | `classification`, `members` |
| 3 | tension — gross yielding, net fracture, block shear | `members` |
| 4 | compression — flexural buckling, Q for slender elements, effective widths | `members` |
| 5 | flexure — compact and non-compact, all three Lb regimes, Cb, web shear | `members` |
| 7 | beam-column interaction eq. 7.1a/7.1b | `members` |
| 8 | bolts — shear, bearing, tension, combined, slip, prying, end plates, base plates | `connections` |
| 9 | welds — fillet, groove, combined stress, minimum sizes | `connections` |
| 10 | concentrated forces — flange bending, web yielding, web crippling | `connections` |
| 14 | serviceability — deflection Table 14.1, drift Table 14.2 | `serviceability` |

**Not covered** — use the code directly: plate girders (ch. 6), flexural-torsional
buckling of tees and angles (cl. 4.3), laced and battened built-up members (cl. 4.4),
fatigue (ch. 11), composite construction (ch. 12), cold-formed sections (ch. 13).

## Worked example — an 8 m floor beam

```python
from ecp205 import (steel, get, search, classify, flexure_I, shear_web,
                    web_crippling, deflection_check, udl_deflection, governing)

L, Fy = 800.0, steel("St 37", 15)[0]          # cm, t/cm2

# service loads on the beam, t/cm
wD, wL = 0.015, 0.009
combo, wu = governing(D=wD, L=wL)             # cl. 1.4.1
Mu, Vu = wu * L**2 / 8, wu * L / 2
print(combo, round(Mu), round(Vu, 1))
# 1.2  1.2D+1.6L+0.5Lr 2592 13.0

# first pass, then check everything
beam = get(search(family="IPE", Zx_min=Mu / (0.85 * Fy))[0]["name"])
print(classify(beam, Fy)["class"])                        # 1
print(flexure_I(beam, Fy, Lb=1.0, Mu=Mu).ratio)           # 0.97 - braced by the slab
print(shear_web(beam, Fy, Vu=Vu).ratio)                   # 0.31
print(web_crippling(beam.tw, beam.tf, beam.d, N=10,
                    Fyw=Fy, Ru=Vu, near_end=True).ratio)  # 1.10 - FAILS
print(deflection_check(udl_deflection(wD + wL, L, beam.Ix), L,
                       "brittle_finish", "delta_max").ratio)
```

Web crippling at a short end bearing is the limit state that catches people out —
the end-zone coefficient in eq. 10.5 is half the interior one. `examples/beam_8m.py`
runs the full sequence.

## Validation

The test suite pins the implementation to values printed **in the code itself**, not
to its own output. The strongest of these is Table 8.3, which tabulates bolt
pretension, installation torque and slip resistance:

| Bolt (grade 10.9) | T, eq. 8.7 | Table 8.3 | Ps, eq. 8.8 | Table 8.3 |
|---|---|---|---|---|
| M16 | 9.89 t | 9.89 | 3.14 t | 3.16 |
| M20 | 15.44 t | 15.43 | 4.94 t | 4.93 |
| M24 | 22.23 t | 22.23 | 7.11 t | 7.11 |
| M30 | 35.34 t | 35.34 | 11.31 t | 11.30 |

Agreement here means the unit system, the material table and the bolt clauses are all
right together. The rest of the suite covers Table 1.1 cell by cell, Tables 2.12a–d,
Table 8.2, Table 9.6, Tables 14.1 and 14.2, continuity across every branch point
(λc = 1.1, Pu/φPn = 0.20, N/d = 0.2), and the catalogue's internal consistency.

```bash
pytest          # includes doctests
ruff check .
```

## Section catalogue

53 European hot-rolled profiles (IPE 100–600, HEA 100–600, HEB 100–600) with
`get`, `search` and `weight`.

```python
from ecp205 import search
search(family="IPE", Zx_min=1500)[0]
# {'name': 'IPE 450', 'mass': 77.6, 'd': 45.0, 'A': 98.8, ...}
```

`verify()` re-derives A, Ix and Zx from each profile's plate geometry, root fillets
included, and reports anything more than 3 % adrift — so a typo in the table cannot
pass silently. It returns `[]`.

**Catalogue values are nominal.** Confirm against the mill certificate or the
supplier's table before a section goes on a drawing.

### Angles, double angles and tees

101 EN 10056-1 angles (59 equal, 42 unequal), double angles back to back, and tees
cut from any I section — the sections cl. 4.3 needs. They carry the geometry
flexural-torsional buckling uses: principal axes (Iu, Iv, α), the shear centre
(x0, y0 and u0, v0), J, Cw and r̄o².

```python
from ecp205 import get_angle, double_angle, tee_from, search_angles

a = get_angle("L 100x100x10")     # a.rv = 1.95, a.Iu = 280, a.Iv = 73.0
p = double_angle("L 100x75x8", gap=1.0, legs="long")
t = tee_from("IPE 400")           # "1/2 IPE 400", y0 = -3.85
search_angles(A_min=12, rv_min=1.9)[0]["name"]
```

Every angle property is computed from the outline with the root and toe radii
included, so the set is self-consistent. The catalogue's own A, centroid, Ix and
Iy are kept as witnesses, and `verify()` checks them against the geometry the
same way it checks the I sections.

These are section properties only. The cl. 4.3 buckling check itself is not
implemented yet.

## Custom and built-up sections

```python
from ecp205 import ISection, flexure_I

pg = ISection("PG 900x300", d=90, bf=30, tf=2.5, tw=1.0, r=0.0, rolled=False)
print(pg)          # properties computed from the plate geometry
flexure_I(pg, 3.6, Lb=400)     # FL = 0.60 Fy for a welded section
```

Pass `rolled=False` for welded sections — it changes the Table 2.12c flange limits
and FL in cl. 5.1.3.1. Any property left as `None` is computed from the geometry,
fillets included.

## Scope and responsibility

This package implements the code. It is **not** the code, does not reproduce it, and
is not a substitute for owning a copy — you need the original for the clauses,
figures and tables it does not cover, and for anything you sign.

Loads are not in ECP 205. They come from the Egyptian Code for Loads and Forces
(**ECP 201**); this package only supplies the load *factors* of cl. 1.4.1. ECP 205
also excludes shells, tanks, silos, towers, suspension roofs and offshore structures
(cl. 1.2.1), and has no equivalent of AISC 341 for seismic detailing.

Structural design is a licensed activity. A competent engineer must own the result.
Provided as-is under the MIT licence, with no warranty and no liability — see
[LICENSE](LICENSE).

## Contributing

Issues and pull requests welcome, particularly for the uncovered chapters. Two rules:

1. **Cite the clause.** Every provision carries its clause and equation number in the
   docstring and in the `Check` it returns.
2. **Test against the code, not against the implementation.** If the code prints a
   table, test against the table.

## Licence

MIT — see [LICENSE](LICENSE).

The Egyptian Code of Practice for Steel Construction is published by the Housing and
Building National Research Center and is copyright. This repository contains no part
of that document — only an independent implementation of the design provisions, with
references to where each one is found.
