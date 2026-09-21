# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-09-21

First release.

### Added
- Materials and Table 1.1 steel grades (`steel`, `bolt`).
- Load combinations eq. 1.1-1.6 with wind and earthquake reversal, and the
  heavy-live-load factor of cl. 1.4.1 (`load_combinations`, `governing`).
- Section classification, Tables 2.12a-d (`classify`, `web_limits`,
  `flange_limits`, `angle_limit`, `tee_limit`, `tube_limits`).
- Tension, cl. 3.1 - gross yielding, net fracture, shear lag, block shear.
- Compression, cl. 4.2 - flexural buckling, Q factor, effective widths.
- Flexure, cl. 5.1.3 - compact and non-compact, all three Lb regimes, Cb from
  end moments or quarter points; web shear cl. 5.2.2.
- Beam-column interaction cl. 7.1, with B1/B2 second-order magnifiers cl. 2.2.2.
- Bolts, cl. 8.5-8.11 - shear, bearing, tension, combined, slip resistance,
  prying, end plate thickness, long joints, packing, base plate bearing.
- Welds, cl. 9.5-9.6 - fillet, groove, combined stress, minimum sizes.
- Concentrated forces, cl. 10.2-10.4 - flange local bending, web local
  yielding, web crippling.
- Serviceability, ch. 14 - deflection Table 14.1, drift Table 14.2.
- Catalogue of 53 European profiles (IPE, HEA, HEB) with a geometric
  self-consistency check.
- Test suite pinned to the code's own tables, including Table 8.3.
