"""
ecp205 - the Egyptian Code of Practice for Steel Construction, LRFD.

ECP 205-LRFD, 1st edition 2008, Ministerial Decree 359/2007, Housing and
Building National Research Center.

Units are the code's own: **ton, centimetre, t/cm2**, moments in t.cm. Every
tabulated constant in the code is dimensional and only valid there. Convert at
the boundaries of a calculation, never inside a clause::

    from ecp205 import steel, get, classify, flexure_I, shear_web

    Fy, Fu = steel("St 37", t_mm=13.5)     # Table 1.1
    beam = get("IPE 400")
    print(classify(beam, Fy))              # Table 2.12
    print(flexure_I(beam, Fy, Lb=250, Cb=1.14, Mu=2400))   # cl. 5.1.3
    print(shear_web(beam, Fy, Vu=28))                      # cl. 5.2.2

Every function returns a :class:`~ecp205.core.Check` carrying the demand, the
design capacity, the utilisation, the clause, the governing limit state and the
intermediate values.

This package implements the code. It is not a substitute for it, and it does not
reproduce it - a competent engineer must own the result.
"""

from .classification import (
    angle_limit,
    classify,
    flange_limits,
    internal_flange_limits,
    tee_limit,
    tube_limits,
    web_limits,
)
from .connections import (
    base_plate_bearing,
    bolt_bearing,
    bolt_combined,
    bolt_pretension,
    bolt_shear,
    bolt_tension,
    end_plate_exact,
    end_plate_thickness,
    fillet_weld,
    fillet_weld_combined,
    flange_local_bending,
    groove_weld,
    installation_torque,
    long_joint_factor,
    min_fillet_size,
    packing_factor,
    prying_force,
    slip_resistance,
    web_crippling,
    web_local_yielding,
)
from .core import (
    ALPHA_T,
    DENSITY,
    LAMBDA_MAX,
    NU,
    PHI,
    Check,
    E,
    G,
    t_to_kN,
    tcm2_to_MPa,
    tcm_to_kNm,
)
from .loads import COMBINATIONS, governing, load_combinations
from .materials import (
    BOLT_AREAS,
    BOLT_GRADES,
    FRICTION,
    SLIP_FS,
    STEEL_GRADES,
    bolt,
    steel,
)
from .members import (
    B1,
    B2_drift,
    B2_euler,
    Cb_ends,
    Cb_quarter,
    Cm_endmoments,
    Pe,
    Q_factor,
    beam_column,
    block_shear,
    compression,
    effective_width_stiffened,
    effective_width_unstiffened,
    flexure_I,
    shear_lag_U,
    shear_web,
    tension,
)
from .sections import (
    CATALOGUE,
    HEA,
    HEB,
    IPE,
    ISection,
    get,
    search,
    verify,
    weight,
)
from .serviceability import (
    DEFLECTION_LIMITS,
    DRIFT_LIMITS,
    deflection_check,
    drift_check,
    required_Ix,
    udl_deflection,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # core
    "Check", "E", "G", "NU", "DENSITY", "ALPHA_T", "PHI", "LAMBDA_MAX",
    "t_to_kN", "tcm2_to_MPa", "tcm_to_kNm",
    # materials
    "steel", "bolt", "STEEL_GRADES", "BOLT_GRADES", "BOLT_AREAS",
    "FRICTION", "SLIP_FS",
    # loads
    "load_combinations", "governing", "COMBINATIONS",
    # sections
    "ISection", "get", "search", "weight", "verify",
    "CATALOGUE", "IPE", "HEA", "HEB",
    # classification
    "classify", "web_limits", "flange_limits", "internal_flange_limits",
    "angle_limit", "tee_limit", "tube_limits",
    # members
    "tension", "shear_lag_U", "block_shear",
    "compression", "Q_factor",
    "effective_width_unstiffened", "effective_width_stiffened",
    "Cb_ends", "Cb_quarter", "flexure_I", "shear_web",
    "Pe", "Cm_endmoments", "B1", "B2_drift", "B2_euler", "beam_column",
    # connections
    "bolt_shear", "bolt_bearing", "bolt_tension", "bolt_combined",
    "bolt_pretension", "installation_torque", "slip_resistance",
    "prying_force", "end_plate_thickness", "end_plate_exact",
    "long_joint_factor", "packing_factor", "base_plate_bearing",
    "min_fillet_size", "fillet_weld", "fillet_weld_combined", "groove_weld",
    "flange_local_bending", "web_local_yielding", "web_crippling",
    # serviceability
    "deflection_check", "drift_check", "udl_deflection", "required_Ix",
    "DEFLECTION_LIMITS", "DRIFT_LIMITS",
]
