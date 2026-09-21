"""Load factors and combinations - ECP 205-LRFD cl. 1.4.1, eq. 1.1 to 1.6.

The nominal loads themselves are **not** in this code. They come from the
Egyptian Code for Loads and Forces (ECP 201).
"""

from __future__ import annotations

__all__ = ["COMBINATIONS", "load_combinations", "governing"]

#: The six combinations of cl. 1.4.1, expanded for the "or" and the plus/minus
#: cases. Each entry maps a label to the factors on (D, L, Lr, W, EQ).
#: The factor on L shown as ``None`` is 0.5 normally and 1.0 for garages, places
#: of public assembly and any area with L > 500 kg/m2.
COMBINATIONS = {
    "1.1  1.4D":                  (1.4, 0.0, 0.0, 0.0, 0.0),
    "1.2  1.2D+1.6L+0.5Lr":       (1.2, 1.6, 0.5, 0.0, 0.0),
    "1.3a 1.2D+1.6Lr+0.5L":       (1.2, None, 1.6, 0.0, 0.0),
    "1.3b 1.2D+1.6Lr+0.8W":       (1.2, 0.0, 1.6, 0.8, 0.0),
    "1.3c 1.2D+1.6Lr-0.8W":       (1.2, 0.0, 1.6, -0.8, 0.0),
    "1.4a 1.2D+1.3W+0.5L+0.5Lr":  (1.2, None, 0.5, 1.3, 0.0),
    "1.4b 1.2D-1.3W+0.5L+0.5Lr":  (1.2, None, 0.5, -1.3, 0.0),
    "1.5a 1.2D+1.0EQ+0.5L":       (1.2, None, 0.0, 0.0, 1.0),
    "1.5b 1.2D-1.0EQ+0.5L":       (1.2, None, 0.0, 0.0, -1.0),
    "1.6a 0.9D+1.3W":             (0.9, 0.0, 0.0, 1.3, 0.0),
    "1.6b 0.9D-1.3W":             (0.9, 0.0, 0.0, -1.3, 0.0),
    "1.6c 0.9D+1.0EQ":            (0.9, 0.0, 0.0, 0.0, 1.0),
    "1.6d 0.9D-1.0EQ":            (0.9, 0.0, 0.0, 0.0, -1.0),
}


def load_combinations(D: float = 0.0, L: float = 0.0, Lr: float = 0.0,
                      W: float = 0.0, EQ: float = 0.0,
                      heavy_live: bool = False) -> dict[str, float]:
    """Factored value of every combination, for one load effect.

    Apply to a single effect at a time - an axial force, a moment, a UDL - with
    consistent units. Wind and earthquake are expanded both ways, so the result
    already contains the reversal cases.

    Args:
        D, L, Lr, W, EQ: the nominal (service) load effects.
        heavy_live: True for garages, places of public assembly and any area
            with L > 500 kg/m2, where the factor on L in combinations 1.3, 1.4
            and 1.5 becomes 1.0 instead of 0.5 (cl. 1.4.1).

    Returns:
        ``{label: factored effect}``.

    >>> c = load_combinations(D=1.5, L=0.9)
    >>> round(c["1.2  1.2D+1.6L+0.5Lr"], 3)
    3.24
    >>> heavy = load_combinations(D=1.5, L=0.9, heavy_live=True)
    >>> round(heavy["1.5a 1.2D+1.0EQ+0.5L"], 2)
    2.7
    """
    lf = 1.0 if heavy_live else 0.5
    out = {}
    for label, (fd, fl, flr, fw, feq) in COMBINATIONS.items():
        fl = lf if fl is None else fl
        out[label] = fd * D + fl * L + flr * Lr + fw * W + feq * EQ
    return out


def governing(D: float = 0.0, L: float = 0.0, Lr: float = 0.0, W: float = 0.0,
              EQ: float = 0.0, heavy_live: bool = False,
              minimum: bool = False) -> tuple[str, float]:
    """The governing combination and its value.

    Args:
        minimum: True to return the most negative value instead of the largest -
            what you want when checking uplift or stress reversal, where
            combination 1.6 usually governs.

    >>> governing(D=1.5, L=0.9)
    ('1.2  1.2D+1.6L+0.5Lr', 3.24...)
    """
    combos = load_combinations(D, L, Lr, W, EQ, heavy_live)
    pick = min if minimum else max
    label = pick(combos, key=combos.get)
    return label, combos[label]
