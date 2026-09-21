"""Section classification - ECP 205-LRFD cl. 2.3.1, Tables 2.12a-d.

Class 1 (compact) reaches the plastic moment, class 2 (non-compact) reaches the
yield moment, class 3 (slender) reaches neither and must be designed on effective
widths (cl. 2.3.1.3). If any element is class 3 the whole section is class 3.

All limits take Fy in **t/cm2**.
"""

from __future__ import annotations

import math

from .sections import ISection

__all__ = ["web_limits", "flange_limits", "angle_limit", "tee_limit",
           "tube_limits", "classify"]


def web_limits(Fy: float, alpha: float = 0.5,
               psi: float = -1.0) -> tuple[float, float]:
    """Table 2.12a - ``(lambda_p, lambda_r)`` for a web, ratio dw/tw.

    Args:
        Fy: yield strength, t/cm2.
        alpha: depth ratio of the compressed part, for the **compact** limit.
            0.5 = pure bending, 1.0 = pure compression.
        psi: stress ratio f2/f1, for the **non-compact** limit.
            -1 = pure bending, +1 = pure compression.

    >>> lp, lr = web_limits(2.4)
    >>> round(lp, 1), round(lr, 1)
    (82.0, 143.3)
    >>> round(web_limits(2.4, alpha=1.0, psi=1.0)[0], 1)
    37.4
    """
    s = math.sqrt(Fy)
    if abs(alpha - 0.5) < 1e-9:
        lp = 127 / s
    elif abs(alpha - 1.0) < 1e-9:
        lp = 58 / s
    elif alpha > 0.5:
        lp = (699 / s) / (13 * alpha - 1)
    else:
        lp = (63.6 / alpha) / s

    if abs(psi + 1) < 1e-9:
        lr = 222 / s
    elif abs(psi - 1) < 1e-9:
        lr = 64 / s
    elif psi > -1:
        lr = (222 / s) / (2 + psi)
    else:
        lr = 111 * (1 - psi) * math.sqrt(-psi) / s
    return lp, lr


def flange_limits(Fy: float, rolled: bool = True) -> tuple[float, float]:
    """Table 2.12c - ``(lambda_p, lambda_r)`` for an outstand flange in uniform
    compression due to Mx, ratio c/tf.

    >>> lp, lr = flange_limits(2.4, rolled=True)
    >>> round(lp, 2), round(lr, 2)
    (10.91, 21.3)
    >>> round(flange_limits(2.4, rolled=False)[0], 2)
    9.88
    """
    s = math.sqrt(Fy)
    return (16.9 if rolled else 15.3) / s, (33 if rolled else 28) / s


def internal_flange_limits(Fy: float) -> tuple[float, float]:
    """Table 2.12b - internal flange elements (box flanges), ratio b/tf."""
    s = math.sqrt(Fy)
    return 58 / s, 64 / s


def angle_limit(Fy: float, unequal: bool = False) -> float:
    """Table 2.12d - non-compact limit for an angle in compression.

    ``b/t <= 23/sqrt(Fy)``, or ``(b+h)/2t <= 17/sqrt(Fy)`` for unequal angles.
    Does not apply to angles in continuous contact with other components.
    """
    return (17 if unequal else 23) / math.sqrt(Fy)


def tee_limit(Fy: float) -> float:
    """Table 2.12d - non-compact limit for a T-section stem, ``b/t <= 30/sqrt(Fy)``."""
    return 30 / math.sqrt(Fy)


def tube_limits(Fy: float) -> tuple[float, float]:
    """Table 2.12d - ``(lambda_p, lambda_r)`` for a circular tube, ratio D/t.

    Note these are divided by **Fy**, not sqrt(Fy).

    >>> lp, lr = tube_limits(2.4)
    >>> round(lp, 1), round(lr, 1)
    (68.8, 87.9)
    """
    return 165 / Fy, 211 / Fy


def classify(sec: ISection, Fy: float, action: str = "bending") -> dict:
    """Classify a doubly symmetric I section (cl. 2.3.1, Table 2.12).

    Args:
        sec: the section.
        Fy: yield strength, t/cm2.
        action: ``"bending"`` (alpha = 0.5, psi = -1) or ``"compression"``
            (alpha = 1.0, psi = +1).

    Returns:
        A dict with a ``web`` and ``flange`` entry - each carrying the actual
        ratio, its two limits and its class - plus the governing ``class``.

    >>> from .sections import get
    >>> classify(get("IPE 400"), 2.4)["class"]
    1
    """
    if action not in ("bending", "compression"):
        raise ValueError("action must be 'bending' or 'compression'")
    alpha, psi = (0.5, -1.0) if action == "bending" else (1.0, 1.0)
    lpw, lrw = web_limits(Fy, alpha, psi)
    lpf, lrf = flange_limits(Fy, sec.rolled)
    lam_w, lam_f = sec.dw / sec.tw, sec.c / sec.tf

    def cls(lam: float, lp: float, lr: float) -> int:
        return 1 if lam <= lp else (2 if lam <= lr else 3)

    cw, cf = cls(lam_w, lpw, lrw), cls(lam_f, lpf, lrf)
    return {
        "web": {"lambda": lam_w, "lp": lpw, "lr": lrw, "class": cw},
        "flange": {"lambda": lam_f, "lp": lpf, "lr": lrf, "class": cf},
        "class": max(cw, cf),
        "action": action,
        "clause": "2.3.1 / Table 2.12",
    }
