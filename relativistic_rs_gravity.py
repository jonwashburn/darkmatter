#!/usr/bin/env python3
"""
Minimal placeholder for prospective relativistic computations referenced in the paper.
Provides a function to compute an illustrative lensing enhancement using w(r).
"""
from __future__ import annotations
import numpy as np


def recognition_weight_relativistic(r_kpc: np.ndarray, w_nr: np.ndarray, phi_over_c2: np.ndarray | float = 0.0) -> np.ndarray:
    """Apply a simple relativistic correction factor ~ (1 + 2 Phi / c^2).
    r_kpc: radii [kpc]
    w_nr:  non-relativistic weight array
    phi_over_c2: dimensionless potential term (array or scalar)
    """
    corr = 1.0 + 2.0 * np.asarray(phi_over_c2)
    return w_nr * corr


def enhance_convergence(kappa_N: np.ndarray, w_rel: np.ndarray) -> np.ndarray:
    """Compute kappa_model = w_rel * kappa_N."""
    return w_rel * kappa_N


if __name__ == "__main__":
    r = np.logspace(-1, 2, 50)  # 0.1 .. 100 kpc
    w_nr = 1.0 + 0.5 * (1 - np.exp(-r / 3.0))
    phi_over_c2 = -1e-6
    w_rel = recognition_weight_relativistic(r, w_nr, phi_over_c2)
    kappa_N = 1.0 / (1.0 + (r / 35.0) ** 2)
    kappa_model = enhance_convergence(kappa_N, w_rel)
    print("demo:", float(kappa_model[10]))
