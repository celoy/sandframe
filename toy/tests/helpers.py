"""Shared helpers for the toy tests (initial conditions and interface measurement)."""

from __future__ import annotations

import dataclasses

import numpy as np

from sandframe_toy.params import Params, load_params
from sandframe_toy.state import cell_centres, fraction_above


def small_params(nx: int, ny: int, **kw) -> Params:
    """Default parameters on a smaller frame with the same 1 mm cells."""
    p = load_params()
    return dataclasses.replace(p, nx=nx, ny=ny, W_m=nx * p.dx_m, H_m=ny * p.dy_m, **kw)


def interface_height(total: np.ndarray, p: Params, level: float) -> np.ndarray:
    """Per column, lowest y where total phi crosses `level` upwards (linear interpolation)."""
    _, y = cell_centres(p)
    j = np.argmax(total >= level, axis=0)
    k = np.arange(p.nx)
    p0, p1 = total[j - 1, k], total[j, k]
    return y[j - 1] + (level - p0) / (p1 - p0) * p.dy_m


def cosine_amplitude(eta: np.ndarray, p: Params, k: float) -> float:
    x, _ = cell_centres(p)
    return float(abs(2.0 / p.nx * np.sum((eta - eta.mean()) * np.cos(k * x))))


def heavy_top_half(p: Params, cls: int, eps_m: float, k: float) -> np.ndarray:
    """Class `cls` at phi_bed above y = H/2 + eps cos(k x), clear liquid below."""
    phi = np.zeros((p.n_classes, p.ny, p.nx))
    phi[cls] = p.phi_bed * fraction_above(p, lambda x: p.H_m / 2 + eps_m * np.cos(k * x))
    return phi
