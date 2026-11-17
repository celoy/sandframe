"""Suspension state of the toy: volume fractions on cell centres.

Array convention: phi[i, j, k] is the volume fraction [-] of class i in the cell of row j
(y index, j = 0 at the bottom) and column k (x index, k = 0 at the left wall). Cell centres
sit at x = (k + 1/2) dx, y = (j + 1/2) dy. The bed (step 2) is not stored here.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from sandframe_toy.params import Params


@dataclass
class SuspensionState:
    """Volume fractions phi (n_classes, ny, nx) [-] and simulation time t_s [s]."""

    phi: np.ndarray
    t_s: float = 0.0

    def copy(self) -> SuspensionState:
        return SuspensionState(self.phi.copy(), self.t_s)

    @property
    def total(self) -> np.ndarray:
        """phi = sum_i phi_i, shape (ny, nx)."""
        return self.phi.sum(axis=0)

    def density_anomaly(self, p: Params) -> np.ndarray:
        """Delta rho = sum_i phi_i (rho_i - rho_f) [kg/m^3], shape (ny, nx)."""
        return np.tensordot(p.drho_classes_kg_m3, self.phi, axes=1)

    def volume_per_class(self, p: Params) -> np.ndarray:
        """Sand volume per unit gap width, sum over cells of phi_i dx dy [m^2]."""
        return self.phi.sum(axis=(1, 2)) * p.dx_m * p.dy_m


def cell_centres(p: Params) -> tuple[np.ndarray, np.ndarray]:
    """1-D arrays of cell-centre coordinates x (nx,) and y (ny,) [m]."""
    return (np.arange(p.nx) + 0.5) * p.dx_m, (np.arange(p.ny) + 0.5) * p.dy_m


def fraction_above(p: Params, y_interface: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Fraction (ny, nx) of each cell lying above the curve y = y_interface(x).

    The curve is evaluated at the cell-centre x, so the fraction is exact for curves that vary
    slowly on the scale dx. Used to build layered initial conditions with sub-cell interfaces.
    """
    x, _ = cell_centres(p)
    y_face_lo = np.arange(p.ny)[:, None] * p.dy_m
    eta = np.asarray(y_interface(x))[None, :]
    return np.clip((y_face_lo + p.dy_m - eta) / p.dy_m, 0.0, 1.0)
