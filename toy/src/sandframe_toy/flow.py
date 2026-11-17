"""Gap-averaged Darcy flow of the mixture, in streamfunction form.

Assumptions: Hele-Shaw (Darcy) flow with permeability K = b^2/12 and constant mixture viscosity
mu_m; Boussinesq; no inertia; the bed is impermeable (psi = 0 inside it). Valid only for flow
structures larger than the gap b.

Darcy, with gravity g_vec = (g_x, g_y) in frame coordinates:
    u = -(K/mu_m)(grad p - Delta_rho g_vec),   div u = 0,   u = d psi/dy, v = -d psi/dx,
so that
    lap psi = (K/mu_m) (g_x d(Delta_rho)/dy - g_y d(Delta_rho)/dx),   psi = 0 on walls and bed.

Discretisation (staggered, exactly divergence-free):
- psi lives on cell corners (nodes), array (ny+1, nx+1); node (j, k) is at (k dx, j dy).
- u on vertical faces, array (ny, nx+1): u[j, k] = (psi[j+1, k] - psi[j, k]) / dy.
- v on horizontal faces, array (ny+1, nx): v[j, k] = -(psi[j, k+1] - psi[j, k]) / dx.
- The source at a node uses the four surrounding cells. This is exactly the discrete curl of
  the face-wise Darcy law with face densities averaged from the two adjacent cells, so the
  pressure is eliminated exactly and the discrete divergence of (u, v) vanishes in every cell.
- A node is an unknown iff it is interior and its four cells are all fluid; other nodes are 0.
The 5-point Laplacian on the unknown nodes is factorised once with scipy's splu and reused
until the solid (bed) mask changes.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu


class StreamFunctionSolver:
    """Poisson solver for psi on the nodes of an nx x ny cell grid with an optional solid mask."""

    def __init__(self, nx: int, ny: int, dx_m: float, dy_m: float) -> None:
        self.nx, self.ny, self.dx, self.dy = nx, ny, dx_m, dy_m
        self._solid: np.ndarray | None = None
        self._assemble(np.zeros((ny, nx), dtype=bool))

    def _assemble(self, solid: np.ndarray) -> None:
        fluid = ~solid
        free = np.zeros((self.ny + 1, self.nx + 1), dtype=bool)
        free[1:-1, 1:-1] = fluid[:-1, :-1] & fluid[:-1, 1:] & fluid[1:, :-1] & fluid[1:, 1:]
        idx = -np.ones(free.shape, dtype=np.int64)
        n = int(free.sum())
        idx[free] = np.arange(n)
        jj, kk = np.nonzero(free)
        me = idx[jj, kk]
        cx, cy = 1.0 / self.dx**2, 1.0 / self.dy**2
        rows, cols, vals = [me], [me], [np.full(n, -2.0 * (cx + cy))]
        for dj, dk, c in ((0, 1, cx), (0, -1, cx), (1, 0, cy), (-1, 0, cy)):
            nb = idx[jj + dj, kk + dk]
            ok = nb >= 0
            rows.append(me[ok])
            cols.append(nb[ok])
            vals.append(np.full(int(ok.sum()), c))
        a = csc_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        )
        self._lu = splu(a) if n > 0 else None
        self._free = free
        self._solid = solid.copy()

    def solve(self, rhs_nodes: np.ndarray, solid: np.ndarray | None = None) -> np.ndarray:
        """Solve lap psi = rhs on free nodes, psi = 0 elsewhere. Refactorises if the mask changed.

        rhs_nodes: (ny+1, nx+1) [1/s]; solid: (ny, nx) bool, True in bed cells.
        """
        if solid is None:
            solid = np.zeros((self.ny, self.nx), dtype=bool)
        if not np.array_equal(solid, self._solid):
            self._assemble(solid)
        psi = np.zeros((self.ny + 1, self.nx + 1))
        if self._lu is not None:
            psi[self._free] = self._lu.solve(rhs_nodes[self._free])
        return psi


def buoyancy_source(
    drho: np.ndarray, mobility: float, g_vec: tuple[float, float], dx: float, dy: float
) -> np.ndarray:
    """Node source (K/mu_m)(g_x d(drho)/dy - g_y d(drho)/dx), shape (ny+1, nx+1) [1/s].

    drho: (ny, nx) density anomaly [kg/m^3], zero in bed cells; mobility = K/mu_m.
    Boundary nodes get 0 (they are fixed anyway).
    """
    gx, gy = g_vec
    c00, c01 = drho[:-1, :-1], drho[:-1, 1:]  # cells below the node: left, right
    c10, c11 = drho[1:, :-1], drho[1:, 1:]  # cells above the node: left, right
    ddx = (c01 + c11 - c00 - c10) / (2.0 * dx)
    ddy = (c10 + c11 - c00 - c01) / (2.0 * dy)
    rhs = np.zeros((drho.shape[0] + 1, drho.shape[1] + 1))
    rhs[1:-1, 1:-1] = mobility * (gx * ddy - gy * ddx)
    return rhs


def face_velocities(psi: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    """u (ny, nx+1) on vertical faces and v (ny+1, nx) on horizontal faces [m/s]."""
    u = (psi[1:, :] - psi[:-1, :]) / dy
    v = -(psi[:, 1:] - psi[:, :-1]) / dx
    return u, v


def cell_velocities(u: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Face velocities averaged to cell centres, each (ny, nx) [m/s]."""
    return 0.5 * (u[:, 1:] + u[:, :-1]), 0.5 * (v[1:, :] + v[:-1, :])
