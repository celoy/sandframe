"""Transport of the sand classes: advection by the mixture, dispersion, hindered settling.

Each class obeys  d phi_i/dt + div(phi_i (u + w_i)) = div(D grad phi_i),
with w_i = V_i (1 - phi)^n along gravity (assumption 4) and isotropic dispersion D
(assumption 5). Finite volumes on cell centres, conservative face fluxes, zero flux through the
walls. One time step is split (Lie splitting, first order in time) into three sub-steps, each of
which keeps 0 <= phi_i and phi <= phi_bed on its own:

1. advect(): mixture velocity u (discretely divergence-free, from flow.py). The TOTAL phi is
   reconstructed at faces with a minmod-limited MUSCL slope and upwinded; each class then takes
   the share phi_i/phi of the upwind (donor) cell. Because minmod commutes with phi -> C - phi
   and div u = 0, the scheme maps phi_bed - phi the same way it maps phi, so positivity of both
   (for (|u|/dx + |v|/dy) dt <= 1/2) gives 0 <= phi <= phi_bed. Class shares are first-order
   upwind, so composition contrasts carried by the flow smear at first order (limitation).
2. diffuse(): explicit centred diffusion of each class; with D dt (1/dx^2 + 1/dy^2) <= 1/2 the
   update is a convex combination, so bounds hold for every class and for the total.
3. settle(): vertical settling (gravity along -y; the shift of Phase 2 will add the x part).
   The total settling flux through a face is the Godunov (sending/receiving) flux of the
   Kynch flux f(phi) = Vbar phi (1 - phi)^n, which is unimodal with maximum at
   phi* = 1/(n+1):  F = min(S(phi_above), R(phi_below)),  S = f(min(phi, phi*)),
   R = f(max(phi, phi*)), with Vbar = sum phi_i V_i / phi of the donor cell above.
   Inflow into nearly full cells is limited by the capacity  F <= (phi_bed - phi_below) dy/dt,
   so the receiving cell cannot exceed phi_bed. The total is shared among classes in proportion
   to phi_i V_i in the donor, so coarse grains leave a cell before fines.
   Consequence: a layer at exactly phi_bed is jammed (no internal settling), consistent with
   assumption 6 (no compaction), but its lower edge rains grains into clear liquid at the
   maximum Kynch flux f(phi*).
"""

from __future__ import annotations

import numpy as np


def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where(a * b > 0.0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def _shares(phi: np.ndarray, total: np.ndarray) -> np.ndarray:
    return np.divide(phi, total, out=np.zeros_like(phi), where=total > 0.0)


def advect(phi: np.ndarray, u: np.ndarray, v: np.ndarray, dx: float, dy: float, dt: float):
    """One explicit MUSCL-minmod step of phi (nc, ny, nx) by face velocities u, v [m/s]."""
    nc, ny, nx = phi.shape
    tot = phi.sum(axis=0)
    share = _shares(phi, tot)

    sx = np.zeros((ny, nx))
    sx[:, 1:-1] = _minmod(tot[:, 1:-1] - tot[:, :-2], tot[:, 2:] - tot[:, 1:-1])
    uf = u[:, 1:-1]
    fx_tot = np.where(uf > 0.0, uf * (tot[:, :-1] + 0.5 * sx[:, :-1]),
                      uf * (tot[:, 1:] - 0.5 * sx[:, 1:]))  # fmt: skip
    fx = np.zeros((nc, ny, nx + 1))
    fx[:, :, 1:-1] = fx_tot * np.where(uf > 0.0, share[:, :, :-1], share[:, :, 1:])

    sy = np.zeros((ny, nx))
    sy[1:-1, :] = _minmod(tot[1:-1, :] - tot[:-2, :], tot[2:, :] - tot[1:-1, :])
    vf = v[1:-1, :]
    fy_tot = np.where(vf > 0.0, vf * (tot[:-1, :] + 0.5 * sy[:-1, :]),
                      vf * (tot[1:, :] - 0.5 * sy[1:, :]))  # fmt: skip
    fy = np.zeros((nc, ny + 1, nx))
    fy[:, 1:-1, :] = fy_tot * np.where(vf > 0.0, share[:, :-1, :], share[:, 1:, :])

    return phi - dt * ((fx[:, :, 1:] - fx[:, :, :-1]) / dx + (fy[:, 1:, :] - fy[:, :-1, :]) / dy)


def diffuse(phi: np.ndarray, D: float, dx: float, dy: float, dt: float) -> np.ndarray:
    """One explicit centred diffusion step of every class, zero flux at the walls."""
    nc, ny, nx = phi.shape
    fx = np.zeros((nc, ny, nx + 1))
    fx[:, :, 1:-1] = -D * (phi[:, :, 1:] - phi[:, :, :-1]) / dx
    fy = np.zeros((nc, ny + 1, nx))
    fy[:, 1:-1, :] = -D * (phi[:, 1:, :] - phi[:, :-1, :]) / dy
    return phi - dt * ((fx[:, :, 1:] - fx[:, :, :-1]) / dx + (fy[:, 1:, :] - fy[:, :-1, :]) / dy)


def settling_fluxes(
    phi: np.ndarray, V: np.ndarray, n_rz: float, phi_bed: float, dy: float, dt: float
) -> np.ndarray:
    """Downward settling flux of each class through the horizontal faces, (nc, ny+1, nx) [m/s].

    Face j separates cell row j-1 (below, receiver) from row j (above, donor). Walls carry 0.
    """
    nc, ny, nx = phi.shape
    tot = phi.sum(axis=0)
    q = np.tensordot(V, phi, axes=1)  # sum_i phi_i V_i
    tot_up, tot_lo, q_up = tot[1:, :], tot[:-1, :], q[1:, :]
    vbar = np.divide(q_up, tot_up, out=np.zeros_like(q_up), where=tot_up > 0.0)
    phi_star = 1.0 / (n_rz + 1.0)

    def kynch(p: np.ndarray) -> np.ndarray:
        return p * np.clip(1.0 - p, 0.0, None) ** n_rz

    send = vbar * kynch(np.minimum(tot_up, phi_star))
    receive = vbar * kynch(np.maximum(tot_lo, phi_star))
    capacity = np.clip(phi_bed - tot_lo, 0.0, None) * dy / dt
    f_tot = np.minimum(np.minimum(send, receive), capacity)

    weights = V[:, None, None] * phi[:, 1:, :]
    f = np.zeros((nc, ny + 1, nx))
    f[:, 1:-1, :] = f_tot * np.divide(weights, q_up, out=np.zeros_like(weights), where=q_up > 0.0)
    return f


def settle(
    phi: np.ndarray, V: np.ndarray, n_rz: float, phi_bed: float, dy: float, dt: float
) -> np.ndarray:
    """One settling step (gravity along -y). V: Stokes velocities (nc,) [m/s]."""
    f = settling_fluxes(phi, V, n_rz, phi_bed, dy, dt)
    return phi + dt * (f[:, 1:, :] - f[:, :-1, :]) / dy
