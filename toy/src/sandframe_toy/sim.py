"""Time stepping of the toy: flow solve, then advection, dispersion and settling sub-steps.

Phase 1 step 1: suspension only (no bed, no flip yet). Gravity is fixed along -y.
"""

from __future__ import annotations

import numpy as np

from sandframe_toy import transport
from sandframe_toy.flow import StreamFunctionSolver, buoyancy_source, face_velocities
from sandframe_toy.params import Params
from sandframe_toy.state import SuspensionState


class Sim:
    """Owns the state and the flow solver. The latest psi, u, v are kept for diagnostics."""

    def __init__(self, params: Params, phi0: np.ndarray) -> None:
        p = params
        if phi0.shape != (p.n_classes, p.ny, p.nx):
            raise ValueError(f"phi0 shape {phi0.shape} != {(p.n_classes, p.ny, p.nx)}")
        self.p = p
        self.state = SuspensionState(phi0.astype(float).copy(), 0.0)
        self.g_vec = (0.0, -p.g_m_s2)
        self.V = p.stokes_velocity_m_s
        self.solver = StreamFunctionSolver(p.nx, p.ny, p.dx_m, p.dy_m)
        self.update_flow()

    def update_flow(self) -> None:
        """Solve for psi and the face velocities from the current state."""
        p = self.p
        rhs = buoyancy_source(
            self.state.density_anomaly(p), p.mobility_m3_s_kg, self.g_vec, p.dx_m, p.dy_m
        )
        self.psi = self.solver.solve(rhs)
        self.u, self.v = face_velocities(self.psi, p.dx_m, p.dy_m)

    def stable_dt(self) -> float:
        """Largest dt satisfying the advective, diffusive and settling bounds and dt_max."""
        p = self.p
        courant_rate = (
            np.maximum(np.abs(self.u[:, 1:]), np.abs(self.u[:, :-1])) / p.dx_m
            + np.maximum(np.abs(self.v[1:, :]), np.abs(self.v[:-1, :])) / p.dy_m
        ).max()
        dt = min(p.dt_max_s, p.dt_diffusion_s, p.dt_settling_s)
        if courant_rate > 0.0:
            dt = min(dt, p.cfl / courant_rate)
        return dt

    def step_once(self, dt_max: float | None = None) -> float:
        """Advance by one stable step (at most dt_max); returns the dt used [s]."""
        p = self.p
        self.update_flow()
        dt = self.stable_dt() if dt_max is None else min(self.stable_dt(), dt_max)
        phi = transport.advect(self.state.phi, self.u, self.v, p.dx_m, p.dy_m, dt)
        phi = transport.diffuse(phi, p.D_m2_s, p.dx_m, p.dy_m, dt)
        phi = transport.settle(phi, self.V, p.n_rz, p.phi_bed, p.dy_m, dt)
        self.state.phi = phi
        self.state.t_s += dt
        return dt

    def step(self, dt: float) -> None:
        """Advance by exactly dt [s], sub-cycling with stable steps."""
        t_end = self.state.t_s + dt
        while self.state.t_s < t_end - 1e-12:
            self.step_once(t_end - self.state.t_s)

    def run(self, t_end_s: float) -> None:
        self.step(t_end_s - self.state.t_s)
