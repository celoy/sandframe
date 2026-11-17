"""Parameters of the sand-frame toy: TOML loader, frozen dataclasses, derived quantities.

Modelling assumptions (repeated from CLAUDE.md):

1. Gap-averaged Darcy flow (Hele-Shaw), permeability K = b^2/12. Valid for structures much
   larger than b and negligible inertia.
2. The suspension is a Newtonian mixture of constant viscosity mu_m = mu_r * mu_f. The real
   viscosity diverges near packing and dense layers can have a yield stress; ignored.
3. Boussinesq: density differences enter only the buoyancy.
4. Relative to the mixture, grains settle along gravity at the hindered Stokes velocity
   V_i (1 - phi)^n, V_i = (rho_i - rho_f) g d_i^2 / (18 mu_f), n = n_rz (Richardson-Zaki).
5. Isotropic diffusion D regularises the finger scale. For a sharp interface in Darcy flow with
   equal viscosities, n(k) = K g drho k / (2 mu_m) - D k^2, maximal at lambda_star when
   D = lambda_star K g drho / (8 pi mu_m).
6. The deposit (bed) is rigid, impermeable, packed at phi_bed; no compaction.
7. The flip is instantaneous.

Run ``python -m sandframe_toy.params`` to print the derived quantities.
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_PATH = Path(__file__).resolve().parents[2] / "params" / "default.toml"


@dataclass(frozen=True)
class SandClass:
    """One grain class: diameter d_m [m], grain density rho_kg_m3 [kg/m^3], display colour."""

    name: str
    d_m: float
    rho_kg_m3: float
    color_rgb: tuple[int, int, int]


@dataclass(frozen=True)
class Params:
    """All toy parameters (SI). Derived quantities are properties."""

    W_m: float
    H_m: float
    b_m: float
    nx: int
    ny: int
    rho_f_kg_m3: float
    mu_f_pa_s: float
    g_m_s2: float
    mu_r: float
    phi_bed: float
    n_rz: float
    theta_r_deg: float
    lambda_star_over_b: float
    interface_noise: float
    cfl: float
    diffusion_number: float
    dt_max_s: float
    seed: int
    classes: tuple[SandClass, ...]

    # --- grid -------------------------------------------------------------------------------
    @property
    def dx_m(self) -> float:
        return self.W_m / self.nx

    @property
    def dy_m(self) -> float:
        return self.H_m / self.ny

    @property
    def n_classes(self) -> int:
        return len(self.classes)

    # --- flow -------------------------------------------------------------------------------
    @property
    def K_m2(self) -> float:
        """Hele-Shaw permeability b^2/12 [m^2]."""
        return self.b_m**2 / 12.0

    @property
    def mu_m_pa_s(self) -> float:
        """Effective mixture viscosity mu_r mu_f [Pa s]."""
        return self.mu_r * self.mu_f_pa_s

    @property
    def mobility_m3_s_kg(self) -> float:
        """K / mu_m [m^2 / (Pa s)]."""
        return self.K_m2 / self.mu_m_pa_s

    @property
    def drho_classes_kg_m3(self) -> np.ndarray:
        """rho_i - rho_f for each class [kg/m^3]."""
        return np.array([c.rho_kg_m3 - self.rho_f_kg_m3 for c in self.classes])

    @property
    def drho_ref_kg_m3(self) -> float:
        """Reference density excess: packed reservoir (phi_bed) of equal parts of each class."""
        return self.phi_bed * float(self.drho_classes_kg_m3.mean())

    @property
    def U_m_s(self) -> float:
        """Plume velocity scale U = K g drho_ref / mu_m [m/s]."""
        return self.mobility_m3_s_kg * self.g_m_s2 * self.drho_ref_kg_m3

    @property
    def lambda_star_m(self) -> float:
        return self.lambda_star_over_b * self.b_m

    @property
    def D_m2_s(self) -> float:
        """Dispersion D = lambda_star K g drho_ref / (8 pi mu_m) [m^2/s] (derived)."""
        return self.lambda_star_m * self.U_m_s / (8.0 * math.pi)

    def growth_rate_s(self, k_per_m: float | np.ndarray, drho_kg_m3: float | None = None):
        """Sharp-interface growth rate n(k) = K g drho k / (2 mu_m) - D k^2 [1/s]."""
        drho = self.drho_ref_kg_m3 if drho_kg_m3 is None else drho_kg_m3
        a = self.mobility_m3_s_kg * self.g_m_s2 * drho / 2.0
        return a * k_per_m - self.D_m2_s * k_per_m**2

    # --- settling ---------------------------------------------------------------------------
    @property
    def stokes_velocity_m_s(self) -> np.ndarray:
        """Isolated-grain Stokes velocity V_i [m/s] (positive = along gravity)."""
        d = np.array([c.d_m for c in self.classes])
        return self.drho_classes_kg_m3 * self.g_m_s2 * d**2 / (18.0 * self.mu_f_pa_s)

    @property
    def reynolds_particle(self) -> np.ndarray:
        """rho_f V_i d_i / mu_f [-]."""
        d = np.array([c.d_m for c in self.classes])
        return self.rho_f_kg_m3 * self.stokes_velocity_m_s * d / self.mu_f_pa_s

    @property
    def fall_time_s(self) -> np.ndarray:
        """H / V_i [s]: unhindered fall time over the frame height."""
        return self.H_m / self.stokes_velocity_m_s

    # --- time step bounds -------------------------------------------------------------------
    @property
    def dt_diffusion_s(self) -> float:
        return 2.0 * self.diffusion_number / (self.D_m2_s * (self.dx_m**-2 + self.dy_m**-2))

    @property
    def dt_settling_s(self) -> float:
        return self.cfl * self.dy_m / float(self.stokes_velocity_m_s.max())

    def dt_advection_s(self, speed_m_s: float) -> float:
        """Advective bound for a speed along a diagonal (worst case)."""
        return self.cfl / (speed_m_s * (1.0 / self.dx_m + 1.0 / self.dy_m))

    def summary(self) -> str:
        """Human-readable table of the derived quantities (printed at start-up)."""
        k_star = 2.0 * math.pi / self.lambda_star_m
        lines = [
            f"grid            {self.nx} x {self.ny}, dx = {self.dx_m * 1e3:.3g} mm, "
            f"dy = {self.dy_m * 1e3:.3g} mm",
            f"K = b^2/12      {self.K_m2:.3e} m^2",
            f"mu_m            {self.mu_m_pa_s:.3g} Pa s",
            f"drho_ref        {self.drho_ref_kg_m3:.1f} kg/m^3 (phi_bed, equal parts)",
            f"U = K g drho/mu {self.U_m_s * 1e3:.2f} mm/s",
            f"lambda_star     {self.lambda_star_m * 1e3:.2f} mm "
            f"({self.lambda_star_m / self.dx_m:.1f} dx)",
            f"D               {self.D_m2_s:.3e} m^2/s",
            f"n(lambda_star)  {self.growth_rate_s(k_star):.3f} 1/s "
            f"(e-folding {1.0 / self.growth_rate_s(k_star):.2f} s)",
        ]
        for c, v, t, re in zip(
            self.classes,
            self.stokes_velocity_m_s,
            self.fall_time_s,
            self.reynolds_particle,
            strict=True,
        ):
            lines.append(
                f"{c.name:<8s} V = {v * 1e3:.4f} mm/s, H/V = {t:7.1f} s ({t / 60:5.1f} min), "
                f"Re_p = {re:.2e}"
            )
        lines += [
            f"dt bounds       diffusion {self.dt_diffusion_s:.4f} s, settling "
            f"{self.dt_settling_s:.4f} s, advection at U {self.dt_advection_s(self.U_m_s):.4f} s,"
            f" cap {self.dt_max_s} s",
        ]
        return "\n".join(lines)


def load_params(path: str | Path | None = None) -> Params:
    """Load a TOML parameter file (default: toy/params/default.toml)."""
    with open(DEFAULT_PATH if path is None else path, "rb") as f:
        t = tomllib.load(f)
    classes = tuple(
        SandClass(c["name"], float(c["d_m"]), float(c["rho_kg_m3"]), tuple(c["color_rgb"]))
        for c in t["classes"]
    )
    return Params(
        **{k: float(v) for k, v in t["cell"].items()},
        nx=int(t["grid"]["nx"]),
        ny=int(t["grid"]["ny"]),
        **{k: float(v) for k, v in t["fluid"].items()},
        **{k: float(v) for k, v in t["mixture"].items()},
        **{k: float(v) for k, v in t["numerics"].items()},
        seed=int(t["seed"]),
        classes=classes,
    )


if __name__ == "__main__":
    print(load_params().summary())
