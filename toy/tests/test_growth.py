"""Toy test 7: linear growth of a sinusoidal perturbation of a flat heavy-over-light interface.

Two versions:
- with dispersion switched off (D ~ 0), the measured rate must match the sharp-interface Darcy
  rate A k, A = K g drho / (2 mu_m), within 30 % for resolved wavelengths (>= 10 dx). This
  checks the flow solver and the advection scheme.
- the brief's literal test with the default D against n(k) = A k - D k^2. It fails by a large
  factor (measured about 0.3 1/s with a maximum near 40 mm, against 1.7 1/s at 9 mm), because
  the interface diffuses to a thickness ~ 1/k_star within one growth time, so the sharp-
  interface formula does not apply. Kept as a strict xfail until the user decides on D.
"""

import dataclasses

import numpy as np
import pytest

from helpers import cosine_amplitude, heavy_top_half, interface_height
from sandframe_toy.params import load_params
from sandframe_toy.sim import Sim

FINE = 2  # fine class: settling negligible over the test


def _rate_until(p, lam: float, ka_max: float = 0.3, t_max: float = 12.0) -> float:
    """Fit of log amplitude over the last 2/3 of the run, stopped when k a > ka_max."""
    k = 2 * np.pi / lam
    sim = Sim(p, heavy_top_half(p, FINE, 5e-5, k))
    ts, la = [], []
    while sim.state.t_s < t_max:
        a = cosine_amplitude(interface_height(sim.state.total, p, p.phi_bed / 2), p, k)
        ts.append(sim.state.t_s)
        la.append(np.log(a))
        if k * a > ka_max:
            break
        sim.step_once()
    n = len(ts)
    return float(np.polyfit(ts[n // 3 :], la[n // 3 :], 1)[0])


@pytest.mark.parametrize("lam_mm", [10.0, 12.5, 20.0, 25.0, 50.0])
def test_growth_rate_without_dispersion(lam_mm: float) -> None:
    p = dataclasses.replace(load_params(), lambda_star_over_b=1e-4)  # D ~ 1e-10 m^2/s
    drho = p.phi_bed * p.drho_classes_kg_m3[FINE]
    k = 2 * np.pi / (lam_mm * 1e-3)
    expected = p.mobility_m3_s_kg * p.g_m_s2 * drho / 2 * k
    measured = _rate_until(p, lam_mm * 1e-3)
    assert abs(measured / expected - 1.0) < 0.3, (measured, expected)


@pytest.mark.xfail(strict=True, reason="sharp-interface n(k) invalid at default D; see docstring")
@pytest.mark.parametrize("lam_mm", [20.0, 25.0, 40.0])
def test_growth_rate_brief_formula_default_D(lam_mm: float) -> None:
    p = load_params()
    drho = p.phi_bed * p.drho_classes_kg_m3[FINE]
    k = 2 * np.pi / (lam_mm * 1e-3)
    sim = Sim(p, heavy_top_half(p, FINE, 5e-5, k))
    la = []
    for t in (1.0, 3.0):
        sim.run(t)
        la.append(
            np.log(cosine_amplitude(interface_height(sim.state.total, p, p.phi_bed / 2), p, k))
        )
    measured = (la[1] - la[0]) / 2.0
    expected = p.growth_rate_s(k, drho)
    assert abs(measured / expected - 1.0) < 0.3, (measured, expected)
