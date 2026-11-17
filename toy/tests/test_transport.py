"""Toy tests 3, 4 and 8: mass conservation, bounds, hindered settling of a uniform suspension."""

import numpy as np

from helpers import small_params
from sandframe_toy.sim import Sim
from sandframe_toy.state import fraction_above


def _layered_reservoir_sim() -> Sim:
    """Three bands at phi_bed in the upper half (coarse lowest), rough lower interface."""
    p = small_params(60, 30)
    rng = np.random.default_rng(p.seed)
    noise = rng.uniform(-1.5e-3, 1.5e-3, p.nx)
    y0 = p.H_m / 2
    lower = fraction_above(p, lambda x: y0 + noise)
    b1 = fraction_above(p, lambda x: np.full_like(x, y0 + 5e-3))
    b2 = fraction_above(p, lambda x: np.full_like(x, y0 + 10e-3))
    phi = np.zeros((p.n_classes, p.ny, p.nx))
    phi[0] = p.phi_bed * (lower - b1)
    phi[1] = p.phi_bed * (b1 - b2)
    phi[2] = p.phi_bed * b2
    return Sim(p, phi)


def test_mass_conserved_and_bounds_over_1000_steps() -> None:
    sim = _layered_reservoir_sim()
    p = sim.p
    m0 = sim.state.volume_per_class(p)
    lo, hi, vmax = 0.0, 0.0, 0.0
    for _ in range(1000):
        sim.step_once()
        vmax = max(vmax, np.abs(sim.v).max())
        lo = min(lo, sim.state.phi.min())
        hi = max(hi, sim.state.total.max() - p.phi_bed)
    rel = np.abs(sim.state.volume_per_class(p) / m0 - 1.0)
    assert np.all(rel < 1e-9), rel
    assert lo >= -1e-12, lo  # test 4: phi_i >= 0 (round-off tolerance)
    assert hi <= 1e-12, hi  # test 4: phi <= phi_bed
    assert vmax > 1e-3  # the run did develop plumes


def test_uniform_suspension_front_descends_at_hindered_velocity() -> None:
    p = small_params(8, 100)
    phi0 = 0.05
    phi = np.zeros((p.n_classes, p.ny, p.nx))
    phi[0] = phi0  # coarse only, no lateral gradient -> no flow
    sim = Sim(p, phi)
    y = (np.arange(p.ny) + 0.5) * p.dy_m
    times, fronts = [], []
    for t in np.arange(3.0, 15.01, 1.0):
        sim.run(t)
        col = sim.state.total[:, 0]
        j = np.nonzero(col >= phi0 / 2)[0].max()  # highest cell above half level
        fronts.append(y[j] + (col[j] - phi0 / 2) / (col[j] - col[j + 1]) * p.dy_m)
        times.append(t)
    speed = -np.polyfit(times, fronts, 1)[0]
    expected = p.stokes_velocity_m_s[0] * (1 - phi0) ** p.n_rz
    assert np.abs(sim.v).max() < 1e-12
    assert abs(speed / expected - 1.0) < 0.02, (speed, expected)
