"""Derived quantities against the values quoted in CLAUDE.md."""

import math

import numpy as np

from sandframe_toy.params import load_params


def test_derived_values() -> None:
    p = load_params()
    assert math.isclose(p.K_m2, 7.5e-7)
    assert 0.009 < p.U_m_s < 0.011  # "about 1 cm/s"
    np.testing.assert_allclose(p.fall_time_s[:2], [27.2, 108.7], rtol=0.01)
    assert np.all(p.reynolds_particle < 0.1)


def test_D_puts_maximum_of_n_at_lambda_star() -> None:
    p = load_params()
    k = np.linspace(10.0, 3000.0, 200001)
    k_max = k[np.argmax(p.growth_rate_s(k))]
    assert math.isclose(2 * math.pi / k_max, p.lambda_star_m, rel_tol=1e-3)
