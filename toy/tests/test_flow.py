"""Toy tests 1 and 2: Poisson solver against a sine mode; sign of the buoyancy-driven flow."""

import numpy as np
import pytest

from sandframe_toy.flow import (
    StreamFunctionSolver,
    buoyancy_source,
    cell_velocities,
    face_velocities,
)


def _sine_error(n: int, with_bed: bool) -> float:
    """Max error for psi = sin(pi x/W) sin(2 pi (y - y0)/(H - y0)) on an n x n/2 grid."""
    W, H = 0.2, 0.1
    nx, ny = n, n // 2
    dx, dy = W / nx, H / ny
    solid = np.zeros((ny, nx), dtype=bool)
    jb = ny // 4 if with_bed else 0
    solid[:jb, :] = True  # flat bed: fluid domain is [0, W] x [jb dy, H]
    y0 = jb * dy
    a, b = np.pi / W, 2 * np.pi / (H - y0)
    X, Y = np.meshgrid(np.arange(nx + 1) * dx, np.arange(ny + 1) * dy)
    exact = np.sin(a * X) * np.sin(b * (Y - y0))
    exact[Y < y0] = 0.0
    psi = StreamFunctionSolver(nx, ny, dx, dy).solve(-(a**2 + b**2) * exact, solid)
    return float(np.abs(psi - exact).max())


@pytest.mark.parametrize("with_bed", [False, True])
def test_poisson_sine_mode_second_order(with_bed: bool) -> None:
    e1, e2 = _sine_error(40, with_bed), _sine_error(80, with_bed)
    assert e2 < 5e-3
    assert 3.5 < e1 / e2 < 4.5  # second order


def test_buoyancy_chain_sine_mode() -> None:
    """drho = C cos(a x) sin(b y) must give psi = sin(a x) sin(b y) through the full chain."""
    W, H, nx, ny, mob, g = 0.2, 0.1, 100, 50, 1.25e-6, 9.81
    dx, dy = W / nx, H / ny
    a, b = np.pi / W, np.pi / H
    xc, yc = (np.arange(nx) + 0.5) * dx, (np.arange(ny) + 0.5) * dy
    drho = (a**2 + b**2) / (mob * g * a) * np.cos(a * xc)[None, :] * np.sin(b * yc)[:, None]
    # lap psi = mob g d(drho)/dx = -(a^2 + b^2) sin(a x) sin(b y)
    psi = StreamFunctionSolver(nx, ny, dx, dy).solve(buoyancy_source(drho, mob, (0.0, -g), dx, dy))
    X, Y = np.meshgrid(np.arange(nx + 1) * dx, np.arange(ny + 1) * dy)
    assert np.abs(psi - np.sin(a * X) * np.sin(b * Y)).max() < 2e-3


def test_velocity_is_discretely_divergence_free() -> None:
    rng = np.random.default_rng(1)
    nx, ny, dx, dy = 30, 20, 1e-3, 1e-3
    psi = StreamFunctionSolver(nx, ny, dx, dy).solve(rng.standard_normal((ny + 1, nx + 1)))
    u, v = face_velocities(psi, dx, dy)
    div = (u[:, 1:] - u[:, :-1]) / dx + (v[1:, :] - v[:-1, :]) / dy
    assert np.abs(div).max() < 1e-9 * np.abs(u).max() / dx
    assert np.all(u[:, [0, -1]] == 0.0) and np.all(v[[0, -1], :] == 0.0)


def test_sign_heavy_right_half_sinks_clockwise() -> None:
    """Theta = 0: heavy right half sinks, light left half rises, clockwise cell (psi < 0)."""
    nx, ny, dx, dy = 40, 20, 5e-3, 5e-3
    drho = np.zeros((ny, nx))
    drho[:, nx // 2 :] = 100.0
    rhs = buoyancy_source(drho, 1.25e-6, (0.0, -9.81), dx, dy)
    psi = StreamFunctionSolver(nx, ny, dx, dy).solve(rhs)
    _, vc = cell_velocities(*face_velocities(psi, dx, dy))
    assert vc[ny // 2, 3 * nx // 4] < 0.0  # sinks on the right
    assert vc[ny // 2, nx // 4] > 0.0  # rises on the left
    assert psi.min() < 0.0 and psi.max() <= 0.0  # clockwise everywhere
