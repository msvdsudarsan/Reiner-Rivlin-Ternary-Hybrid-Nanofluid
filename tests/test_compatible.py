import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from compatible_model import solve


def test_newtonian_limit():
    s = solve(0.0, 0.0, 0.0, "compatible", None, 1e-10, L=30.0)
    assert abs(s.y[1, 0] - 0.5102) < 5e-5
    assert abs(s.y[3, 0] + 0.6159) < 5e-5
    assert abs(s.y[4, -1] + 0.8845) < 5e-5


def test_closures_coincide_at_zero_K():
    a = solve(0.0, 0.5, 0.5, "compatible", None, 1e-9)
    b = solve(0.0, 0.5, 0.5, "zero-gradient", None, 1e-9)
    assert abs(a.y[1, 0] - b.y[1, 0]) < 1e-8


def test_closures_differ_for_positive_K():
    a = solve(0.3, 0.5, 0.5, "compatible", None, 1e-9)
    b = solve(0.3, 0.5, 0.5, "zero-gradient", None, 1e-9)
    assert abs(a.y[1, 0] - b.y[1, 0]) > 0.3


def test_negative_K_wall_minimum_of_regularity_factor():
    K, S, lam = -0.3, 1.5, -1.6
    sol = None
    for l in np.arange(-0.1, lam - 1e-9, -0.1):
        sol = solve(K, S, float(l), "compatible", sol, 1e-8)
    eta = np.linspace(0, 10, 400)
    d = 1 - 2 * K * sol.sol(eta)[0]
    assert np.argmin(d) == 0
    assert d[0] > 0
