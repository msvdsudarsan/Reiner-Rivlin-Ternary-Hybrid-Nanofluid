"""Pressure-compatible closure versus the zero-pressure-gradient closure.

Solves the reduced Reiner-Rivlin momentum system for both closures and prints the
wall derivatives and the wall shear resultant (1 - 2*K*lam)*sqrt(F'(0)^2 + G'(0)^2).
"""
import numpy as np
from scipy.integrate import solve_bvp

L = 15.0


def rhs(x, y, K, kind):
    F, Fp, G, Gp, H = y
    D = 1 - 2 * K * F
    if kind == "compatible":
        num = F * F - G * G + H * Fp + K * (Fp * Fp + 3 * Gp * Gp)
    else:
        num = F * F - G * G + H * Fp - K * (Fp * Fp - Gp * Gp)
    return np.vstack((Fp, num / D, Gp, (2 * F * G + H * Gp - 2 * K * Fp * Gp) / D, -2 * F))


def solve(K, S, lam, kind, prev=None, tol=1e-8):
    if prev is None:
        x = np.linspace(0, L, 400)
        d = np.exp(-x)
        y = np.vstack((lam * d, -lam * d, d, -d, np.full_like(x, -S)))
    else:
        x, y = prev.x, prev.y
    bc = lambda a, b: np.array([a[0] - lam, a[2] - 1, a[4] + S, b[0], b[2]])
    return solve_bvp(lambda t, u: rhs(t, u, K, kind), bc, x, y, tol=tol, max_nodes=100000)


if __name__ == "__main__":
    K, S = -0.2, 1.5
    print("closure, lambda, F'(0), G'(0), wall shear resultant")
    for kind in ("zero-gradient", "compatible"):
        sol = None
        for lam in np.arange(-0.1, -2.4001, -0.1):
            sol = solve(K, S, float(lam), kind, sol)
            if sol.status:
                raise SystemExit("continuation failed at lambda = %.2f" % lam)
            if round(abs(lam) * 10) in (5, 10, 15, 20, 24):
                shear = (1 - 2 * K * lam) * np.hypot(sol.y[1, 0], sol.y[3, 0])
                print(kind, round(lam, 2), round(sol.y[1, 0], 4), round(sol.y[3, 0], 4), round(shear, 4))
