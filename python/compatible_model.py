"""Pressure-compatible Reiner-Rivlin rotating-disk similarity system.

p = P0(eta) + K r^2 (F'^2 + G'^2) in units rho = Omega = nu = 1, so that all three
momentum components hold. kind = 'compatible' (primary model) or 'zero-gradient'
(closure p_r = 0 used in earlier literature).
"""
import numpy as np
import scipy.linalg as sl
from scipy.integrate import solve_bvp

L_DEFAULT = 25.0


def numerators(F, Fp, G, Gp, H, K, kind):
    if kind == "compatible":
        nr = F * F - G * G + H * Fp + K * (Fp * Fp + 3 * Gp * Gp)
    else:
        nr = F * F - G * G + H * Fp - K * (Fp * Fp - Gp * Gp)
    ng = 2 * F * G + H * Gp - 2 * K * Fp * Gp
    return nr, ng


def rhs(x, y, K, kind):
    F, Fp, G, Gp, H = y
    D = 1 - 2 * K * F
    nr, ng = numerators(F, Fp, G, Gp, H, K, kind)
    return np.vstack((Fp, nr / D, Gp, ng / D, -2 * F))


def solve(K, S, lam, kind, prev=None, tol=1e-8, L=L_DEFAULT, max_nodes=200000):
    if prev is None:
        x = np.linspace(0, L, 400)
        d = np.exp(-x)
        y = np.vstack((lam * d, -lam * d, d, -d, np.full_like(x, -S)))
    else:
        x, y = prev.x, prev.y
    bc = lambda a, b: np.array([a[0] - lam, a[2] - 1, a[4] + S, b[0], b[2]])
    return solve_bvp(lambda t, u: rhs(t, u, K, kind), bc, x, y, tol=tol, max_nodes=max_nodes)


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    return D - np.diag(D.sum(1)), x


def growth_rate(K, kind, sol, N=80, L=30.0):
    """Largest real part of the spectrum of the linearised similarity-preserving
    operator, perturbations decaying at both ends (F, G zero at wall and far field,
    H zero at the wall). Negative means the base state is linearly stable."""
    D, x = cheb(N)
    eta = L * (1 - x) / 2
    D1 = -2 / L * D
    D2 = D1 @ D1
    n = N + 1
    Y = sol.sol(eta)
    q0 = np.hstack([Y[0], Y[2], Y[4]])

    def R(q):
        f, g, h = q[:n], q[n:2 * n], q[2 * n:]
        Dn = 1 - 2 * K * f
        nr, ng = numerators(f, D1 @ f, g, D1 @ g, h, K, kind)
        return np.hstack([Dn * (D2 @ f) - nr, Dn * (D2 @ g) - ng, D1 @ h + 2 * f])

    J = np.zeros((3 * n, 3 * n))
    e = 1e-6
    for j in range(3 * n):
        qp = q0.copy(); qp[j] += e
        qm = q0.copy(); qm[j] -= e
        J[:, j] = (R(qp) - R(qm)) / (2 * e)
    M = np.zeros((3 * n, 3 * n))
    M[:2 * n, :2 * n] = np.eye(2 * n)
    for k, rows in enumerate(([0, N], [0, N], [0])):
        for r in rows:
            row = k * n + r
            J[row, :] = 0
            J[row, row] = 1
            M[row, :] = 0
    w = sl.eig(J, M, right=False)
    w = w[np.isfinite(w)]
    return float(np.max(w.real))
