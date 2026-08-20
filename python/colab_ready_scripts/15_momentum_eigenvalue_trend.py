"""
15_momentum_eigenvalue_trend.py

Momentum-eigenvalue TREND along the shrinking branch at K=0.3, S=1.5,
using: (a) an adaptive-step, warm-started base-state continuation in
lambda (matching Algorithm 1's own step-halving philosophy), and
(b) the corrected Chebyshev-collocation eigenvalue solve of script 14
(with the V26 Gfpp bug fix), cross-checked at THREE independent
collocation resolutions at every lambda tested.

Honest result: resolved (three resolutions agree to within a few
percent) over lambda in [0, -3], revealing a non-monotonic dip toward
marginal stability near lambda ~ -2. Beyond lambda ~ -3, resolutions
disagree by 10% to qualitatively (spurious near-zero modes appear at
some resolutions and not others) -- traced to the accuracy of the
cubic-spline base-state representation used to build the collocation
matrices, which degrades as the base-state curvature grows approaching
the regularity boundary. This script reports both the resolved trend
AND the unresolved region explicitly; it does not extrapolate or guess
values in the unresolved region.
"""
import pickle
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import interp1d
from scipy.linalg import eig


# ---------------- base-state solver (adaptive continuation in lambda) ----------------
def mom_rhs(eta, y, Kv):
    F, Fp, G, Gp, H = y
    Delta = 1 - 2*Kv*F
    Delta = np.where(np.abs(Delta) < 1e-10, np.sign(Delta)*1e-10 + 1e-14, Delta)
    Fpp = (F**2 - G**2 + H*Fp - Kv*(Fp**2 - Gp**2)) / Delta
    Gpp = (2*F*G + H*Gp - 2*Kv*Fp*Gp) / Delta
    return np.vstack([Fp, Fpp, Gp, Gpp, -2*F])


def mom_bc(ya, yb, lamv, Sv):
    return np.array([ya[0]-lamv, ya[2]-1.0, ya[4]+Sv, yb[0], yb[2]])


def mom_guess(eta, lamv, Sv):
    d = np.exp(-eta)
    return np.vstack([lamv*d, -lamv*d, d, -d, np.full_like(eta, -Sv)])


def solve_one(sol_prev, K, lam, S, tol=1e-10, max_nodes=40000):
    if sol_prev is None:
        eta = np.linspace(0, 40.0, 800)
        y0 = mom_guess(eta, lam, S)
    else:
        eta, y0 = sol_prev.x, sol_prev.y
    return solve_bvp(lambda e, y: mom_rhs(e, y, K), lambda ya, yb: mom_bc(ya, yb, lam, S),
                      eta, y0, tol=tol, max_nodes=max_nodes, verbose=0)


def continue_in_K(lam, S, K_target, n=41):
    sol = None
    for Kstep in np.linspace(0, K_target, n)[1:]:
        s = solve_one(sol, Kstep, lam, S)
        if s.status != 0:
            return None
        sol = s
    return sol


def adaptive_continue_lambda(K, S, lam_start, lam_target, sol0, dlam0=0.02, dlam_min=1e-5):
    sol = sol0
    lam = lam_start
    dlam = -abs(dlam0) if lam_target < lam_start else abs(dlam0)
    while (dlam < 0 and lam > lam_target) or (dlam > 0 and lam < lam_target):
        lam_next = lam + dlam
        if (dlam < 0 and lam_next < lam_target) or (dlam > 0 and lam_next > lam_target):
            lam_next = lam_target
        s = solve_one(sol, K, lam_next, S)
        if s.status == 0:
            sol, lam = s, lam_next
            if abs(dlam) < dlam0:
                dlam = np.sign(dlam) * min(abs(dlam) * 1.5, dlam0)
        else:
            dlam /= 2.0
            if abs(dlam) < dlam_min:
                break
    return sol, lam


# ---------------- corrected Chebyshev momentum eigenvalue (from script 14) ----------------
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2., np.ones(N - 1), 2.]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1. / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def build_H_operator(Deta, n):
    D2 = Deta.copy()
    D2[0, :] = 0.0
    D2[0, 0] = 1.0
    Zsel = np.eye(n)
    Zsel[0, 0] = 0.0
    return -2.0 * np.linalg.solve(D2, Zsel)


def momentum_eigs_cheb_from_interp(interp_fn, K, eta_max, N):
    D, xi = cheb(N)
    eta = eta_max * (1 - xi) / 2
    Deta = D * (-2 / eta_max)
    D2eta = Deta @ Deta
    n = N + 1

    y0 = interp_fn(eta)
    F0, F0p, G0, G0p, H0 = y0[0], y0[1], y0[2], y0[3], y0[4]
    Delta0 = 1 - 2 * K * F0
    F0pp = (F0**2 - G0**2 + H0 * F0p - K * (F0p**2 - G0p**2)) / Delta0
    G0pp = (2 * F0 * G0 + H0 * G0p - 2 * K * F0p * G0p) / Delta0

    P = build_H_operator(Deta, n)
    I, Z = np.eye(n), np.zeros((n, n))

    A_FF = (np.diag(Delta0) @ D2eta + 2*K*np.diag(F0pp) - 2*K*np.diag(F0p) @ Deta
            - np.diag(H0) @ Deta + 2*np.diag(F0) - np.diag(F0p) @ P)
    A_FG = 2*K*np.diag(G0p) @ Deta - 2*np.diag(G0)
    A_GF = -4*K*np.diag(G0pp) + 2*K*np.diag(G0p) @ Deta - 2*np.diag(G0) - np.diag(G0p) @ P
    A_GG = (np.diag(Delta0) @ D2eta + 2*K*np.diag(F0pp) + 2*K*np.diag(F0p) @ Deta
            - np.diag(H0) @ Deta - 2*np.diag(F0))

    A = np.block([[A_FF, A_FG], [A_GF, A_GG]])
    B = np.block([[I, Z], [Z, I]])
    Afull, Bfull = -A.copy(), B.copy()
    m = n
    Afull[0, :] = 0; Afull[0, 0] = 1; Bfull[0, :] = 0
    Afull[m, :] = 0; Afull[m, m] = 1; Bfull[m, :] = 0
    Afull[N, :] = 0; Afull[N, 0:m] = Deta[N, :]; Bfull[N, :] = 0
    Afull[m+N, :] = 0; Afull[m+N, m+N] = 1; Bfull[m+N, :] = 0

    evals = eig(Afull, Bfull, right=False)
    evals = evals[np.isfinite(evals)]
    real_evals = np.sort(evals[np.abs(evals.imag) < 1e-6*np.maximum(1, np.abs(evals.real))].real)
    return real_evals


if __name__ == "__main__":
    K, S = 0.3, 1.5
    lam_targets = [0.5, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, -4.2]

    print("=== Base-state continuation (all should converge, status=0) ===")
    sol = continue_in_K(0.5, S, K)
    lam_now = 0.5
    base_states = {0.5: (sol.x.copy(), sol.y.copy())}
    for lt in lam_targets[1:]:
        sol, lam_reached = adaptive_continue_lambda(K, S, lam_now, lt, sol)
        print(f"  lambda={lam_reached:8.4f}  status={sol.status}  F'(0)={sol.y[1,0]:.6f}")
        base_states[lt] = (sol.x.copy(), sol.y.copy())
        lam_now = lam_reached

    print()
    print("=== Momentum eigenvalue at each lambda, 3-resolution cross-check ===")
    resolutions = [(100, 30.0), (120, 35.0), (150, 40.0)]
    for lt in lam_targets:
        x, y = base_states[lt]
        interp = interp1d(x, y, kind='cubic', axis=1, fill_value='extrapolate')
        vals = []
        for N, eta_max in resolutions:
            evs = momentum_eigs_cheb_from_interp(interp, K, eta_max, N)
            pos = evs[(evs > 1e-3) & (evs < 10)]
            vals.append(pos[0] if len(pos) else np.nan)
        spread = (max(vals) - min(vals)) / max(abs(np.nanmean(vals)), 1e-9)
        status = "RESOLVED" if spread < 0.05 else "NOT YET RESOLVED"
        print(f"  lambda={lt:6.2f}  vals={np.round(vals,6)}  spread={spread*100:5.1f}%  [{status}]")
