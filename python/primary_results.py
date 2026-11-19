"""Numerical results of the pressure-compatible model used in the manuscript.

Writes CSV files to ../data/ and prints each block. Run from the python directory.
"""
import csv
import os
import numpy as np
from compatible_model import solve, growth_rate

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def continue_to(K, S, kind, lam_list, tol=1e-8):
    sol, last = None, None
    for lam in lam_list:
        trial = solve(K, S, float(lam), kind, sol, tol)
        if trial.status:
            break
        sol, last = trial, float(lam)
    return sol, last


def wall_shear(K, lam, sol):
    return (1 - 2 * K * lam) * np.hypot(sol.y[1, 0], sol.y[3, 0])


def table_critical():
    cases = [(-.05, 1.5), (-.1, .5), (-.1, 1.5), (-.15, 1), (-.2, .5), (-.2, 1.5),
             (-.2, 3), (-.3, .5), (-.3, 1.5), (-.3, 2), (-.4, 1.5)]
    rows = [("K", "S", "lambda_c_exact", "lambda_last_converged", "relative_gap", "min_D", "eta_at_min_D")]
    for K, S in cases:
        lc = 1 / (2 * K)
        fr = list(np.arange(0.05, 0.9, 0.05)) + list(np.arange(0.9, 0.995, 0.005)) + [0.996, 0.997, 0.998, 0.999]
        sol, last = continue_to(K, S, "compatible", [f * lc for f in fr])
        g = np.linspace(0, 25, 4001)
        D = 1 - 2 * K * sol.sol(g)[0]
        i = int(np.argmin(D))
        rows.append((K, S, lc, last, abs(last - lc) / abs(lc), float(D[i]), float(g[i])))
        print(rows[-1], flush=True)
    csv.writer(open(os.path.join(OUT, "compatible_critical_rates.csv"), "w", newline="")).writerows(rows)


def table_positive_and_baseline():
    rows = [("quantity", "closure", "value")]
    for kind in ("compatible", "zero-gradient"):
        for K in (0.3, 0.0):
            sol, _ = continue_to(K, 0.5, kind, np.arange(0.05, 0.5001, 0.05))
            rows.append(("Fp0 K=%g lam=0.5 S=0.5" % K, kind, sol.y[1, 0]))
            rows.append(("Gp0 K=%g lam=0.5 S=0.5" % K, kind, sol.y[3, 0]))
            print(rows[-2:], flush=True)
        sol, last = continue_to(0.3, 1.5, kind, -np.arange(0.1, 7.01, 0.1))
        g = np.linspace(0, 25, 4001)
        Fm = float(np.max(sol.sol(g)[0]))
        rows.append(("last converged lambda K=0.3 S=1.5", kind, last))
        rows.append(("max F at last lambda", kind, Fm))
        print(rows[-2:], flush=True)
    csv.writer(open(os.path.join(OUT, "compatible_baseline_and_positive_K.csv"), "w", newline="")).writerows(rows)


def table_closure():
    rows = [("lambda", "Fp0_zero_gradient", "Gp0_zero_gradient", "shear_zero_gradient",
             "Fp0_compatible", "Gp0_compatible", "shear_compatible")]
    K, S = -0.2, 1.5
    lams = [-0.5, -1.0, -1.5, -2.0, -2.4]
    store = {}
    for kind in ("zero-gradient", "compatible"):
        sol = None
        for lam in np.arange(-0.1, -2.4001, -0.1):
            sol = solve(K, S, float(lam), kind, sol)
            assert sol.status == 0
            if round(abs(lam) * 10) in (5, 10, 15, 20, 24):
                store[(kind, round(float(lam), 1))] = (sol.y[1, 0], sol.y[3, 0], wall_shear(K, lam, sol))
    for lam in lams:
        a, b = store[("zero-gradient", lam)], store[("compatible", lam)]
        rows.append((lam,) + tuple(a) + tuple(b))
        print(rows[-1], flush=True)
    csv.writer(open(os.path.join(OUT, "compatible_closure_comparison.csv"), "w", newline="")).writerows(rows)


def table_stability():
    rows = [("K", "S", "closure", "N", "L", "lambda_crossing_or_last", "status")]
    for K in (0.0, 0.3, -0.3):
        lams = np.arange(0, -4.01, -0.05) if K == 0 else (np.arange(0, -2.51, -0.05) if K > 0 else np.arange(0, -1.651, -0.05))
        for kind in ("zero-gradient", "compatible"):
            if K == 0 and kind == "compatible":
                continue
            for N, L in ((80, 30.0), (100, 40.0)):
                sol, prev, out = None, None, None
                for lam in lams:
                    sol = solve(K, 1.5, float(lam), kind, sol, 1e-9, L=L)
                    if sol.status:
                        out = (float(lam), "base solve failed")
                        break
                    m = growth_rate(K, kind, sol, N=N, L=L)
                    if prev is not None and prev[1] < 0 <= m:
                        out = (prev[0] + (0 - prev[1]) * (lam - prev[0]) / (m - prev[1]), "zero crossing")
                        break
                    prev = (float(lam), m)
                if out is None:
                    out = (prev[0], "no crossing; growth rate %.3f at last lambda" % prev[1])
                rows.append((K, 1.5, kind, N, L) + out)
                print(rows[-1], flush=True)
    csv.writer(open(os.path.join(OUT, "compatible_stability.csv"), "w", newline="")).writerows(rows)


if __name__ == "__main__":
    table_closure()
    table_positive_and_baseline()
    table_critical()
    table_stability()
