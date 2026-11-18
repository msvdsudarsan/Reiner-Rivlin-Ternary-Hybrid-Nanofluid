"""Continuation of the pressure-compatible system to 0.99*lambda_c for K < 0."""
import csv
import numpy as np
from compatible_pressure_check import solve, L

cases = [(-.05, 1.5), (-.1, .5), (-.1, 1.5), (-.15, 1), (-.2, .5), (-.2, 1.5),
         (-.2, 3), (-.3, .5), (-.3, 1.5), (-.3, 2), (-.4, 1.5)]
fractions = list(np.arange(0.05, 0.9, 0.05)) + list(np.arange(0.9, 0.9901, 0.01))
rows = [("K", "S", "lambda_c", "lambda_last", "min_1_minus_2KF", "eta_at_min", "Fp0", "Gp0", "max_rms_residual")]
for K, S in cases:
    lc = 1 / (2 * K)
    sol, last = None, None
    for f in fractions:
        trial = solve(K, S, f * lc, "compatible", sol, 1e-8)
        if trial.status:
            break
        sol, last = trial, f * lc
    grid = np.linspace(0, L, 4001)
    F = sol.sol(grid)[0]
    D = 1 - 2 * K * F
    i = int(np.argmin(D))
    rows.append((K, S, lc, last, D[i], grid[i], sol.y[1, 0], sol.y[3, 0], float(np.max(sol.rms_residuals))))
    print(rows[-1], flush=True)
csv.writer(open("compatible_pressure_sweep.csv", "w", newline="")).writerows(rows)
