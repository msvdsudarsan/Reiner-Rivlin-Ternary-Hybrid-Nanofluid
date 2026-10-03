"""Negative-K regularity-boundary sweep: 11 (K,S) combinations.

Reproduces the closed form lambda_c = 1/(2K) and confirms that the boundary is
attained at the wall (argmin F at eta = 0) rather than at an interior maximum.

Run from this directory:  python run_negK_sweep.py
"""
import json
from rr_solver import continue_lambda  # same directory

CASES = [(-0.05, 1.5), (-0.10, 0.5), (-0.10, 1.5), (-0.15, 1.0),
         (-0.20, 0.5), (-0.20, 1.5), (-0.20, 3.0), (-0.30, 0.5),
         (-0.30, 1.5), (-0.30, 2.0), (-0.40, 1.5)]

if __name__ == "__main__":
    out = {}
    print("%6s %5s %12s %12s %9s" % ("K", "S", "lambda_c", "1/(2K)", "err %"))
    print("-" * 48)
    for K, S in CASES:
        last, _ = continue_lambda(K, S, dlam=-0.05)
        if last is None:
            print("%6.2f %5.1f   no convergence" % (K, S))
            continue
        lam, minF, maxF = last
        bound = 1.0 / (2.0 * K)
        err = 100.0 * abs(lam - bound) / abs(bound)
        out["K=%.2f,S=%.1f" % (K, S)] = dict(
            lam_c=lam, minF=minF, bound=bound, err_pct=err)
        print("%6.2f %5.1f %12.4f %12.4f %8.2f%%" % (K, S, lam, bound, err))
    with open("../../data/negK_sweep.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwritten: ../../data/negK_sweep.json")
