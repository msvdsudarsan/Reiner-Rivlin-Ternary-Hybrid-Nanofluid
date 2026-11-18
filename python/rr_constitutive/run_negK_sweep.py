"""Negative-K regularity-boundary sweep: 11 (K,S) combinations.

Reproduces the closed form lambda_c = 1/(2K) and confirms that the boundary is
attained at the wall (argmin F at eta = 0) rather than at an interior maximum.

Run from this directory:  python run_negK_sweep.py
"""
import json, sys
from rr_solver import continue_lambda  # same directory

CASES = [(-0.05, 1.5), (-0.10, 0.5), (-0.10, 1.5), (-0.15, 1.0),
         (-0.20, 0.5), (-0.20, 1.5), (-0.20, 3.0), (-0.30, 0.5),
         (-0.30, 1.5), (-0.30, 2.0), (-0.40, 1.5)]

if __name__ == "__main__":
    # Optional: python run_negK_sweep.py 0,3,4 part_a  -> run only those cases and write
    # ../../data/negK_sweep_part_a.json; without arguments all 11 cases run in sequence.
    idx = [int(v) for v in sys.argv[1].split(",")] if len(sys.argv) > 1 else range(len(CASES))
    tag = ("_" + sys.argv[2]) if len(sys.argv) > 2 else ""
    out = {}
    print("%6s %5s %12s %12s %9s" % ("K", "S", "lambda_c", "1/(2K)", "err %"))
    print("-" * 48)
    for i in idx:
        K, S = CASES[i]
        last, _ = continue_lambda(K, S, dlam=-0.05)
        if last is None:
            print("%6.2f %5.1f   no convergence" % (K, S), flush=True)
            continue
        lam, minF, maxF = last
        bound = 1.0 / (2.0 * K)
        err = 100.0 * abs(lam - bound) / abs(bound)
        out["K=%.2f,S=%.1f" % (K, S)] = dict(
            lam_c=lam, minF=minF, bound=bound, err_pct=err)
        print("%6.2f %5.1f %12.4f %12.4f %8.2f%%" % (K, S, lam, bound, err), flush=True)
        with open("../../data/negK_sweep%s.json" % tag, "w") as f:
            json.dump(out, f, indent=2)
    print("\nwritten: ../../data/negK_sweep%s.json" % tag)
