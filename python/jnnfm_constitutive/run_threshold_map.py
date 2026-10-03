"""Stability threshold lambda* against the regularity boundary lambda_c, negative branch.

For each K the crossing of the smallest momentum eigenvalue through zero is located
by bisection, using the base state solved on the same Chebyshev grid as the
eigenproblem (cheb_base.smallest_real_eig). On the negative branch lambda_c = 1/(2K).

Expected output (S = 1.5):
     K    lambda*   lambda_c    R*
   0.00  -0.5332     none     0.000
  -0.05  -0.5471   -10.0000   0.055
  -0.10  -0.5605    -5.0000   0.112
  -0.15  -0.5733    -3.3333   0.172
  -0.20  -0.5849    -2.5000   0.234
  -0.30  -0.6036    -1.6667   0.362
  -0.40  -0.6117    -1.2500   0.489
  -0.60  -0.5777    -0.8333   0.693
  -0.80  -0.4912    -0.6250   0.786
  -1.00  -0.4039    -0.5000   0.808

The Newtonian disk (K = 0) already loses stability, so the instability is
hydrodynamic in origin; lambda* stays ahead of lambda_c over the whole range.
Full run: roughly 15 minutes on one core.
"""
from cheb_base import smallest_real_eig

BRACKETS = {0.0: (-0.3, -0.9), -0.6: (-0.55, -0.60), -0.8: (-0.45, -0.55), -1.0: (-0.30, -0.45)}

def crossing(K, S=1.5, iters=13):
    a, b = BRACKETS.get(K, (-0.5, -0.9))
    for _ in range(iters):
        m = 0.5 * (a + b)
        if smallest_real_eig(K, m, S) > 0: a = m
        else: b = m
    return 0.5 * (a + b)

if __name__ == "__main__":
    print("%6s %9s %9s %7s" % ("K", "lambda*", "lambda_c", "R*"))
    for K in (0.0, -0.05, -0.10, -0.15, -0.20, -0.30, -0.40, -0.60, -0.80, -1.00):
        ls = crossing(K)
        lc = 1 / (2 * K) if K else None
        print("%6.2f %9.4f %9s %7.3f" % (K, ls, ("%.4f" % lc) if lc else "none",
                                          2 * abs(K) * abs(ls)))


def suction_sweep(K=-1.0, S_values=(0.5, 1.0, 2.0, 3.0, 4.0, 6.0)):
    """Stability crossing at fixed K for several suction values (data/suction_sweep_K-1.csv).
    Expected at K=-1.0: S=0.5 -0.2462, 1.0 -0.3379, 2.0 -0.4351, 3.0 -0.4518, 4.0 -0.4549, 6.0 -0.4554;
    R* saturates near 0.91, so the instability still precedes lambda_c = -0.5."""
    brackets = {0.5: (-0.2, -0.35), 1.0: (-0.2, -0.35), 2.0: (-0.35, -0.45)}
    out = {}
    for S in S_values:
        a, b = brackets.get(S, (-0.45, -0.48))
        for _ in range(11):
            m = 0.5 * (a + b)
            if smallest_real_eig(K, m, S) > 0: a = m
            else: b = m
        out[S] = 0.5 * (a + b)
    return out
