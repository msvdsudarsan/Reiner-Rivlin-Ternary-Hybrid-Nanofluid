"""Linear stability of the shrinking branch (similarity-preserving perturbations only).

For K = +0.3 and K = -0.3 at S = 1.5 this script
  1. locates the zero crossing lambda* of the smallest reliable eigenvalue by bisection, at
     several collocation sizes N and domain lengths einf;
  2. integrates the semi-discrete linearised system forward in time (implicit Euler) at
     lambda = -1.0 and compares the measured growth rate with the rate that the eigenvalue
     predicts.  Both calculations use the same spatial operator, so this is a consistency
     check of the eigenvalue extraction and not an independent validation.

Only the two representative cases above are covered.  No stability map over K or S is
produced: near the crossing at larger |K| the small eigenvalues move with the domain length
(the far-field operator has its spectrum edge at the origin), so a threshold map is not
claimed.
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from cheb_base import base_state_from_bvp
from stability_tracked import pencil, crossing, smallest_reliable

def growth_rate(K, lam, S=1.5, N=100, einf=30.0, dt=0.01, steps=1500, seed=1):
    eta, F, G, H, De, D2 = base_state_from_bvp(K, lam, S, N=N, einf=einf)
    A, B = pencil(K, eta, F, G, H, De, D2)
    lu = lu_factor(B + dt*A)
    rng = np.random.default_rng(seed); n = N + 1
    q = rng.standard_normal(2*n)
    for r in (0, n, n + N, N): q[r] = 0
    t, ln = [], []
    for s in range(steps):
        q = lu_solve(lu, B@q)
        if s % 10 == 0: t.append(s*dt); ln.append(np.log(np.linalg.norm(q)))
    i = len(t)//2
    return np.polyfit(t[i:], ln[i:], 1)[0]

if __name__ == "__main__":
    for K, lo, hi in ((0.3, -0.40, -0.50), (-0.3, -0.55, -0.65)):
        print("K=%+.1f, S=1.5" % K)
        for N, L in ((100, 30), (160, 50), (200, 60)):
            print("   N=%3d, einf=%2d:  lambda* = %.5f" % (N, L, crossing(K, 1.5, N, L, lo, hi)))
    print("\nimplicit-Euler consistency check (dt=0.01, same spatial operator)")
    for K, lam in ((0.3, -1.0), (-0.3, -0.8), (-0.3, -1.0)):
        g = smallest_reliable(K, lam)[0]
        print("  K=%+.1f lambda=%.1f  gamma_1=%+.6f  predicted rate=%+.5f  measured rate=%+.5f"
              % (K, lam, g, -np.log(1 + g*0.01)/0.01, growth_rate(K, lam)))
