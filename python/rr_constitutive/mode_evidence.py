"""Eigenpair residual, outer-half amplitude share and domain/resolution convergence for the
momentum-block modes. Prints the smallest reliable eigenvalue at the Example 1 baseline
(K=0.3, lambda=0.5, S=0.5) and at the two stability crossings (K=+-0.3, S=1.5)."""
import numpy as np
from stability_tracked import real_modes
cases = (("baseline K=0.3 lam=0.5 S=0.5", 0.3, 0.5, 0.5),
         ("crossing K=+0.3 lam=-0.446905 S=1.5", 0.3, -0.446905, 1.5),
         ("crossing K=-0.3 lam=-0.603570 S=1.5", -0.3, -0.60357, 1.5))
for name, K, lam, S in cases:
    print(name)
    print("   N   einf   gamma_1      residual   outer-half share   next modes")
    for N, L in ((60, 20), (80, 30), (100, 30), (140, 40), (200, 60)):
        eta, out = real_modes(K, lam, S, N, L, gmax=60.0)
        if out is None: print("  base state failed", N, L); continue
        c = [o for o in out if abs(o[0]) > 1e-6 and o[2] < 1e-4]
        o = c[0]; nxt = ", ".join("%.4f" % x[0] for x in c[1:3])
        print("%4d  %4d  %+.6f   %.1e     %.2e          %s" % (N, L, o[0], o[2], o[3], nxt))
