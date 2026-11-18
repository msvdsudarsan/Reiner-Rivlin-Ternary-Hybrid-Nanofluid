"""Smallest reliable momentum eigenvalue along the shrinking branch (K=+0.3 and K=-0.3, S=1.5).
Writes ../../data/stability_branch.csv and ../../figures/stability_corrected.pdf"""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stability_tracked import smallest_reliable
here = os.path.dirname(os.path.abspath(__file__))
cases = ((0.3, np.linspace(0, -4, 41), -0.4469, 1/(2*0.3)*-1 if False else None),
         (-0.3, np.linspace(0, -1.6, 33), -0.6036, 1/(2*-0.3)))
rows = []; res = {}
for K, lams, lstar, lc in cases:
    g = []
    for lam in lams:
        r = smallest_reliable(K, lam, 1.5, 100, 30.0)
        g.append(np.nan if r is None else r[0]); rows.append((K, lam, g[-1]))
    res[K] = (lams, np.array(g))
    print(K, [("%.2f" % l, "%.4f" % v) for l, v in zip(lams[::8], g[::8])], "last", g[-1])
np.savetxt(os.path.join(here, "..", "..", "data", "stability_branch.csv"), rows, delimiter=",",
           header="K,lambda,gamma_smallest_reliable", comments="", fmt="%.6g")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for a, (K, lams, lstar, lc) in zip(ax, cases):
    l, g = res[K]
    a.plot(l, g, "-", color="#1f4e8c", lw=1.8)
    a.axhline(0, color="k", lw=0.6); a.axvline(lstar, color="#1f4e8c", ls="--", lw=1)
    a.axvspan(l.min(), lstar, color="#1f4e8c", alpha=0.08)
    if lc is not None and lc >= l.min(): a.axvline(lc, color="#dc143c", ls=":", lw=1.4)
    a.set_xlabel(r"$\lambda$"); a.set_ylabel(r"smallest eigenvalue $\gamma_1$")
    a.set_title(r"$K=%+.1f,\ S=1.5$" % K, fontsize=10); a.invert_xaxis()
fig.tight_layout(); fig.savefig(os.path.join(here, "..", "..", "figures", "stability_corrected.pdf"))
