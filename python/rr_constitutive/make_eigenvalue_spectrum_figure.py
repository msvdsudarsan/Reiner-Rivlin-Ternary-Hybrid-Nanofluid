"""Real momentum-block eigenvalues at the Example 1 baseline (K=0.3, lambda=0.5, S=0.5) from a
direct generalised eigensolve, plotted against the share of the squared eigenfunction amplitude
in the outer half of the domain. Writes ../../figures/eigenvalue_spectrum.pdf and .png"""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stability_tracked import real_modes
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(here, "..", "..", "figures")
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
for (N, L), mk, col in (((100, 30), "o", "#1f4e8c"), ((140, 40), "s", "#dc143c")):
    eta, modes = real_modes(0.3, 0.5, 0.5, N, L, gmax=3.0)
    m = [o for o in modes if abs(o[0]) > 1e-6 and o[2] < 1e-4]
    g = np.array([o[0] for o in m]); t = np.maximum([o[3] for o in m], 1e-18)
    for a, lim in zip(ax, (3.0, 0.8)):
        s = g < lim
        a.semilogy(g[s], t[s], mk, ms=4, mfc="none", color=col, label=r"$N=%d,\ \eta_\infty=%d$" % (N, L))
for a in ax:
    a.axvline(0.41857, color="k", ls="--", lw=0.9)
    a.set_xlabel(r"real eigenvalue $\gamma$"); a.set_ylabel("outer-half amplitude share")
ax[0].set_title(r"(a) $0<\gamma<3$", fontsize=10); ax[1].set_title(r"(b) zoom, $0<\gamma<0.8$", fontsize=10)
ax[1].legend(fontsize=8, loc="center right")
fig.tight_layout(); fig.savefig(os.path.join(out, "eigenvalue_spectrum.pdf")); fig.savefig(os.path.join(out, "eigenvalue_spectrum.png"), dpi=300)
