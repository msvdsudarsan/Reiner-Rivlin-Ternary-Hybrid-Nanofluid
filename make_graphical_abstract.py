"""Graphical abstract for the JNNFM manuscript, drawn from data/lamstar_map.csv.

Usage:  python make_graphical_abstract.py   (writes Graphical_Abstract.png/.pdf here)
Size 10 x 4 in at 243 dpi = 2430 x 972 px (ratio 2.5, above the journal minimum).
K = 0 is the Newtonian control: it has a stability threshold but no finite lambda_c.
"""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
d = np.genfromtxt(os.path.join(here, "..", "data", "lamstar_map.csv"), delimiter=",", names=True)
absK, lstar = np.abs(d["K"]), np.abs(d["lambda_star"])
fig = plt.figure(figsize=(10, 4))
fig.text(0.035, 0.84, "Reiner−Rivlin rotating disk", fontsize=15.5, fontweight="bold")
kw = dict(fontsize=14.5)
fig.text(0.035, 0.67, r"$\mu_c=\Psi_2\ \ \Rightarrow\ \ K=\Psi_2\Omega/\eta_0$", **kw)
fig.text(0.035, 0.52, r"reported $\Psi_2<0\ \ \Rightarrow\ \ K<0$", **kw)
fig.text(0.035, 0.37, r"boundary at the wall:  $\lambda_c=1/(2K)$", color="#dc143c", **kw)
fig.text(0.035, 0.22, r"instability first:  $|\lambda^*|<|\lambda_c|$", color="#1f4e8c", **kw)
fig.text(0.035, 0.105, r"for all nonzero $K$ tested, $-1\leq K<0$;  $\mathrm{Wi}_2=|K|$", fontsize=11, style="italic")
fig.text(0.035, 0.04, r"$K=0$ is the Newtonian control (no finite $\lambda_c$)", fontsize=11, style="italic")
ax = fig.add_axes([0.53, 0.14, 0.44, 0.80])
x = np.linspace(0.045, 1.05, 400)
ax.plot(x, 1/(2*x), color="#dc143c", lw=2.2, label=r"regularity boundary $|\lambda_c|$")
ax.plot(absK, lstar, "-o", color="#1f4e8c", lw=2.2, ms=5.5, label=r"stability threshold $|\lambda^*|$")
xs = np.linspace(0.045, 1.0, 300)
ax.fill_between(xs, np.interp(xs, absK, lstar), 1/(2*xs), color="#dc143c", alpha=0.13, lw=0, label="steady but unstable")
ax.annotate("Newtonian control ($K=0$)", xy=(0, lstar[0]), xytext=(0.09, 0.33), fontsize=8.5,
            arrowprops=dict(arrowstyle="-", lw=0.8, color="0.3"), va="center")
ax.set_yscale("log"); ax.set_xlim(0, 1.05); ax.set_ylim(0.27, 12)
ax.set_xlabel(r"$|K|=\mathrm{Wi}_2$", fontsize=11); ax.set_ylabel(r"critical $|\lambda|$", fontsize=11)
ax.grid(True, which="both", alpha=0.3); ax.legend(loc="upper right", fontsize=9)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(here, "Graphical_Abstract." + ext), dpi=243)
