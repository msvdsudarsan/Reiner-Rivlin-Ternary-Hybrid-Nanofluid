"""Figures for the pressure-compatible analysis (reads ../data, writes ../figures)."""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compatible_model import solve, growth_rate

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, FIG = os.path.join(HERE, "..", "data"), os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9.5, "legend.fontsize": 8, "xtick.labelsize": 8.5,
                     "ytick.labelsize": 8.5, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "serif", "mathtext.fontset": "stix"})
BLUE, ORANGE, GRAY = "#1f5fa8", "#d1661c", "#555555"
rows = list(csv.DictReader(open(os.path.join(DATA, "compatible_critical_rates.csv"))))
K = np.array([float(r["K"]) for r in rows])
lc = np.array([float(r["lambda_c_exact"]) for r in rows])
ll = np.array([float(r["lambda_last_converged"]) for r in rows])


def fig_boundary():
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    kk = np.linspace(-0.42, -0.04, 100)
    ax[0].plot(-kk, 1 / (2 * kk), color=GRAY, lw=1.4, label=r"closed form $1/(2K)$")
    ax[0].plot(-K, ll, "o", mfc="white", mec=BLUE, ms=6, label="continuation limit")
    ax[0].set_xlabel(r"$|K|=\mathrm{Wi}_2$")
    ax[0].set_ylabel(r"$\lambda_c$")
    ax[0].legend(frameon=False, loc="lower right")
    ax[0].set_title("(a) wall-attained boundary, $K<0$", fontsize=9)
    for Kv, lam, col, st, lab in ((0.3, -3.5, BLUE, "-", r"$K=+0.3,\ \lambda=-3.50$"),
                                  (-0.3, -1.65, ORANGE, "--", r"$K=-0.3,\ \lambda=-1.65$")):
        sol = None
        for l in -np.arange(0.1, abs(lam) + 1e-9, 0.05 if Kv < 0 else 0.1):
            sol = solve(Kv, 1.5, float(l), "compatible", sol)
        g = np.linspace(0, 8, 600)
        ax[1].plot(g, sol.sol(g)[0], st, color=col, lw=1.6, label=lab)
        ax[1].axhline(1 / (2 * Kv), color=col, lw=0.8, ls=":")
    ax[1].set_xlabel(r"$\eta$")
    ax[1].set_ylabel(r"$F(\eta)$")
    ax[1].set_title(r"(b) $F$ near termination, $S=1.5$", fontsize=9)
    ax[1].legend(frameon=False, loc="center right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "negative_K_boundary.pdf"))


def fig_collapse():
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    wi = -K
    ax[0].plot(wi, 2 * wi * np.abs(ll), "o", mfc="white", mec=BLUE, ms=6)
    ax[0].axhline(1, color=GRAY, lw=1.2)
    ax[0].set_ylim(0.99, 1.01)
    ax[0].set_xlabel(r"$\mathrm{Wi}_2=|\Psi_2|\Omega/\eta_0$")
    ax[0].set_ylabel(r"$2\,\mathrm{Wi}_2\,|\lambda_c|$")
    ax[0].set_title("(a) eleven computed cases", fontsize=9)
    e = np.linspace(0, 0.3, 50)
    for r, st, col in ((0.0, "-", BLUE), (0.05, "--", ORANGE), (0.10, ":", GRAY)):
        ax[1].plot(100 * e, 100 * np.sqrt(e ** 2 + r ** 2), st, color=col, lw=1.6, label=r"$\delta\eta_0/\eta_0=%d\%%$" % round(100 * r))
    ax[1].set_xlabel(r"$\delta\Psi_2/\Psi_2$ (%)")
    ax[1].set_ylabel(r"$\delta|\lambda_c|/|\lambda_c|$ (%)")
    ax[1].legend(frameon=False)
    ax[1].set_title("(b) propagated uncertainty", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "regularity_collapse.pdf"))


def fig_stability():
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=False)
    spec = ((0.3, ax[0], -2.5), (-0.3, ax[1], -1.65))
    out = []
    for Kv, a, end in spec:
        for kind, col, st, lab in (("compatible", BLUE, "-", "pressure-compatible"),
                                   ("zero-gradient", ORANGE, "--", r"zero gradient, $p_r=0$")):
            lim = end if (kind == "compatible" or Kv < 0) else -2.5
            sol, xs, ys = None, [], []
            for l in np.arange(0, lim - 1e-9, -0.1):
                sol = solve(Kv, 1.5, float(l), kind, sol, 1e-8, L=30.0)
                if sol.status != 0:
                    break
                xs.append(l)
                ys.append(growth_rate(Kv, kind, sol, N=80, L=30.0))
                out.append((Kv, kind, float(l), ys[-1]))
            a.plot(xs, ys, st, color=col, lw=1.6, label=lab)
        a.axhline(0, color=GRAY, lw=0.8)
        a.axvline(1 / (2 * Kv) if Kv < 0 else -3.5, color=GRAY, lw=0.8, ls=":")
        a.set_xlabel(r"$\lambda$")
        a.set_ylabel("largest growth rate")
        a.set_title(r"$K=%+.1f,\ S=1.5$" % Kv, fontsize=9)
    ax[0].legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "stability_corrected.pdf"))
    csv.writer(open(os.path.join(DATA, "compatible_growth_rate_curves.csv"), "w", newline="")).writerows(
        [("K", "closure", "lambda", "largest_growth_rate")] + out)


def fig_regime():
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    k = np.linspace(-0.5, -0.03, 200)
    ax.plot(k, 1 / (2 * k), color=BLUE, lw=1.6)
    ax.fill_between(k, 1 / (2 * k), -14, color=GRAY, alpha=0.18)
    ax.set_ylim(-14, 0)
    ax.set_xlabel(r"$K$")
    ax.set_ylabel(r"$\lambda$")
    ax.text(-0.3, -8, "wall screen fails", fontsize=8)
    ax.text(-0.15, -1.2, "wall screen passes", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "regime_map.pdf"))


if __name__ == "__main__":
    fig_regime()
    fig_boundary()
    fig_collapse()
    fig_stability()
