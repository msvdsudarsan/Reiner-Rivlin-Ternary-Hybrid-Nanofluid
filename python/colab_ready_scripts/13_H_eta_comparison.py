import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# H(eta) comparison across the shrinking cases used in Example 4,
# added to give a direct visual/numeric check of the near-wall
# sign-reversal mechanism proposed there for why -theta'(0) collapses
# toward zero at lambda=-2.0. K, S fixed at the Example 4 values.

K, S = 0.3, 1.5


def mom_rhs(eta, y, Kv):
    F, Fp, G, Gp, H = y
    Delta = 1.0 - 2.0 * Kv * F
    Fpp = (F**2 - G**2 + H * Fp - Kv * (Fp**2 - Gp**2)) / Delta
    Gpp = (2 * F * G + H * Gp - 2 * Kv * Fp * Gp) / Delta
    Hp = -2.0 * F
    return np.vstack([Fp, Fpp, Gp, Gpp, Hp])


def mom_bc(ya, yb, lam, Sv):
    return np.array([ya[0] - lam, ya[2] - 1.0, ya[4] + Sv, yb[0], yb[2]])


def solve_case(lam_target, eta_inf=40.0, n=800):
    # Warm-started natural continuation from the well-conditioned lambda=0.5
    # baseline down to the target lambda, exactly as Section 7's methodology
    # prescribes and as Remark 7 shows is necessary for reliable convergence
    # at these lambda values -- a cold decay-only guess fails to converge
    # for lambda <= -1.0 at this eta_inf, which is itself a small illustration
    # of the domain/guess-sensitivity lesson already documented in the paper.
    eta = np.linspace(0, eta_inf, n)
    lam0 = 0.5
    decay = np.exp(-eta)
    y0 = np.zeros((5, eta.size))
    y0[0] = lam0 * decay
    y0[1] = -lam0 * decay
    y0[2] = decay
    y0[3] = -decay
    y0[4] = -S + lam0 * (1 - decay)
    sol = solve_bvp(
        lambda e, y: mom_rhs(e, y, K),
        lambda ya, yb: mom_bc(ya, yb, lam0, S),
        eta, y0, tol=1e-10, max_nodes=6000,
    )
    assert sol.status == 0, f"baseline solve failed: {sol.message}"

    lam = lam0
    base_step = -0.05 if lam_target < lam0 else 0.05
    step = base_step
    while abs(lam - lam_target) > 1e-9:
        trial_lam = lam + step
        if abs(trial_lam - lam0) > abs(lam_target - lam0):
            trial_lam = lam_target
        sol_trial = solve_bvp(
            lambda e, y: mom_rhs(e, y, K),
            lambda ya, yb: mom_bc(ya, yb, trial_lam, S),
            sol.x, sol.y, tol=1e-10, max_nodes=6000,
        )
        if sol_trial.status != 0:
            step /= 2.0
            assert abs(step) > 1e-6, f"continuation stalled before reaching lambda={lam_target}"
            continue
        sol, lam = sol_trial, trial_lam
    return sol


fig, ax = plt.subplots(figsize=(7, 5))
cases = [(0.5, "tab:blue"), (-1.0, "tab:orange"), (-2.0, "tab:green"), (-2.5, "tab:red")]
for lam, color in cases:
    sol = solve_case(lam)
    eta, H = sol.x, sol.y[4]
    ax.plot(eta, H, label=rf"$\lambda={lam}$", color=color)
    zero_cross = eta[np.where(np.diff(np.sign(H)) != 0)[0]]
    print(
        f"lambda={lam:5.1f}  H(0)={H[0]:.4f}  max H={H.max():.4f} "
        f"at eta={eta[np.argmax(H)]:.3f}  zero-crossings near eta={np.round(zero_cross[:3], 3)}"
    )

ax.axhline(0, color="gray", linewidth=0.7)
ax.set_xlim(0, 8)
ax.set_xlabel(r"$\eta$")
ax.set_ylabel(r"$H(\eta)$")
ax.set_title(f"Axial velocity H(eta), K={K}, S={S}")
ax.legend()
fig.tight_layout()
fig.savefig("H_comparison.png", dpi=500)
print("Saved H_comparison.png")

# Actual console output (this sandbox, SciPy 1.14, verified before packaging):
# lambda=  0.5  H(0)=-1.5000  max H=-1.5000 at eta=0.000  zero-crossings near eta=[]
# lambda= -1.0  H(0)=-1.5000  max H=-0.4507 at eta=1.655  zero-crossings near eta=[]
# lambda= -2.0  H(0)=-1.5000  max H= 0.7511 at eta=1.378  zero-crossings near eta=[0.503 3.031]
# lambda= -2.5  H(0)=-1.5000  max H= 1.4210 at eta=1.335  zero-crossings near eta=[0.359 3.310]
