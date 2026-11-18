"""Rheometry-to-model workflow for the restricted Pi=0 Reiner-Rivlin disk model.

IMPORTANT: INPUT below is explicitly synthetic. Replace it only with rheometry measurements for the same characterised sample.
The eta0 input must be the constant-viscosity parameter consistent with the model;
if viscosity is rate-dependent, do not substitute an apparent viscosity without
reformulating the closure. Report the shear-rate range, protocol and uncertainty.

The wall factor is D_w = 1 - 2 K lambda. D_w > 0 is only a wall screen.
It is NOT a proof of global regularity: the full BVP must also satisfy
min_eta(1 - 2 K F(eta)) > 0. For K>0 an interior maximum can be decisive.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"
FIG_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

INPUT = dict(
    R=0.025, rpm=60.0, rho=1000.0,
    eta0=1.0, Psi2=-0.02,                 # synthetic values, not measurements
    rel_sd_eta0=0.05, rel_sd_Psi2=0.15,   # independent 1-sigma CVs
    rel_sd_rpm=0.0, rel_sd_R=0.0,
    lambda_disk=0.0, S=0.0,
    label="synthetic illustration",
)

def rhs_factory(K):
    def rhs(eta, y):
        F, Fp, G, Gp, H = y
        D = 1.0 - 2.0*K*F
        return np.vstack((
            Fp,
            (F**2-G**2+H*Fp-K*(Fp**2-Gp**2))/D,
            Gp,
            (2*F*G+H*Gp-2*K*Fp*Gp)/D,
            -2*F,
        ))
    return rhs

def solve_state(K, lam=0.0, suction=0.0, guess=None, Z=40.0):
    """Solve the Pi=0 similarity BVP; returns scipy solve_bvp result."""
    def bc(a, b):
        return np.array([a[0]-lam, a[2]-1.0, a[4]+suction, b[0], b[2]])
    if guess is None:
        x = np.linspace(0.0, Z, 800)
        d = np.exp(-x)
        y = np.vstack((0.2*x*d, 0.2*(1-x)*d, d, -d, -suction*(1-d)))
    else:
        x, y = guess.x, guess.y
    return solve_bvp(rhs_factory(K), bc, x, y, tol=1e-9, max_nodes=200000)

def torque_coefficient_table(Kmax=0.6, n=25):
    """Continuation table for |G'(0)| at lambda=S=0, K<=0."""
    ks = np.linspace(0.0, -Kmax, n)
    sol = solve_state(0.0)
    if sol.status != 0:
        raise RuntimeError(f"Newtonian seed BVP failed: {sol.message}")
    rows = [(0.0, abs(sol.y[3, 0]))]
    for k in ks[1:]:
        trial = solve_state(k, guess=sol)
        if trial.status != 0:
            break
        Dmin = float(np.min(1.0-2.0*k*trial.y[0]))
        if Dmin <= 1e-6:
            break
        sol = trial
        rows.append((float(k), float(abs(sol.y[3, 0]))))
    if len(rows) < 5:
        raise RuntimeError("Insufficient converged K-grid for uncertainty interpolation.")
    return np.asarray(rows)

def normal_positive(mean, cv, size, rng):
    """Normal measurement-error model, rejecting nonphysical eta<=0 draws."""
    x = mean*(1.0 + cv*rng.standard_normal(size))
    bad = x <= 0
    while np.any(bad):
        x[bad] = mean*(1.0 + cv*rng.standard_normal(np.count_nonzero(bad)))
        bad = x <= 0
    return x

def q025_50_975(x):
    return np.percentile(np.asarray(x), [2.5, 50.0, 97.5]).tolist()

def main():
    p = INPUT
    if p["eta0"] <= 0 or p["rho"] <= 0 or p["R"] <= 0 or p["rpm"] <= 0:
        raise ValueError("R, rpm, rho and eta0 must be positive.")
    if p["rel_sd_eta0"] < 0 or p["rel_sd_Psi2"] < 0:
        raise ValueError("Relative standard deviations must be non-negative.")

    Om = 2*np.pi*p["rpm"]/60.0
    eta0, psi2, rho, R = p["eta0"], p["Psi2"], p["rho"], p["R"]
    K = psi2*Om/eta0
    lam = p["lambda_disk"]
    if not (K < 0):
        raise ValueError("This workflow currently implements the negative-K branch only.")
    if abs(lam) > 1e-14 or abs(p["S"]) > 1e-14:
        raise ValueError("Torque continuation table is validated only for lambda=0 and S=0; do not change these inputs.")
    Dwall = 1.0 - 2.0*K*lam
    lamc = np.inf if K == 0 else 1.0/(2.0*K)

    print(f"Input label: {p['label']} (all values must be classified as measured or synthetic)")
    print(f"Omega={Om:.8g} rad/s; K=Psi2*Omega/eta0={K:.8g}")
    print(f"lambda={lam:.8g}; wall factor D_w=1-2Klambda={Dwall:.8g}")
    if K < 0:
        print(f"Wall singularity value lambda_c=1/(2K)={lamc:.8g}.")
        print("Wall screen:", "PASS" if Dwall > 0 else "FAIL")
        print("PASS is necessary only; global regularity requires the full BVP denominator check.")
    elif K > 0:
        print("K>0: wall screen alone cannot detect an interior singularity; solve and inspect the full profile.")
    else:
        print("K=0: Newtonian limit; no cross-viscosity denominator singularity.")

    # Solve base state and explicitly check the whole denominator profile.
    tab = torque_coefficient_table()
    # tab is generated in decreasing K order; sort ascending for np.interp.
    tab = tab[np.argsort(tab[:, 0])]
    kvals, cgvals = tab[:, 0], tab[:, 1]
    if K < kvals.min() or K > kvals.max():
        raise ValueError(f"Nominal K={K:.5g} lies outside computed torque table [{kvals.min():.3g},{kvals.max():.3g}].")
    cg_nom = float(np.interp(K, kvals, cgvals))
    cg_newton = float(cgvals[np.argmin(np.abs(kvals))])
    nominal_state = solve_state(K, lam=lam, suction=p["S"])
    if nominal_state.status != 0:
        raise RuntimeError(f"Nominal BVP failed: {nominal_state.message}")
    Dmin_nominal = float(np.min(1.0-2.0*K*nominal_state.y[0]))
    if Dmin_nominal <= 1e-6:
        raise RuntimeError("Nominal BVP reaches the singular denominator; torque output withheld.")
    print(f"Full-profile denominator minimum: {Dmin_nominal:.8g} (global numerical check for this computed state)")
    nu0 = eta0/rho
    M_nom = 0.5*np.pi*eta0*R**4*Om*np.sqrt(Om/nu0)*cg_nom
    M_newton_same_eta = 0.5*np.pi*eta0*R**4*Om*np.sqrt(Om/nu0)*cg_newton
    print(f"|G'(0)|: Newtonian={cg_newton:.8f}; at K={K:.8g}: {cg_nom:.8f}")
    print(f"Torque per face at nominal eta: {M_nom:.8e} N m")
    print(f"Torque change relative to K=0 at same eta: {100*(cg_nom/cg_newton-1):+.3f}%")

    # Independent normal measurement errors; fixed seed for reproducibility.
    rng = np.random.default_rng(20261116)
    N = 50000
    eta_s = normal_positive(eta0, p["rel_sd_eta0"], N, rng)
    psi_s = psi2*(1.0 + p["rel_sd_Psi2"]*rng.standard_normal(N))
    rpm_s = normal_positive(p["rpm"], p["rel_sd_rpm"], N, rng)
    R_s = normal_positive(R, p["rel_sd_R"], N, rng)
    Om_s = 2*np.pi*rpm_s/60.0
    K_s = psi_s*Om_s/eta_s
    lamc_s = 1.0/(2.0*K_s)
    wall_s = 1.0 - 2.0*K_s*lam
    # K-dependent torque coefficient is only available on the computed negative-K grid.
    in_grid = (K_s >= kvals.min()) & (K_s <= kvals.max()) & (K_s < 0)
    cg_s = np.interp(K_s[in_grid], kvals, cgvals)
    # Correct sample-wise kinematic viscosity: nu_i=eta_i/rho.
    M_s = 0.5*np.pi*eta_s[in_grid]*R_s[in_grid]**4*Om_s[in_grid] * \
          np.sqrt(Om_s[in_grid]*rho/eta_s[in_grid]) * cg_s
    valid_fraction = float(np.mean(in_grid))
    print(f"Monte Carlo N={N}; negative-K samples within solved grid={valid_fraction:.3%}")
    print("K 2.5/50/97.5 percentiles:", np.round(q025_50_975(K_s), 6))
    print("lambda_c 2.5/50/97.5 percentiles:", np.round(q025_50_975(lamc_s), 6))
    print("Torque [N m] 2.5/50/97.5 percentiles on in-grid draws:", np.array2string(
        np.asarray(q025_50_975(M_s)), formatter={"float_kind":lambda x:f"{x:.6e}"}))
    print(f"Wall-screen pass fraction at lambda={lam}: {np.mean(wall_s>0):.3%}")
    # First-order sensitivity, with Omega and R held exact unless their CV is supplied.
    cvK = np.sqrt(p["rel_sd_Psi2"]**2 + p["rel_sd_rpm"]**2 + p["rel_sd_eta0"]**2)
    print(f"First-order relative SD(K) and |lambda_c|: {cvK:.3%} (independent errors)")

    # Regime map: wall singularity locus, not a proved global existence boundary.
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    kk = np.linspace(-1.0, -0.05, 400)
    lc = 1.0/(2.0*kk)
    ax.plot(kk, lc, lw=1.8, label=r"wall singularity: $1-2K\lambda=0$")
    ax.fill_between(kk, -6.0, lc, alpha=0.15, label="wall-screen failure region")
    ax.scatter([0.3, -0.3], [-0.4469, -0.6036], marker="o",
               facecolors="none", edgecolors="black",
               label=r"reported numerical stability crossings ($S=1.5$)")
    ax.axvline(0.0, color="0.4", lw=0.8, ls=":")
    ax.set_xlim(-1.0, 0.5); ax.set_ylim(-6.0, 1.0)
    ax.set_xlabel(r"$K=\Psi_2\Omega/\eta_0$")
    ax.set_ylabel(r"shrinking parameter $\lambda$")
    ax.text(0.53, 0.20, "$K>0$: inspect the full computed profile", 
            transform=ax.transAxes, fontsize=7, ha="left", va="bottom")
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG_DIR/"regime_map.png", dpi=300)
    fig.savefig(FIG_DIR/"regime_map.pdf")

    out = dict(
        input=p, Omega_rad_s=Om, K_nominal=K, lambda_wall=lamc,
        wall_factor_at_requested_lambda=Dwall,
        full_profile_denominator_min_nominal=Dmin_nominal,
        torque_coefficient_newtonian=cg_newton, torque_coefficient_nominal=cg_nom,
        torque_nominal_Nm=M_nom, torque_newtonian_same_eta_Nm=M_newton_same_eta,
        mc_seed=20261116, mc_draws=N, mc_in_grid_fraction=valid_fraction,
        mc_K_percentiles=q025_50_975(K_s),
        mc_lambda_wall_percentiles=q025_50_975(lamc_s),
        mc_torque_percentiles_Nm=q025_50_975(M_s),
        mc_wall_screen_pass_fraction=float(np.mean(wall_s>0)),
        first_order_relative_sd_K_and_lambda_wall=cvK,
        warning="Synthetic inputs; wall screen is not a proof of global regularity; torque MC uses sample-specific nu=eta/rho."
    )
    (DATA_DIR/"rheometry_workflow_output.json").write_text(json.dumps(out, indent=2))
    print(f"Wrote {FIG_DIR/'regime_map.png'} and {DATA_DIR/'rheometry_workflow_output.json'}")

if __name__ == "__main__":
    main()
