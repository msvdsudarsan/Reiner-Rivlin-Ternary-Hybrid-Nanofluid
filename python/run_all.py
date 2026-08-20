"""
run_all.py

Top-level driver reproducing the CORE numerical workflow reported in
Main.tex: baseline solution (Example 1), validation checks, a natural-
continuation search (Example 2, illustrative only -- see note below),
stability eigenvalues (Example 3, partial -- see notes below), and
thermal/species response (Example 4).

This is the core workflow only, not the full set of analyses in the
manuscript. The extended analyses -- the systematic 14-point (K,S)
broadened search and its bisection refinement of lambda_c, the
Puspanathan et al. parameter-matched and K-continuation comparison, the
corrected Chebyshev momentum-eigenvalue solve and its trend along the
branch, the viscosity/thermal-conductivity closure sensitivity checks,
the H(eta) mechanism comparison, and the reference-verification pass --
are each reproduced by their own dedicated script in
colab_ready_scripts/ (see that directory's README for the full list),
not by this driver.

THIS SCRIPT USED to produce every numerical value in the manuscript,
independently re-executed and confirmed by the corresponding author
(see Main.tex Section "Computational Environment").

Run: python3 run_all.py
(Takes several minutes; the stability eigenvalue searches in particular
are slow because a robust-but-slow shooting method was needed -- see
stability_shooting.py and the notes in README_python_notes.md.)
"""

import json
import numpy as np

from parameters import Parameters
from effective_properties import effective_properties
from base_state import solve_momentum, solve_energy, solve_species, solve_base_state_decoupled
from validation import newtonian_check, reiner_rivlin_singlephase_check, mesh_domain_independence, residual_norms
from continuation import solve_momentum_natural, _init_guess
from stability import solve_thermal_eigenvalue
from stability_shooting import find_momentum_eigenvalue, find_scalar_eigenvalue


def example1():
    print("\n=== Example 1: baseline solution ===")
    p = Parameters()
    props = effective_properties(p)
    bs = solve_base_state_decoupled(p, props, tol_mom=1e-10, tol_scalar=1e-11)
    print(f"F'(0)      = {bs.Fp0:.6f}")
    print(f"-G'(0)     = {-bs.Gp0:.6f}")
    print(f"-theta'(0) = {-bs.thp0:.6f}")
    print(f"-phi'(0)   = {-bs.php0:.6f}")
    return bs, props, p


def validation_suite():
    print("\n=== Validation ===")
    nc = newtonian_check()
    print("Newtonian check:", json.dumps(nc, indent=2))
    rr = reiner_rivlin_singlephase_check()
    print("Reiner-Rivlin single-phase check:", json.dumps(rr, indent=2))
    mesh = mesh_domain_independence()
    print("Mesh/domain independence:")
    for row in mesh:
        print(" ", row)
    res = residual_norms()
    print("Residual norms:", json.dumps(res, indent=2))
    return nc, rr, mesh, res


def example2(K=0.3, S=1.5):
    print(f"\n=== Example 2: continuation search, K={K}, S={S} (ILLUSTRATIVE ONLY) ===")
    print("  NOTE: this natural-continuation sweep uses a fixed step and no bisection")
    print("  refinement, so it approaches but does not pin down lambda_c to the precision")
    print("  reported in the manuscript (-5.026262); for that, and for the systematic")
    print("  14-point (K,S) broadened search, see 03_example2_lambda_c.py.")
    sing_val = 1 / (2 * K)
    lam = 0.0
    eta_mesh = np.linspace(0, 15, 400)
    y_guess = _init_guess(eta_mesh, lam, S)
    dlam = -0.05
    lams, Fp0s, maxFs = [], [], []
    while lam > -6.0:
        s = solve_momentum_natural(K, lam, S, eta_mesh=eta_mesh, y_guess=y_guess)
        if s.status != 0:
            print(f"  natural continuation FAILED near lambda={lam:.4f} "
                  f"(last converged max F={maxFs[-1]:.4f}, singular F={sing_val:.4f})")
            break
        eta_mesh, y_guess = s.x, s.y
        lams.append(lam); Fp0s.append(s.y[1, 0]); maxFs.append(s.y[0].max())
        lam += dlam
    print(f"  reached {len(lams)} points, last lambda={lams[-1]:.4f}, "
          f"max F(eta)={maxFs[-1]:.4f} (singular at {sing_val:.4f})")
    return np.array(lams), np.array(Fp0s), np.array(maxFs)


def example3_baseline(bs, p):
    print("\n=== Example 3: stability eigenvalues at Example 1 baseline ===")
    print("(slow -- shooting method; momentum ~1-2 min, thermal ~seconds)")
    print("  NOTE: this shooting search uses a small eta_max=6.0 for speed and is an")
    print("  illustrative cross-check, not the manuscript's authoritative momentum-")
    print("  eigenvalue calculation. The corrected system has several closely spaced")
    print("  modes (~0.42, ~0.47-0.55), and this coarse shooting search is not")
    print("  guaranteed to land on the smallest one -- see")
    print("  14_momentum_eigenvalue_chebyshev.py for the resolution-independent value")
    print("  (gamma_1^(M) = 0.41857) actually reported in the manuscript.")
    roots_M, _, _ = find_momentum_eigenvalue(bs.mom_sol, K=p.K, eta_max=6.0,
                                              gamma_bracket=(-2, 2), n_scan=13, alpha_bracket=(-3, 3))
    print("  momentum eigenvalue root(s) found by this illustrative search:", roots_M)

    eig_T = solve_thermal_eigenvalue(bs.mom_sol, effective_properties(p), gamma_guess=0.6, n=150)
    gamma_T = float(eig_T.p[0]) if eig_T.status == 0 else None
    print("  thermal eigenvalue (solve_bvp):", gamma_T)

    print("  species eigenvalue: NOT obtained by direct shooting/solve_bvp search here")
    print("  (Sc_thnf too large for either method within a practical time budget). This")
    print("  is not simply an open item: a corrected analytical far-field threshold")
    print("  (gamma_thr = +Sc_thnf*beta, corrected in V33 from an earlier sign error)")
    print("  and an independent Chebyshev collocation cross-check both find no")
    print("  resolution-stable discrete eigenvalue in the regime tested; see")
    print("  08_broad_search_and_species_threshold.py and 10_chebyshev_species_search.py,")
    print("  and the manuscript's species-threshold Remark for the full account.")
    return roots_M, gamma_T


def example4():
    print("\n=== Example 4: thermal/species response ===")
    K, S = 0.3, 1.5
    lam_values = [0.5, -1.0, -2.0]
    phi_values = [0.0, 0.01, 0.03]
    beta_values = [0.0, 0.5]
    rows = []
    for lam in lam_values:
        p = Parameters(K=K, lam=lam, S=S)
        mom = solve_momentum(p, tol=1e-10)
        for phi in phi_values:
            p2 = Parameters(K=K, lam=lam, S=S, phi1=phi, phi2=phi, phi3=phi)
            props = effective_properties(p2)
            en = solve_energy(mom, props.Pr_thnf, tol=1e-11)
            for beta in beta_values:
                sp = solve_species(mom, props.Sc_thnf, beta, tol=1e-11)
                row = {'lambda': lam, 'phi': phi, 'beta': beta, 'Pr_thnf': props.Pr_thnf,
                       'Sc_thnf': props.Sc_thnf, 'negthetap0': float(-en.y[1, 0]),
                       'negphip0': float(-sp.y[1, 0])}
                rows.append(row)
                print(" ", row)
    return rows


if __name__ == '__main__':
    bs, props, p = example1()
    validation_suite()
    lams, Fp0s, maxFs = example2()
    example3_baseline(bs, p)
    example4()
    print("\nDone. See ../data/ for saved CSV/JSON outputs from the original run, "
          "and ../figures/ for the generated PDF figures.")
