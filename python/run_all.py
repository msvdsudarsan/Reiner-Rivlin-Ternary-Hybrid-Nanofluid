"""
run_all.py

Top-level driver reproducing the full Python/SciPy computational
workflow reported in main.tex: baseline solution (Example 1),
validation checks, continuation search (Example 2), stability
eigenvalues (Example 3, partial -- see notes below), and thermal/
species response (Example 4).

THIS IS THE PYTHON/SCIPY VERIFICATION IMPLEMENTATION, NOT MATLAB R2026a.
Every numerical value it produces is marked [Python/SciPy] in main.tex
and is PROVISIONAL pending the independent MATLAB R2026a cross-check
described in main.tex Section "Computational Environment".

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
    print(f"\n=== Example 2: continuation search, K={K}, S={S} ===")
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
    print("(slow -- shooting method; momentum ~1-2 min, thermal ~seconds, species: see notes)")
    roots_M, _, _ = find_momentum_eigenvalue(bs.mom_sol, K=p.K, eta_max=6.0,
                                              gamma_bracket=(-2, 2), n_scan=13, alpha_bracket=(-3, 3))
    print("  momentum eigenvalue root(s):", roots_M)

    eig_T = solve_thermal_eigenvalue(bs.mom_sol, effective_properties(p), gamma_guess=0.6, n=150)
    gamma_T = float(eig_T.p[0]) if eig_T.status == 0 else None
    print("  thermal eigenvalue (solve_bvp):", gamma_T)

    print("  species eigenvalue: NOT reliably obtained in this run (Sc_thnf too large for "
          "either solve_bvp or the shooting method within the available time budget) -- open item.")
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
