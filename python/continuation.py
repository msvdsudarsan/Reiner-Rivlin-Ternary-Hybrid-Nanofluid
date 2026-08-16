"""
continuation.py

Python (scipy.integrate.solve_bvp) implementation of branch continuation
in the shrinking parameter lambda for the MOMENTUM SUB-SYSTEM ONLY
(states F, F', G, G', H), which by Proposition 1 of main.tex is exactly
independent of the ternary hybrid nanoparticle loading and of the
reaction parameter -- so continuation only needs to be performed once,
not repeated for every (phi, beta) combination.

Two stages are used:
  1. Natural parameter continuation (step lambda, warm-start from the
     previous solution) while it converges cleanly.
  2. Once natural continuation starts struggling (a proxy for approaching
     a fold), switch to a genuine pseudo-arclength continuation (PALC) in
     the reduced (lambda, F'(0)) plane: lambda is treated as an unknown
     parameter (via scipy's solve_bvp `p` mechanism, exactly as used for
     the stability eigenvalue gamma), tied to the state via F(0)=lambda,
     and an arclength constraint on (lambda, F'(0)) replaces one of the
     natural-continuation degrees of freedom, allowing the solver to
     move BACKWARD in lambda if that is what the branch requires to pass
     through the fold.

This IS therefore a real (if reduced-plane, not full-function-space)
pseudo-arclength scheme, an upgrade over the natural-parameter-only
continuation shipped in the first version of the corresponding external reference script
(which is flagged there as unable to cross a fold). The external reference file has
been corrected to match this approach -- see the corresponding external reference script.

PYTHON IMPLEMENTATION -- not part of any other computational environment.
"""

import numpy as np
from scipy.integrate import solve_bvp
from dataclasses import dataclass, field
from typing import List

from parameters import Parameters


def _mom_rhs(eta, y, K):
    F, Fp, G, Gp, H = y
    Delta = 1 - 2*K*F
    Delta = np.where(np.abs(Delta) < 1e-10, np.sign(Delta)*1e-10 + 1e-14, Delta)
    Fpp = (F**2 - G**2 + H*Fp - K*(Fp**2 - Gp**2)) / Delta
    Gpp = (2*F*G + H*Gp - 2*K*Fp*Gp) / Delta
    Hp = -2*F
    return np.vstack([Fp, Fpp, Gp, Gpp, Hp])


def _mom_bc_natural(ya, yb, lam, S):
    return np.array([ya[0]-lam, ya[2]-1.0, ya[4]+S, yb[0], yb[2]])


def _init_guess(eta, lam, S):
    decay = np.exp(-eta)
    F0 = lam*decay
    G0 = decay
    H0 = -S + lam*(1-decay)
    return np.vstack([F0, -lam*decay, G0, -decay, H0])


@dataclass
class BranchPoint:
    lam: float
    Fp0: float
    Gp0: float
    sol: object  # scipy OdeSolution-bearing solve_bvp result
    stage: str   # 'natural' or 'palc'


def solve_momentum_natural(K, lam, S, eta_inf=15.0, n=400, eta_mesh=None, y_guess=None, tol=1e-9):
    if eta_mesh is None:
        eta_mesh = np.linspace(0, eta_inf, n)
    if y_guess is None:
        y_guess = _init_guess(eta_mesh, lam, S)
    sol = solve_bvp(lambda e,y: _mom_rhs(e,y,K), lambda ya,yb: _mom_bc_natural(ya,yb,lam,S),
                     eta_mesh, y_guess, tol=tol, max_nodes=100000, verbose=0)
    return sol


def _mom_rhs_palc(eta, y, params, K):
    return _mom_rhs(eta, y, K)


def _mom_bc_palc(ya, yb, params, S, prev_state, tangent, ds):
    lam = params[0]
    lam_prev, Fp0_prev = prev_state
    t_lam, t_Fp0 = tangent
    arclength = t_lam*(lam - lam_prev) + t_Fp0*(ya[1] - Fp0_prev) - ds
    return np.array([ya[0]-lam, ya[2]-1.0, ya[4]+S, yb[0], yb[2], arclength])


def solve_momentum_palc(K, S, eta_mesh, y_guess, lam_guess, prev_state, tangent, ds, tol=1e-9):
    sol = solve_bvp(
        lambda e,y,prm: _mom_rhs_palc(e,y,prm,K),
        lambda ya,yb,prm: _mom_bc_palc(ya,yb,prm,S,prev_state,tangent,ds),
        eta_mesh, y_guess, p=[lam_guess], tol=tol, max_nodes=100000, verbose=0)
    return sol


def trace_branch(K, S, lam_start=0.0, lam_stop_natural=-2.5, ds0=0.05, max_palc_steps=250, eta_inf=15.0):
    """Run natural continuation from lam_start down to lam_stop_natural,
    then switch to PALC and continue for up to max_palc_steps additional
    points, tracing through any fold encountered."""
    points: List[BranchPoint] = []

    # ---- stage 1: natural continuation ----
    lam = lam_start
    eta_mesh = np.linspace(0, eta_inf, 400)
    y_guess = _init_guess(eta_mesh, lam, S)
    dlam = -0.05
    while lam > lam_stop_natural - 1e-9:
        sol = solve_momentum_natural(K, lam, S, eta_mesh=eta_mesh, y_guess=y_guess)
        if sol.status != 0:
            dlam /= 2
            if abs(dlam) < 1e-4:
                break
            continue
        points.append(BranchPoint(lam, sol.y[1,0], sol.y[3,0], sol, 'natural'))
        eta_mesh, y_guess = sol.x, sol.y
        lam += dlam
        dlam = max(dlam*1.1, -0.08) if dlam < 0 else dlam  # keep step modest, negative

    # ---- stage 2: PALC through the fold and beyond ----
    if len(points) >= 2:
        p_prev, p_curr = points[-2], points[-1]
        t_lam = p_curr.lam - p_prev.lam
        t_Fp0 = p_curr.Fp0 - p_prev.Fp0
        norm = np.hypot(t_lam, t_Fp0)
        tangent = (t_lam/norm, t_Fp0/norm)
        prev_state = (p_curr.lam, p_curr.Fp0)
        eta_mesh, y_guess = p_curr.sol.x, p_curr.sol.y
        ds = ds0

        for _ in range(max_palc_steps):
            lam_pred = prev_state[0] + tangent[0]*ds
            sol = solve_momentum_palc(K, S, eta_mesh, y_guess, lam_pred, prev_state, tangent, ds)
            if sol.status != 0:
                ds /= 2
                if ds < 1e-5:
                    break
                continue
            lam_new = sol.p[0]
            Fp0_new = sol.y[1,0]
            points.append(BranchPoint(lam_new, Fp0_new, sol.y[3,0], sol, 'palc'))

            t_lam_new = lam_new - prev_state[0]
            t_Fp0_new = Fp0_new - prev_state[1]
            norm = np.hypot(t_lam_new, t_Fp0_new)
            if norm > 1e-12:
                tangent = (t_lam_new/norm, t_Fp0_new/norm)
            prev_state = (lam_new, Fp0_new)
            eta_mesh, y_guess = sol.x, sol.y
            ds = min(ds*1.15, 4*ds0)

    return points
