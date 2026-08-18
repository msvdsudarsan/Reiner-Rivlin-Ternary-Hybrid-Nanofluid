"""
stability.py

Python (scipy.integrate.solve_bvp) implementation of the linearized
temporal stability eigenvalue problem.

CORRECTION FOUND DURING THIS IMPLEMENTATION: main.tex's linearized
equations (17)-(20) show that the THETA and PHI perturbation equations
do not couple to each other at all (each is
forced only by the momentum perturbation Hf, not by one another). The
"scalar-transport block" of Proposition 2 therefore further separates
into two INDEPENDENT diagonal blocks -- a thermal block and a species
block -- each with its OWN eigenvalue, rather than a single block
sharing one eigenvalue as originally (incorrectly) implemented. This
module solves three independent eigenvalue problems:

  1. momentum block:  Ff, Ff', Gf, Gf', Hf      -> eigenvalue gamma_M
  2. thermal block:   Theta, Theta'              -> eigenvalue gamma_T
  3. species block:   Phi, Phi'                  -> eigenvalue gamma_C

Each is solved with the wall-vanishing condition (perturbation is zero
at the wall, since wall values are prescribed/not perturbed) plus a
derivative normalization (Weidman-type relaxation), and decay at the
truncated far field.
"""

import numpy as np
from scipy.integrate import solve_bvp

from parameters import Parameters
from effective_properties import EffectiveProperties


# ---------------------------------------------------------------- momentum
def _momentum_rhs(eta, z, params, base_interp, K):
    gamma = params[0]
    y0 = base_interp(eta)
    F0, F0p, G0, G0p, H0 = y0[0], y0[1], y0[2], y0[3], y0[4]
    Delta = 1 - 2 * K * F0
    Delta = np.where(np.abs(Delta) < 1e-10, np.sign(Delta) * 1e-10 + 1e-14, Delta)
    F0pp = (F0**2 - G0**2 + H0 * F0p - K * (F0p**2 - G0p**2)) / Delta
    G0pp = (2 * F0 * G0 + H0 * G0p - 2 * K * F0p * G0p) / Delta

    Ff, Ffp, Gf, Gfp, Hf = z
    Ffpp = (2*K*Ff*F0pp - 2*K*F0p*Ffp + 2*K*G0p*Gfp - gamma*Ff
            + 2*F0*Ff - 2*G0*Gf + H0*Ffp + Hf*F0p) / Delta
    Gfpp = (2*K*Ff*G0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p - gamma*Gf
            + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p) / Delta
    Hfp = -2*Ff
    return np.vstack([Ffp, Ffpp, Gfp, Gfpp, Hfp])


def _momentum_bc(za, zb, params):
    return np.array([za[0], za[1] - 1.0, za[2], za[4], zb[0], zb[2]])


def solve_momentum_eigenvalue(base_sol, p: Parameters, gamma_guess=0.5, n=300, tol=1e-6, max_nodes=500000):
    eta = np.linspace(base_sol.x[0], base_sol.x[-1], n)
    base_interp = base_sol.sol
    # initial guess consistent with the BCs Ff(0)=0, Ff'(0)=1, Gf(0)=0, Hf(0)=0:
    # Ff(eta)=eta*exp(-eta) satisfies Ff(0)=0, Ff'(0)=1 exactly (previous version
    # used -decay for Ff', which had the WRONG SIGN relative to the Ff'(0)=1
    # normalization BC and was the actual cause of the solve_bvp mesh-blowup
    # failures seen initially -- not numerical stiffness as first suspected).
    Ff0 = eta*np.exp(-eta)
    Ffp0 = (1-eta)*np.exp(-eta)
    decay = np.exp(-eta)
    z0 = np.vstack([Ff0, Ffp0, 0.3*eta*decay, 0.3*(1-eta)*decay, np.zeros_like(eta)])
    sol = solve_bvp(
        lambda e, z, prm: _momentum_rhs(e, z, prm, base_interp, p.K),
        _momentum_bc, eta, z0, p=[gamma_guess], tol=tol, max_nodes=max_nodes, verbose=0)
    return sol


# ----------------------------------------------------------------- thermal
def _thermal_rhs(eta, z, params, base_interp, Pr):
    gamma = params[0]
    y0 = base_interp(eta)
    H0 = y0[4]
    Th, Thp = z
    Thpp = Pr * H0 * Thp - gamma * Th
    return np.vstack([Thp, Thpp])


def _thermal_bc(za, zb, params):
    return np.array([za[0], za[1] - 1.0, zb[0]])


def solve_thermal_eigenvalue(base_sol, props: EffectiveProperties, gamma_guess=0.5, n=300, tol=1e-6, max_nodes=200000):
    eta = np.linspace(base_sol.x[0], base_sol.x[-1], n)
    base_interp = base_sol.sol
    # guess consistent with Theta(0)=0, Theta'(0)=1 (see fix note in
    # solve_momentum_eigenvalue above -- same bug existed here)
    z0 = np.vstack([eta*np.exp(-eta), (1-eta)*np.exp(-eta)])
    sol = solve_bvp(
        lambda e, z, prm: _thermal_rhs(e, z, prm, base_interp, props.Pr_thnf),
        _thermal_bc, eta, z0, p=[gamma_guess], tol=tol, max_nodes=max_nodes, verbose=0)
    return sol


# ----------------------------------------------------------------- species
def _species_rhs(eta, z, params, base_interp, Sc, beta):
    gamma = params[0]
    y0 = base_interp(eta)
    H0 = y0[4]
    Ph, Php = z
    Phpp = Sc * H0 * Php + Sc * beta * Ph - gamma * Ph
    return np.vstack([Php, Phpp])


def _species_bc(za, zb, params):
    return np.array([za[0], za[1] - 1.0, zb[0]])


def solve_species_eigenvalue(base_sol, props: EffectiveProperties, p: Parameters, gamma_guess=0.5, n=300, tol=1e-6, max_nodes=200000):
    eta = np.linspace(base_sol.x[0], base_sol.x[-1], n)
    base_interp = base_sol.sol
    z0 = np.vstack([eta*np.exp(-eta), (1-eta)*np.exp(-eta)])
    sol = solve_bvp(
        lambda e, z, prm: _species_rhs(e, z, prm, base_interp, props.Sc_thnf, p.beta),
        _species_bc, eta, z0, p=[gamma_guess], tol=tol, max_nodes=max_nodes, verbose=0)
    return sol
