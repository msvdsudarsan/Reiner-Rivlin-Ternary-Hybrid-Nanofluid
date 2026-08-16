"""
stability_shooting.py

A more robust shooting-method implementation of the momentum stability
eigenvalue problem, replacing the scipy.integrate.solve_bvp-with-unknown-
-parameter approach in stability.py, which was found (during this run)
to fail to converge (mesh-node exhaustion) even for the well-behaved
Example 1 baseline case, regardless of initial guess quality.

Root cause understood: solve_bvp's Newton iteration over the unknown
parameter gamma was poorly conditioned for this problem; a direct 2D
shooting/root-find (integrate the linear homogeneous ODE forward with
solve_ivp for trial (alpha, gamma), require the two far-field decay
conditions Ff(eta_max)=0, Gf(eta_max)=0) is far more robust for a LINEAR
homogeneous eigenvalue problem of this kind, since it does not require
solving a nonlinear collocation Jacobian with the eigenvalue as an
unknown parameter simultaneously with the state.

State for the IVP: [Ff, Ffp, Gf, Gfp, Hf], integrated forward from
Ff(0)=0, Ff'(0)=1 (normalization), Gf(0)=0, Gf'(0)=alpha (unknown, to be
shot), Hf(0)=0. For fixed gamma, alpha is chosen (1D root-find) so that
Ff(eta_max) is small; the SECOND far-field condition Gf(eta_max)=0 is
then used as the residual whose root in gamma gives the eigenvalue.

PYTHON IMPLEMENTATION -- not part of any other computational environment.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def _rhs(eta, z, gamma, K, base_interp):
    y0 = base_interp(eta)
    F0, F0p, G0, G0p, H0 = y0[0], y0[1], y0[2], y0[3], y0[4]
    Delta = 1 - 2*K*F0
    if abs(Delta) < 1e-10:
        Delta = 1e-10 if Delta >= 0 else -1e-10
    F0pp = (F0**2 - G0**2 + H0*F0p - K*(F0p**2 - G0p**2)) / Delta
    G0pp = (2*F0*G0 + H0*G0p - 2*K*F0p*G0p) / Delta

    Ff, Ffp, Gf, Gfp, Hf = z
    Ffpp = (2*K*Ff*F0pp - 2*K*F0p*Ffp + 2*K*G0p*Gfp - gamma*Ff
            + 2*F0*Ff - 2*G0*Gf + H0*Ffp + Hf*F0p) / Delta
    Gfpp = (2*K*Ff*G0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p - gamma*Gf
            + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p) / Delta
    Hfp = -2*Ff
    return [Ffp, Ffpp, Gfp, Gfpp, Hfp]


def _integrate(gamma, alpha, K, base_interp, eta_max):
    z0 = [0.0, 1.0, 0.0, alpha, 0.0]
    sol = solve_ivp(_rhs, [0, eta_max], z0, args=(gamma, K, base_interp),
                     method='RK45', rtol=1e-8, atol=1e-10, dense_output=False,
                     max_step=eta_max/50)
    return sol


def _far_field_residuals(gamma, alpha, K, base_interp, eta_max):
    sol = _integrate(gamma, alpha, K, base_interp, eta_max)
    if not sol.success or np.any(np.abs(sol.y) > 1e8):
        return 1e8, 1e8
    return sol.y[0, -1], sol.y[2, -1]  # Ff(eta_max), Gf(eta_max)


def solve_for_alpha(gamma, K, base_interp, eta_max, alpha_bracket=(-5, 5), n_scan=41):
    """For fixed gamma, find alpha such that Ff(eta_max) = 0 (root in alpha)."""
    alphas = np.linspace(*alpha_bracket, n_scan)
    vals = []
    for a in alphas:
        Ff_end, _ = _far_field_residuals(gamma, a, K, base_interp, eta_max)
        vals.append(Ff_end)
    vals = np.array(vals)
    sign_changes = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if len(sign_changes) == 0:
        return None
    i = sign_changes[0]
    try:
        alpha_root = brentq(lambda a: _far_field_residuals(gamma, a, K, base_interp, eta_max)[0],
                             alphas[i], alphas[i+1], xtol=1e-10)
        return alpha_root
    except Exception:
        return None


def Gf_end_residual(gamma, K, base_interp, eta_max, alpha_bracket=(-5, 5)):
    """Residual function whose root in gamma is the momentum eigenvalue:
    for the alpha that zeros Ff(eta_max), report Gf(eta_max)."""
    alpha = solve_for_alpha(gamma, K, base_interp, eta_max, alpha_bracket)
    if alpha is None:
        return None
    _, Gf_end = _far_field_residuals(gamma, alpha, K, base_interp, eta_max)
    return Gf_end


def find_momentum_eigenvalue(base_sol, K, eta_max=None, gamma_bracket=(-3, 3), n_scan=25, alpha_bracket=(-5,5)):
    if eta_max is None:
        eta_max = base_sol.x[-1]
    base_interp = base_sol.sol

    gammas = np.linspace(*gamma_bracket, n_scan)
    residuals = []
    for g in gammas:
        r = Gf_end_residual(g, K, base_interp, eta_max, alpha_bracket)
        residuals.append(r if r is not None else np.nan)
    residuals = np.array(residuals, dtype=float)

    valid = ~np.isnan(residuals)
    sign_changes = []
    idxs = np.where(valid)[0]
    for j in range(len(idxs)-1):
        i0, i1 = idxs[j], idxs[j+1]
        if i1 == i0+1 and np.sign(residuals[i0]) != np.sign(residuals[i1]):
            sign_changes.append(i0)

    roots = []
    for i in sign_changes:
        try:
            root = brentq(lambda g: Gf_end_residual(g, K, base_interp, eta_max, alpha_bracket) or 1e8,
                           gammas[i], gammas[i+1], xtol=1e-8)
            roots.append(root)
        except Exception:
            pass
    return roots, gammas, residuals


# ---------------------------------------------------------- scalar (thermal/species)
def _scalar_rhs(eta, z, gamma, coeff_H, coeff_reaction, base_interp):
    """Generic RHS for Theta''=Pr*H0*Theta'-gamma*Theta (coeff_reaction=0)
    or Phi''=Sc*H0*Phi'+Sc*beta*Phi-gamma*Phi (coeff_reaction=Sc*beta)."""
    y0 = base_interp(eta)
    H0 = y0[4]
    z1, z2 = z
    z2p = coeff_H*H0*z2 + coeff_reaction*z1 - gamma*z1
    return [z2, z2p]


def _scalar_end_value(gamma, coeff_H, coeff_reaction, base_interp, eta_max):
    z0 = [0.0, 1.0]
    sol = solve_ivp(_scalar_rhs, [0, eta_max], z0, args=(gamma, coeff_H, coeff_reaction, base_interp),
                     method='RK45', rtol=1e-9, atol=1e-11, max_step=eta_max/50)
    if not sol.success or np.any(np.abs(sol.y) > 1e10):
        return None
    return sol.y[0, -1]


def find_scalar_eigenvalue(base_sol, coeff_H, coeff_reaction, eta_max=None,
                            gamma_bracket=(0.01, 30), n_scan=40):
    """Find the smallest positive eigenvalue of a scalar (thermal or
    species) perturbation block by scanning gamma and locating sign
    changes in Theta/Phi(eta_max), then refining with brentq."""
    if eta_max is None:
        eta_max = base_sol.x[-1]
    base_interp = base_sol.sol

    gammas = np.linspace(*gamma_bracket, n_scan)
    vals = [_scalar_end_value(g, coeff_H, coeff_reaction, base_interp, eta_max) for g in gammas]
    vals = np.array([v if v is not None else np.nan for v in vals])

    roots = []
    for i in range(len(gammas)-1):
        if np.isnan(vals[i]) or np.isnan(vals[i+1]):
            continue
        if np.sign(vals[i]) != np.sign(vals[i+1]):
            try:
                root = brentq(lambda g: _scalar_end_value(g, coeff_H, coeff_reaction, base_interp, eta_max) or 1e8,
                               gammas[i], gammas[i+1], xtol=1e-9)
                roots.append(root)
            except Exception:
                pass
    return roots, gammas, vals
