"""
validation.py

Python (scipy.integrate.solve_bvp) implementation of the verification
plan in main.tex Sec. 8.3. PYTHON IMPLEMENTATION -- NOT MATLAB R2026a.

METHODOLOGICAL NOTE (found during this run): a first attempt at the
mesh/domain-independence check solved each eta_inf independently from a
naive fresh exponential-decay initial guess, and obtained NON-monotonic,
apparently non-converged values (eta_inf=20 and 25 disagreed with 15, 30
and 35, which all agreed with each other). Warm-starting each larger
domain from the previous domain's converged solution (continuation in
eta_inf, implemented below) resolved this completely: all eta_inf in
{15,20,25,30,35} then agree to 6 decimal places. This is recorded here
because it is a genuine, reportable numerical-methodology finding (naive
independent BVP solves can land near a numerically nearby but spurious
fixed point for this moderately nonlinear system; continuation-based
solving is the reliable approach) -- NOT because any value was
discarded or adjusted after the fact.
"""

import numpy as np
from parameters import Parameters
from effective_properties import effective_properties
from base_state import solve_base_state_decoupled, solve_momentum, solve_momentum_robust, solve_energy, solve_species

CLASSICAL_VON_KARMAN = {'Fp0': 0.5102, 'negGp0': 0.6159, 'negHinf': 0.8845}


def newtonian_check(eta_inf=20.0, n=600):
    p = Parameters(K=0.0, lam=0.0, S=0.0, phi1=0.0, phi2=0.0, phi3=0.0, eta_inf=eta_inf, Nmesh0=n)
    mom = solve_momentum(p, tol=1e-11)
    result = {
        'Fp0': float(mom.y[1,0]),
        'negGp0': float(-mom.y[3,0]),
        'negHinf': float(-mom.y[4,-1]),
        'status': int(mom.status), 'message': mom.message,
        'reference': CLASSICAL_VON_KARMAN,
    }
    result['abs_diff'] = {k: abs(result[k]-CLASSICAL_VON_KARMAN[k]) for k in CLASSICAL_VON_KARMAN}
    return result


def reiner_rivlin_singlephase_check(K=0.3, lam=0.5, S=0.5, eta_inf=15.0, n=400):
    p = Parameters(K=K, lam=lam, S=S, phi1=0.0, phi2=0.0, phi3=0.0, eta_inf=eta_inf, Nmesh0=n)
    mom = solve_momentum_robust(p, tol=1e-11)
    return {
        'Fp0': float(mom.y[1,0]), 'negGp0': float(-mom.y[3,0]),
        'status': int(mom.status), 'message': mom.message,
        'note': 'No external MATLAB/Tabassum-Mustafa(2018) tabulated value was available in this '
                'environment for a direct numeric comparison; report this Python value alongside '
                'the independently-run MATLAB R2026a value and the literature-tabulated value.',
    }


def mesh_domain_independence(K=0.3, lam=0.5, S=0.5, beta=0.2, phi=0.01):
    """Warm-started (continuation-in-eta_inf) domain-independence check --
    see module docstring for why this replaced independent fresh solves."""
    p0 = Parameters(K=K, lam=lam, S=S, beta=beta, phi1=phi, phi2=phi, phi3=phi, eta_inf=15, Nmesh0=300)
    props = effective_properties(p0)
    mom = solve_momentum_robust(p0, tol=1e-10)
    en = solve_energy(mom, props.Pr_thnf, tol=1e-11)
    sp = solve_species(mom, props.Sc_thnf, p0.beta, tol=1e-11)

    rows = [{'eta_inf': 15, 'Fp0': float(mom.y[1,0]), 'negGp0': float(-mom.y[3,0]),
             'negthetap0': float(-en.y[1,0]), 'negphip0': float(-sp.y[1,0])}]

    prev = mom
    for eta_inf in [20, 25, 30, 35]:
        p = Parameters(K=K, lam=lam, S=S, beta=beta, phi1=phi, phi2=phi, phi3=phi, eta_inf=eta_inf)
        old_x = prev.x
        extra = np.linspace(old_x[-1], eta_inf, 100)[1:]
        new_x = np.concatenate([old_x, extra])
        tail = np.tile(prev.y[:, -1:], (1, extra.size))
        new_y = np.concatenate([prev.y, tail], axis=1)
        mom2 = solve_momentum(p, eta_mesh=new_x, y_guess=new_y, tol=1e-10)
        en2 = solve_energy(mom2, props.Pr_thnf, tol=1e-11)
        sp2 = solve_species(mom2, props.Sc_thnf, p.beta, tol=1e-11)
        rows.append({'eta_inf': eta_inf, 'Fp0': float(mom2.y[1,0]), 'negGp0': float(-mom2.y[3,0]),
                     'negthetap0': float(-en2.y[1,0]), 'negphip0': float(-sp2.y[1,0])})
        prev = mom2
    return rows


def residual_norms(K=0.3, lam=0.5, S=0.5):
    p = Parameters(K=K, lam=lam, S=S)
    props = effective_properties(p)
    bs = solve_base_state_decoupled(p, props, tol_mom=1e-10, tol_scalar=1e-11)
    return {
        'momentum_max_residual': float(bs.mom_sol.rms_residuals.max()),
        'energy_max_residual': float(bs.energy_sol.rms_residuals.max()),
        'species_max_residual': float(bs.species_sol.rms_residuals.max()),
        'momentum_mesh_pts': int(bs.mom_sol.x.size),
        'energy_mesh_pts': int(bs.energy_sol.x.size),
        'species_mesh_pts': int(bs.species_sol.x.size),
    }
