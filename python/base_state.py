"""
base_state.py

Python (scipy.integrate.solve_bvp) implementation of the base-state
boundary-value problem (main.tex Eqs. 12, 17, 18).

IMPORTANT NUMERICAL NOTE: an initial monolithic 9-state solve (F,G,H,
theta,phi together) was tried first and found NOT to be reliably
converged -- F'(0) shifted from -0.4430 to -0.4454 between tol=1e-9 and
tol=1e-10, and tightening further exhausted memory, because the large
Schmidt number (Sc~883, driven by a physically ordinary base-fluid mass
diffusivity D_f=1e-9 m^2/s against nu_thnf~8e-7 m^2/s) forces a very
thin concentration boundary layer that drags global mesh refinement
across ALL nine states even though only phi needs it.

This was fixed by exploiting Proposition 1 of main.tex directly: the
momentum sub-system (F,G,H) does not depend on theta or phi at all, so
it is solved on its own (5 states); theta and phi are then each solved
as their own SEPARATE 2-state problems using the already-known H(eta),
letting solve_bvp concentrate mesh refinement exactly where each field
needs it. This is both more numerically robust AND a direct
computational demonstration of Proposition 1's decoupling claim.

PYTHON IMPLEMENTATION -- NOT MATLAB R2026a.

ADDENDUM (found via MATLAB cross-check, see README_python_notes.md Bug
4): eta_inf=15 was subsequently found to be an insufficient domain size
for K=0.3 cold-start solves at these parameters -- solve_momentum_robust
below now defaults to eta_inf=40 (see parameters.py), confirmed stable
to 8 decimal places for eta_inf in {35,40,50}.
"""

import numpy as np
from scipy.integrate import solve_bvp
from dataclasses import dataclass

from parameters import Parameters
from effective_properties import EffectiveProperties


# ------------------------------------------------------------- momentum
def _mom_rhs(eta, y, K):
    F, Fp, G, Gp, H = y
    Delta = 1 - 2*K*F
    Delta = np.where(np.abs(Delta) < 1e-10, np.sign(Delta)*1e-10 + 1e-14, Delta)
    Fpp = (F**2 - G**2 + H*Fp - K*(Fp**2 - Gp**2)) / Delta
    Gpp = (2*F*G + H*Gp - 2*K*Fp*Gp) / Delta
    Hp = -2*F
    return np.vstack([Fp, Fpp, Gp, Gpp, Hp])


def _mom_bc(ya, yb, lam, S):
    return np.array([ya[0]-lam, ya[2]-1.0, ya[4]+S, yb[0], yb[2]])


def _mom_guess(eta, lam, S):
    decay = np.exp(-eta)
    F0 = lam*decay
    G0 = decay
    H0 = -S + lam*(1-decay)
    return np.vstack([F0, -lam*decay, G0, -decay, H0])


def solve_momentum(p: Parameters, eta_mesh=None, y_guess=None, tol=1e-10, max_nodes=200000):
    """Direct (single-shot) momentum solve at the target K. WARNING: this
    is NOT reliably unique -- see solve_momentum_robust below and the
    module docstring addendum for why this matters."""
    if eta_mesh is None:
        eta_mesh = np.linspace(0, p.eta_inf, p.Nmesh0)
    if y_guess is None:
        y_guess = _mom_guess(eta_mesh, p.lam, p.S)
    sol = solve_bvp(lambda e,y: _mom_rhs(e,y,p.K), lambda ya,yb: _mom_bc(ya,yb,p.lam,p.S),
                     eta_mesh, y_guess, tol=tol, max_nodes=max_nodes, verbose=0)
    return sol


def solve_momentum_robust(p: Parameters, tol=1e-10, max_nodes=300000):
    """Robust momentum solve via TWO-STAGE continuation, found necessary
    after a third round of cross-verification with MATLAB (see
    python/README_python_notes.md "Bug 5" for the full story).

    Stage 1 -- DOMAIN continuation at K=0: solve at a small, easily-
    converged eta_inf (8), then incrementally extend eta_inf (warm-
    starting each step from the previous, constant-tail-extended
    solution) up to the target eta_inf. This was found to be necessary
    because even the K=0 (classical Newtonian) problem at lambda=0.5,
    S=0.5 is NOT reliably obtained by a direct cold solve at a large
    eta_inf: different eta_inf values gave different "converged"
    (residual ~1e-10) answers (-0.287, -0.321, -0.313, ... depending on
    eta_inf), and a cold solve at eta_inf=60 failed outright. Only
    incremental domain growth from a small, well-behaved starting domain
    gives a value that is stable across eta_inf from 8 to 80: -0.287470.

    Stage 2 -- PARAMETER continuation in K: from that domain-correct
    K=0 solution, step K up to the target, warm-starting each step from
    the previous.

    Doing stage 2 first at the full target eta_inf directly (an earlier
    version of this function) is what caused MATLAB's bvp4c to require
    tens of thousands of mesh points and still fail to converge
    (residuals ~1e-2, not just slow) when this was cross-checked:
    jumping straight to a large domain with a K=0 cold-start guess is a
    poorly-conditioned problem regardless of solver implementation.
    """
    K_target = p.K

    # ---- stage 1: domain continuation at K=0 ----
    p_small = Parameters(K=0.0, lam=p.lam, S=p.S, eta_inf=8.0, Nmesh0=200,
                          RelTol=p.RelTol, AbsTol=p.AbsTol)
    mom = solve_momentum(p_small, tol=tol, max_nodes=max_nodes)
    if mom.status != 0:
        raise RuntimeError(f"K=0 small-domain solve failed: {mom.message}")

    domain_steps = sorted(set([d for d in [10, 15, 20, 25, 30, p.eta_inf] if 8 < d <= p.eta_inf] + [p.eta_inf]))

    for eta_inf in domain_steps:
        old_x = mom.x
        if eta_inf <= old_x[-1]:
            continue
        extra = np.linspace(old_x[-1], eta_inf, 100)[1:]
        new_x = np.concatenate([old_x, extra])
        tail = np.tile(mom.y[:, -1:], (1, extra.size))
        new_y = np.concatenate([mom.y, tail], axis=1)
        p_ext = Parameters(K=0.0, lam=p.lam, S=p.S, eta_inf=eta_inf,
                            RelTol=p.RelTol, AbsTol=p.AbsTol)
        mom = solve_momentum(p_ext, eta_mesh=new_x, y_guess=new_y, tol=tol, max_nodes=max_nodes)
        if mom.status != 0:
            raise RuntimeError(f"domain continuation failed extending to eta_inf={eta_inf}: {mom.message}")

    if K_target == 0.0:
        return mom

    # ---- stage 2: parameter continuation in K, at the fixed final domain ----
    n_K_steps = 31
    K_steps = np.linspace(0.0, K_target, n_K_steps)[1:]
    eta_mesh, y_guess = mom.x, mom.y
    for K in K_steps:
        pk = Parameters(K=K, lam=p.lam, S=p.S, eta_inf=p.eta_inf, Nmesh0=p.Nmesh0,
                          RelTol=p.RelTol, AbsTol=p.AbsTol)
        mom = solve_momentum(pk, eta_mesh=eta_mesh, y_guess=y_guess, tol=tol, max_nodes=max_nodes)
        if mom.status != 0:
            raise RuntimeError(f"K-continuation failed at K={K}: {mom.message}")
        eta_mesh, y_guess = mom.x, mom.y
    return mom


# --------------------------------------------------------------- energy
def _energy_rhs(eta, y, Pr, H_interp):
    th, thp = y
    H = H_interp(eta)
    thpp = Pr * H * thp
    return np.vstack([thp, thpp])


def _energy_bc(ya, yb):
    return np.array([ya[0]-1.0, yb[0]])


def solve_energy(mom_sol, Pr, eta_mesh=None, tol=1e-11, max_nodes=200000):
    if eta_mesh is None:
        eta_mesh = np.linspace(mom_sol.x[0], mom_sol.x[-1], 300)
    H_interp = lambda e: mom_sol.sol(e)[4]
    decay = np.exp(-eta_mesh)
    y0 = np.vstack([decay, -decay])
    sol = solve_bvp(lambda e,y: _energy_rhs(e,y,Pr,H_interp), _energy_bc,
                     eta_mesh, y0, tol=tol, max_nodes=max_nodes, verbose=0)
    return sol


# -------------------------------------------------------------- species
def _species_rhs(eta, y, Sc, beta, H_interp):
    ph, php = y
    H = H_interp(eta)
    phpp = Sc*H*php + Sc*beta*ph
    return np.vstack([php, phpp])


def _species_bc(ya, yb):
    return np.array([ya[0]-1.0, yb[0]])


def solve_species(mom_sol, Sc, beta, eta_mesh=None, tol=1e-11, max_nodes=400000):
    if eta_mesh is None:
        # concentrate initial mesh near the wall for large Sc (thin layer)
        eta_lin = np.linspace(mom_sol.x[0], mom_sol.x[-1], 200)
        eta_mesh = np.unique(np.concatenate([eta_lin, np.geomspace(1e-4, max(mom_sol.x[-1],1e-3), 400)]))
        eta_mesh = eta_mesh[eta_mesh <= mom_sol.x[-1]]
    H_interp = lambda e: mom_sol.sol(e)[4]
    decay = np.exp(-eta_mesh)
    y0 = np.vstack([decay, -decay])
    sol = solve_bvp(lambda e,y: _species_rhs(e,y,Sc,beta,H_interp), _species_bc,
                     eta_mesh, y0, tol=tol, max_nodes=max_nodes, verbose=0)
    return sol


@dataclass
class BaseState:
    mom_sol: object
    energy_sol: object
    species_sol: object

    @property
    def Fp0(self): return float(self.mom_sol.y[1,0])
    @property
    def Gp0(self): return float(self.mom_sol.y[3,0])
    @property
    def Hinf(self): return float(self.mom_sol.y[4,-1])
    @property
    def thp0(self): return float(self.energy_sol.y[1,0])
    @property
    def php0(self): return float(self.species_sol.y[1,0])


def solve_base_state_decoupled(p: Parameters, props: EffectiveProperties, tol_mom=1e-10, tol_scalar=1e-11, robust=True):
    mom_sol = solve_momentum_robust(p, tol=tol_mom) if robust else solve_momentum(p, tol=tol_mom)
    if mom_sol.status != 0:
        raise RuntimeError(f"momentum solve failed: {mom_sol.message}")
    energy_sol = solve_energy(mom_sol, props.Pr_thnf, tol=tol_scalar)
    if energy_sol.status != 0:
        raise RuntimeError(f"energy solve failed: {energy_sol.message}")
    species_sol = solve_species(mom_sol, props.Sc_thnf, p.beta, tol=tol_scalar)
    if species_sol.status != 0:
        raise RuntimeError(f"species solve failed: {species_sol.message}")
    return BaseState(mom_sol, energy_sol, species_sol)
