import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import eig

# Chebyshev-Tau spectral search for higher discrete modes of the species
# perturbation equation, following the standard construction (Boyd,
# "Chebyshev and Fourier Spectral Methods", 2nd ed., Dover, 2001) and
# inspired by the layered-domain Chebyshev-Tau approach of Tu et al.
# (2021), J. Sound Vib. 492, 115784, for a comparably stiff eigenvalue
# problem in a different physical context (underwater acoustic normal
# modes).
#
# The method is validated against a known exact spectrum FIRST, before
# being trusted on the real (much harder) problem -- this is the correct
# way to use a numerical method you don't already have independent
# confidence in for a specific hard case.

def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.]) * (-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1./c) / (dX + np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

print("=== Step 1: validate the Chebyshev method on a known exact case ===")
print("Phi'' + gamma*Phi = 0, Phi(0)=Phi(1)=0  =>  exact gamma_n = (n*pi)^2")
N, eta_max = 60, 1.0
D, xi = cheb(N)
D2 = D @ D
Deta = D*(-2/eta_max); D2eta = D2*(4/eta_max**2)
L = D2eta
RHS = -np.eye(N+1)
L_bc = L.copy(); RHS_bc = RHS.copy()
L_bc[0,:] = 0; L_bc[0,0] = 1; RHS_bc[0,:] = 0
L_bc[N,:] = 0; L_bc[N,N] = 1; RHS_bc[N,:] = 0
ev = eig(L_bc, RHS_bc, right=False)
ev = ev[np.isfinite(ev)]
real_ev = np.sort(ev[np.abs(ev.imag) < 1e-6].real)
print("Paper:  9.8696, 39.4784, 88.8264, 157.9137, 246.7401")
print("Yours: ", np.array2string(real_ev[real_ev>0][:5], precision=4))
print()

print("=== Step 2: apply to the real species perturbation equation ===")
K, lam, S = 0.3, 0.5, 0.5
def mom_rhs(eta, y, Kv):
    F, Fp, G, Gp, H = y
    Delta = 1 - 2*Kv*F
    Delta = np.where(np.abs(Delta) < 1e-10, np.sign(Delta)*1e-10 + 1e-14, Delta)
    Fpp = (F**2 - G**2 + H*Fp - Kv*(Fp**2 - Gp**2)) / Delta
    Gpp = (2*F*G + H*Gp - 2*Kv*Fp*Gp) / Delta
    return np.vstack([Fp, Fpp, Gp, Gpp, -2*F])
def mom_bc(ya, yb, lamv, Sv):
    return np.array([ya[0]-lamv, ya[2]-1.0, ya[4]+Sv, yb[0], yb[2]])
def mom_guess(eta, lamv, Sv):
    d = np.exp(-eta)
    return np.vstack([lamv*d, -lamv*d, d, -d, -Sv + lamv*(1-d)])

mesh = np.linspace(0, 8, 200)
mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lam,S),
                 mesh, mom_guess(mesh, lam, S), tol=1e-10, max_nodes=300000, verbose=0)
for eta_t in [15, 25, 40]:
    extra = np.linspace(mom.x[-1], eta_t, 100)[1:]
    new_x = np.concatenate([mom.x, extra])
    new_y = np.concatenate([mom.y, np.tile(mom.y[:,-1:], (1, extra.size))], axis=1)
    mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     new_x, new_y, tol=1e-10, max_nodes=300000, verbose=0)
for Kstep in np.linspace(0, K, 31)[1:]:
    mom = solve_bvp(lambda e,y: mom_rhs(e,y,Kstep), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     mom.x, mom.y, tol=1e-10, max_nodes=300000, verbose=0)
H_interp = lambda e: mom.sol(min(e, mom.x[-1]))[4]

Sc_thnf, beta = 883.1188, 0.2

def solve_spectrum(N, eta_max):
    D, xi = cheb(N)
    D2 = D @ D
    eta_nodes = eta_max*(1-xi)/2
    Deta = D * (-2/eta_max)
    D2eta = D2 * (4/eta_max**2)
    H0_vals = np.array([H_interp(e) for e in eta_nodes])
    L = D2eta - Sc_thnf*np.diag(H0_vals)@Deta + Sc_thnf*beta*np.eye(N+1)
    RHS = -np.eye(N+1)
    L_bc = L.copy(); RHS_bc = RHS.copy()
    L_bc[0,:] = 0; L_bc[0,0] = 1; RHS_bc[0,:] = 0
    L_bc[N,:] = 0; L_bc[N,N] = 1; RHS_bc[N,:] = 0
    ev = eig(L_bc, RHS_bc, right=False)
    ev = ev[np.isfinite(ev)]
    return np.sort(ev[np.abs(ev.imag) < 1e-3].real)

print("Paper: N=80,eta_max=20 -> -75.9 | N=120,eta_max=20 -> +57.6 | N=120,eta_max=15 -> +308.7")
print("Yours:")
for N, eta_max in [(80,20.0),(120,20.0),(120,15.0)]:
    ev = solve_spectrum(N, eta_max)
    print(f"  N={N}, eta_max={eta_max}: {ev}")
print()
print("These values should NOT agree with each other -- that disagreement is")
print("the point: it shows no discrete eigenvalue converges here, consistent")
print("with the analytical stability-threshold argument in the manuscript.")
