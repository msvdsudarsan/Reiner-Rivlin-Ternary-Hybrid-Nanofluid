import numpy as np
from scipy.linalg import eig
from rr_solver import solve
from cheb import chebD

def base_state(K=0.3,lam=0.5,S=0.5,einf=25.0):
    sol=None;l=0.0;step=0.05
    while abs(l-lam)>1e-9:
        nxt=lam if abs(lam-l)<abs(step) else l+step
        s=solve(sol,K,nxt,S,einf=einf)
        if s.status==0: sol,l=s,nxt
        else: step/=2
    return sol

def species_spectrum(sol,Sc=883.1188,beta=0.2,N=160,zmax=60.0):
    """Solve in stretched wall coordinate zeta = Sc*eta.
       Phi_zz - H(zeta/Sc) Phi_z + Lam Phi = 0,  Phi(0)=Phi(zmax)=0
       gamma = Sc*beta + Sc^2 * Lam"""
    D,xi=chebD(N)
    z = zmax*(1-xi)/2
    Dz = D*(-2/zmax); D2z = Dz@Dz
    eta = z/Sc
    H = sol.sol(eta)[4]
    A = D2z - np.diag(H)@Dz
    # Dirichlet at both ends -> interior problem
    Ai = A[1:N,1:N]
    lam_vals = np.linalg.eigvals(-Ai)      # A Phi = -Lam Phi
    lam_vals = lam_vals[np.isfinite(lam_vals)]
    gam = Sc*beta + Sc**2*lam_vals
    gam = gam[np.abs(gam.imag)<1e-6*np.abs(gam.real)+1e-9].real
    return np.sort(gam)
