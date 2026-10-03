"""Linear stability of the shrinking branch: the corrected analysis.

Earlier versions located the smallest momentum eigenvalue by scanning the
smallest singular value of the collocation pencil over gamma in [0, 1.2].  A
one-sided scan over non-negative gamma cannot detect a negative eigenvalue, and
it missed one.  This script does the analysis two ways that do not share that
blind spot:

  1. the base state is solved by Chebyshev collocation on the SAME grid as the
     eigenvalue problem (no interpolation of a solve_bvp mesh), and the
     spectrum is obtained by a direct generalised eigensolve;
  2. independently of any eigensolver, the linearised equations are integrated
     forward in time from a random perturbation and the growth rate measured.

Expected output (K=+0.3, S=1.5):
    lambda* = -0.4469          (zero crossing, bisected)
    gamma_1(-1.0) = -0.5886     eigensolve
    growth rate   = +0.5974     time integration; equals -gamma_1 once the
                                implicit-Euler damping ln(1/(1-r dt))/dt is applied
"""
import numpy as np
from cheb_base import base_state_from_bvp, smallest_real_eig

def bisect_lambda_star(K, S=1.5, a=-0.25, b=-0.8, iters=14):
    for _ in range(iters):
        m = 0.5*(a+b)
        if smallest_real_eig(K, m, S) > 0: a = m
        else: b = m
    return 0.5*(a+b)

def growth_rate(K, lam, S=1.5, N=80, einf=25.0, dt=0.05, steps=400):
    eta, F, G, H, De, D2 = base_state_from_bvp(K, lam, S, N=N, einf=einf)
    m = N+1
    F0p, G0p = De@F, De@G; F0pp, G0pp = D2@F, D2@G; Dl = 1-2*K*F
    D2h = De.copy(); D2h[0,:] = 0; D2h[0,0] = 1
    Z = np.eye(m); Z[0,0] = 0; P = -2.0*np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m,m))
    AFF = np.diag(Dl)@D2+2*K*np.diag(F0pp)-2*K*np.diag(F0p)@De-np.diag(H)@De+2*np.diag(F)-np.diag(F0p)@P
    AFG = 2*K*np.diag(G0p)@De-2*np.diag(G)
    AGF = -4*K*np.diag(G0pp)+2*K*np.diag(G0p)@De-2*np.diag(G)-np.diag(G0p)@P
    AGG = np.diag(Dl)@D2+2*K*np.diag(F0pp)+2*K*np.diag(F0p)@De-np.diag(H)@De-2*np.diag(F)
    A = -np.block([[AFF,AFG],[AGF,AGG]]); B = np.block([[I,O],[O,I]])
    bc = [0, m, m+N, N]
    for r in bc: A[r,:] = 0; B[r,:] = 0; A[r,r] = 1.0
    q = np.zeros(2*m); q[:m] = np.sin(np.pi*eta/einf)*np.exp(-eta/3)
    Mi = np.linalg.inv(B + dt*A)
    ts, ns = [], []
    for n in range(1, steps+1):
        q = Mi@(B@q); q[bc] = 0
        if n % 40 == 0: ts.append(n*dt); ns.append(np.linalg.norm(q))
    return np.polyfit(ts[-5:], np.log(ns[-5:]), 1)[0]

if __name__ == "__main__":
    for K in (0.3, -0.3):
        print("K=%+.1f, S=1.5:  lambda* = %.4f" % (K, bisect_lambda_star(K)))
    g = smallest_real_eig(0.3, -1.0, 1.5)
    r = growth_rate(0.3, -1.0)
    print("\nK=+0.3, lambda=-1.0")
    print("  eigensolve      gamma_1 = %+.4f" % g)
    print("  time integration  rate  = %+.4f" % r)
    print("  predicted rate from eigensolve after implicit-Euler damping = %+.4f"
          % (np.log(1/(1 - (-g)*0.05))/0.05))
