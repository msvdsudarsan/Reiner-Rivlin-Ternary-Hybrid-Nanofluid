"""Base state AND eigenvalue problem on one Chebyshev grid.

The earlier trend computation interpolated a solve_bvp mesh with cubic splines to
build the perturbation coefficients.  Near the regularity boundary the base-state
curvature grows and that interpolation degrades, which is why three resolutions
stopped agreeing beyond lambda ~ -3.  Here the base state is itself solved by
Chebyshev collocation on the same grid the eigenvalue problem uses, so the
coefficients are spectrally accurate and no interpolation enters.
"""
import numpy as np
from cheb import chebD

def _grid(N, einf):
    D, xi = chebD(N)
    eta = einf*(1-xi)/2
    De = D*(-2/einf)
    return eta, De, De@De

def base_state_cheb(K, lam, S, N=120, einf=20.0, guess=None, tol=1e-12, itmax=60):
    """Newton solve of the momentum sub-system on a Chebyshev grid.
       Unknowns: F, G, H at the N+1 nodes.  Returns (eta, F, G, H, De, D2)."""
    eta, De, D2 = _grid(N, einf)
    n = N+1
    if guess is None:
        d = np.exp(-eta)
        F = lam*d; G = 1.0*d; H = -S + 0*eta
    else:
        F, G, H = guess
    I = np.eye(n)
    for it in range(itmax):
        Fp, Gp, Hp = De@F, De@G, De@H
        Fpp, Gpp = D2@F, D2@G
        Dl = 1 - 2*K*F
        r1 = Dl*Fpp - (F**2 - G**2 + H*Fp - K*(Fp**2 - Gp**2))
        r2 = Dl*Gpp - (2*F*G + H*Gp - 2*K*Fp*Gp)
        r3 = Hp + 2*F
        # Jacobian blocks
        J11 = np.diag(Dl)@D2 + np.diag(-2*K*Fpp - 2*F) - np.diag(H)@De + 2*K*np.diag(Fp)@De
        J12 = np.diag(2*G) - 2*K*np.diag(Gp)@De
        J13 = -np.diag(Fp)
        J21 = np.diag(-2*K*Gpp - 2*G) + 2*K*np.diag(Gp)@De
        J22 = np.diag(Dl)@D2 - np.diag(2*F) - np.diag(H)@De + 2*K*np.diag(Fp)@De
        J23 = -np.diag(Gp)
        J31 = 2*I; J32 = np.zeros((n,n)); J33 = De
        J = np.block([[J11,J12,J13],[J21,J22,J23],[J31,J32,J33]])
        r = np.concatenate([r1,r2,r3])
        # BCs.  chebD gives x[0]=+1, and eta = einf*(1-x)/2, so node 0 is the
        # wall (eta=0) and node N is the far field (eta=einf).
        w, f = 0, N
        for row, col, val, cur in [(w,      w,       lam,  F[w]),
                                   (n+w,    n+w,     1.0,  G[w]),
                                   (2*n+w,  2*n+w,   -S,   H[w]),
                                   (f,      f,       0.0,  F[f]),
                                   (n+f,    n+f,     0.0,  G[f])]:
            J[row,:] = 0.0
            J[row, col] = 1.0
            r[row] = cur - val
        try:
            dx = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            return None
        F = F + dx[:n]; G = G + dx[n:2*n]; H = H + dx[2*n:]
        if np.max(np.abs(dx)) < tol:
            return eta, F, G, H, De, D2
    return None


def base_state_from_bvp(K, lam, S, N=120, einf=20.0):
    """Initialise Newton from the solve_bvp continuation, then polish on the
       Chebyshev grid.  Newton from a cold start does not converge reliably
       near the boundary; seeded from a converged solution it converges in a
       handful of iterations."""
    from rr_solver import solve
    sol = solve(None, K, 0.0, S, einf=einf); l = 0.0
    step = 0.05 if lam > 0 else -0.05
    while abs(l - lam) > 1e-9:
        nxt = lam if abs(lam - l) < abs(step) else l + step
        s = solve(sol, K, nxt, S, einf=einf)
        if s.status == 0: sol, l = s, nxt
        else:
            step /= 2
            if abs(step) < 1e-7: return None
    eta, De, D2 = _grid(N, einf)
    y = sol.sol(eta)
    return base_state_cheb(K, lam, S, N=N, einf=einf, guess=(y[0], y[2], y[4]))


def gamma1_on_grid(K, eta, F, G, H, De, D2, gmax=1.2, ng=900):
    """Smallest momentum eigenvalue using base-state values ON the collocation
       grid -- no interpolation.  Same operator as cheb.gamma1_momentum."""
    N = len(eta) - 1; m = N + 1
    F0, G0, H0 = F, G, H
    F0p, G0p = De@F, De@G
    F0pp, G0pp = D2@F, D2@G
    Dl = 1 - 2*K*F0
    D2h = De.copy(); D2h[0, :] = 0; D2h[0, 0] = 1
    Z = np.eye(m); Z[0, 0] = 0
    P = -2.0*np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m, m))
    A_FF = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) - 2*K*np.diag(F0p)@De
            - np.diag(H0)@De + 2*np.diag(F0) - np.diag(F0p)@P)
    A_FG = 2*K*np.diag(G0p)@De - 2*np.diag(G0)
    A_GF = -4*K*np.diag(G0pp) + 2*K*np.diag(G0p)@De - 2*np.diag(G0) - np.diag(G0p)@P
    A_GG = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) + 2*K*np.diag(F0p)@De
            - np.diag(H0)@De - 2*np.diag(F0))
    A = np.block([[A_FF, A_FG], [A_GF, A_GG]]); B = np.block([[I, O], [O, I]])
    Af = -A.copy(); Bf = B.copy()
    for rr in (0, m, m+N):
        Af[rr, :] = 0; Af[rr, rr] = 1; Bf[rr, :] = 0
    Af[N, :] = 0; Af[N, 0:m] = De[N, :]; Bf[N, :] = 0
    gs = np.linspace(0.0, gmax, ng)
    sv = np.array([np.linalg.svd(Af - g*Bf, compute_uv=False)[-1] for g in gs])
    idx = [i for i in range(1, len(gs)-1) if sv[i] < sv[i-1] and sv[i] < sv[i+1]]
    return gs[idx[0]] if idx else None


def gamma_spectrum_direct(K, eta, F, G, H, De, D2, nkeep=6):
    """Direct generalized eigensolve  A q = gamma B q  (no scan, no local-minimum
       ambiguity).  Returns the smallest real, positive, finite eigenvalues."""
    from scipy.linalg import eig
    N = len(eta)-1; m = N+1
    F0p, G0p = De@F, De@G; F0pp, G0pp = D2@F, D2@G
    Dl = 1-2*K*F
    D2h = De.copy(); D2h[0,:] = 0; D2h[0,0] = 1
    Z = np.eye(m); Z[0,0] = 0
    P = -2.0*np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m,m))
    A_FF = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) - 2*K*np.diag(F0p)@De
            - np.diag(H)@De + 2*np.diag(F) - np.diag(F0p)@P)
    A_FG = 2*K*np.diag(G0p)@De - 2*np.diag(G)
    A_GF = -4*K*np.diag(G0pp) + 2*K*np.diag(G0p)@De - 2*np.diag(G) - np.diag(G0p)@P
    A_GG = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) + 2*K*np.diag(F0p)@De
            - np.diag(H)@De - 2*np.diag(F))
    A = -np.block([[A_FF,A_FG],[A_GF,A_GG]]); B = np.block([[I,O],[O,I]])
    for rr in (0, m, m+N):
        A[rr,:] = 0; A[rr,rr] = 1; B[rr,:] = 0
    A[N,:] = 0; A[N,0:m] = De[N,:]; B[N,:] = 0
    w = eig(A, B, right=False)
    w = w[np.isfinite(w)]
    w = w[np.abs(w.imag) < 1e-6*np.maximum(1,np.abs(w.real))].real
    w = np.sort(w[w > 1e-6])
    return w[:nkeep]


def smallest_real_eig(K, lam, S, N=100, einf=30.0):
    """True smallest real eigenvalue (negative allowed) with the physical
       far-field condition F(inf)=0.  Perturbations ~ exp(-gamma*tau)."""
    from scipy.linalg import eig
    r = base_state_from_bvp(K, lam, S, N=N, einf=einf)
    if r is None: return None
    eta, F, G, H, De, D2 = r; m = N+1
    F0p, G0p = De@F, De@G; F0pp, G0pp = D2@F, D2@G; Dl = 1-2*K*F
    D2h = De.copy(); D2h[0,:] = 0; D2h[0,0] = 1
    Z = np.eye(m); Z[0,0] = 0; P = -2.0*np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m,m))
    A_FF = np.diag(Dl)@D2+2*K*np.diag(F0pp)-2*K*np.diag(F0p)@De-np.diag(H)@De+2*np.diag(F)-np.diag(F0p)@P
    A_FG = 2*K*np.diag(G0p)@De-2*np.diag(G)
    A_GF = -4*K*np.diag(G0pp)+2*K*np.diag(G0p)@De-2*np.diag(G)-np.diag(G0p)@P
    A_GG = np.diag(Dl)@D2+2*K*np.diag(F0pp)+2*K*np.diag(F0p)@De-np.diag(H)@De-2*np.diag(F)
    A = -np.block([[A_FF,A_FG],[A_GF,A_GG]]); B = np.block([[I,O],[O,I]])
    for rr in (0, m, m+N, N):
        A[rr,:] = 0; B[rr,:] = 0; A[rr,rr] = 1.0
    w = eig(A, B, right=False); w = w[np.isfinite(w)]
    re = np.abs(w.imag) < 1e-6*np.maximum(1, np.abs(w.real))
    w = np.sort(w[re].real)
    w = w[np.abs(w) > 1e-6]          # drop the shrinking near-zero spurious pair
    return w[0] if len(w) else None
