"""Perturbation eigenproblem for the similarity-preserving momentum modes.

Boundary conditions follow the Supplementary Material (S5): f(0)=g(0)=0, g(inf)=0 and the
relaxed far-field condition f'(inf)=0, imposed with the differentiation-matrix row at the
last node.  The whole real spectrum is kept (negative values included).  Eigenvalues are
returned with their relative eigenpair residual and the share of the eigenfunction amplitude
lying in the outer half of the domain, so that localised modes can be told apart from
modes tied to the truncated far field.

Convention: perturbations vary as exp(-gamma*tau); gamma > 0 decays, gamma < 0 grows.
"""
import numpy as np
from scipy.linalg import eig
from cheb_base import base_state_from_bvp


def pencil(K, eta, F, G, H, De, D2):
    N = len(eta) - 1; m = N + 1
    F0p, G0p = De @ F, De @ G
    F0pp, G0pp = D2 @ F, D2 @ G
    Dl = 1 - 2*K*F
    D2h = De.copy(); D2h[0, :] = 0; D2h[0, 0] = 1
    Z = np.eye(m); Z[0, 0] = 0
    P = -2.0*np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m, m))
    AFF = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) - 2*K*np.diag(F0p)@De
           - np.diag(H)@De + 2*np.diag(F) - np.diag(F0p)@P)
    AFG = 2*K*np.diag(G0p)@De - 2*np.diag(G)
    AGF = -4*K*np.diag(G0pp) + 2*K*np.diag(G0p)@De - 2*np.diag(G) - np.diag(G0p)@P
    AGG = (np.diag(Dl)@D2 + 2*K*np.diag(F0pp) + 2*K*np.diag(F0p)@De
           - np.diag(H)@De - 2*np.diag(F))
    A = -np.block([[AFF, AFG], [AGF, AGG]]); B = np.block([[I, O], [O, I]])
    for r in (0, m, m + N):          # f(0)=0, g(0)=0, g(inf)=0
        A[r, :] = 0; A[r, r] = 1; B[r, :] = 0
    A[N, :] = 0; A[N, 0:m] = De[N, :]; B[N, :] = 0      # f'(inf)=0
    return A, B


def real_modes(K, lam, S, N=100, einf=30.0, gmax=60.0):
    """All real eigenpairs with |gamma| < gmax: list of (gamma, vector, residual, tail)."""
    r = base_state_from_bvp(K, lam, S, N=N, einf=einf)
    if r is None:
        return None, None
    eta, F, G, H, De, D2 = r
    A, B = pencil(K, eta, F, G, H, De, D2)
    w, V = eig(A, B)
    ok = np.isfinite(w) & (np.abs(w.imag) < 1e-8*np.maximum(1, np.abs(w.real)))
    w = w[ok].real; V = V[:, ok].real
    n = len(eta); out = []
    for i in range(len(w)):
        if abs(w[i]) > gmax:
            continue
        v = V[:, i]/np.linalg.norm(V[:, i])
        res = np.linalg.norm(A@v - w[i]*(B@v))/(np.linalg.norm(A@v) + abs(w[i]))
        amp = v[:n]**2 + v[n:]**2
        out.append((w[i], v, res, amp[eta >= 0.5*einf].sum()/amp.sum()))
    return eta, sorted(out, key=lambda o: o[0])


def smallest_reliable(K, lam, S=1.5, N=100, einf=30.0, rtol=1e-4):
    """Smallest real eigenvalue that passes the residual test.  The exactly-zero value that
    the singular pencil always carries is excluded (|gamma| < 1e-6)."""
    eta, out = real_modes(K, lam, S, N, einf)
    if out is None:
        return None
    c = [o for o in out if abs(o[0]) > 1e-6 and o[2] < rtol]
    return c[0] if c else None


def crossing(K, S, N, einf, lo, hi, iters=40):
    """Zero crossing of the smallest reliable eigenvalue in [hi, lo] (lo > hi), by bisection."""
    glo = smallest_reliable(K, lo, S, N, einf)[0]
    for _ in range(iters):
        mid = 0.5*(lo + hi)
        g = smallest_reliable(K, mid, S, N, einf)
        if g is None:
            break
        if (g[0] > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
        if abs(hi - lo) < 1e-7:
            break
    return 0.5*(lo + hi)
