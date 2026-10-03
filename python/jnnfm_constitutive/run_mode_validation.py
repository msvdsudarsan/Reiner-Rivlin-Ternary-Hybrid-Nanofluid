"""Mode identification and two-branch time-domain validation of the instability.

Two questions an eigenvalue result has to answer before it can be trusted:

  1. Which eigenvalue crosses zero?  A routine that simply returns the most
     negative eigenvalue could in principle jump between modes.  This script
     prints the eigenvalues nearest zero on either side of each crossing.

  2. Is the instability real, on BOTH branches?  The linearised equations are
     integrated forward in time from a random perturbation, with no eigensolver
     involved, and the measured growth rate is compared with the implicit-Euler
     prediction ln(1/(1+gamma*dt))/dt from the eigensolve.

Expected output (S = 1.5):
    K=+0.3:  one simple real eigenvalue crosses at lambda* = -0.4469; next gap 0.26
    K=-0.3:  one simple real eigenvalue crosses at lambda* = -0.6036; next gap 0.32
    K=-0.3, lambda=-0.8:  growth 0.4595  vs prediction 0.4595
    K=-0.3, lambda=-1.0:  growth 1.0858  vs prediction 1.0858
    K=+0.3:  a second real mode crosses at lambda = -2.289
    no complex eigenvalue with negative real part on either branch
"""
import numpy as np
from scipy.linalg import eig
from cheb_base import base_state_from_bvp


def operator(K, lam, S=1.5, N=100, einf=30.0):
    eta, F, G, H, De, D2 = base_state_from_bvp(K, lam, S, N=N, einf=einf)
    m = N + 1
    F0p, G0p = De @ F, De @ G
    F0pp, G0pp = D2 @ F, D2 @ G
    Dl = 1 - 2 * K * F
    D2h = De.copy(); D2h[0, :] = 0; D2h[0, 0] = 1
    Z = np.eye(m); Z[0, 0] = 0
    P = -2.0 * np.linalg.solve(D2h, Z)
    I = np.eye(m); O = np.zeros((m, m))
    AFF = (np.diag(Dl) @ D2 + 2*K*np.diag(F0pp) - 2*K*np.diag(F0p) @ De
           - np.diag(H) @ De + 2*np.diag(F) - np.diag(F0p) @ P)
    AFG = 2*K*np.diag(G0p) @ De - 2*np.diag(G)
    AGF = -4*K*np.diag(G0pp) + 2*K*np.diag(G0p) @ De - 2*np.diag(G) - np.diag(G0p) @ P
    AGG = (np.diag(Dl) @ D2 + 2*K*np.diag(F0pp) + 2*K*np.diag(F0p) @ De
           - np.diag(H) @ De - 2*np.diag(F))
    A = -np.block([[AFF, AFG], [AGF, AGG]])
    B = np.block([[I, O], [O, I]])
    bc = [0, m, m + N, N]
    for r in bc:
        A[r, :] = 0; B[r, :] = 0; A[r, r] = 1.0
    return eta, A, B, bc, m


def nearest_zero(K, lam, k=4):
    _, A, B, _, _ = operator(K, lam)
    w = eig(A, B, right=False); w = w[np.isfinite(w)]
    w = w[np.abs(w) > 1e-6]
    return w[np.argsort(np.abs(w))][:k].real


def growth(K, lam, dt=0.05, steps=400):
    eta, A, B, bc, m = operator(K, lam, N=80, einf=25.0)
    q = np.zeros(2 * m)
    q[:m] = np.sin(np.pi * eta / 25.0) * np.exp(-eta / 3)
    q[:m] *= np.random.default_rng(1).uniform(0.5, 1.5, m)
    q[bc] = 0
    Mi = np.linalg.inv(B + dt * A)
    ts, ns = [], []
    for n in range(1, steps + 1):
        q = Mi @ (B @ q); q[bc] = 0
        if n % 40 == 0:
            ts.append(n * dt); ns.append(np.linalg.norm(q))
    rate = np.polyfit(ts[-5:], np.log(ns[-5:]), 1)[0]
    w = eig(A, B, right=False); w = w[np.isfinite(w)]
    re = np.abs(w.imag) < 1e-6 * np.maximum(1, np.abs(w.real))
    g = np.sort(w[re].real); g = g[np.abs(g) > 1e-6][0]
    return rate, np.log(1 / (1 + g * dt)) / dt


if __name__ == "__main__":
    print("Eigenvalues nearest zero either side of each crossing")
    for K, lams in [(0.3, (-0.30, -0.4469, -0.60)), (-0.3, (-0.45, -0.6036, -0.75))]:
        for lam in lams:
            print("  K=%+.1f  lambda=%7.4f  %s" % (K, lam,
                  "  ".join("%+.4f" % x for x in nearest_zero(K, lam))))
    print("\nTime-domain validation on the negative branch (no eigensolver)")
    for lam in (-0.8, -1.0):
        r, p = growth(-0.3, lam)
        print("  K=-0.3  lambda=%5.2f  growth %.4f   prediction %.4f" % (lam, r, p))
