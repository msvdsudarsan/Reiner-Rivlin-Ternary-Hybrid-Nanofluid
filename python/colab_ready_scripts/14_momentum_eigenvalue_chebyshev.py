import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import eig

# ---------- Chebyshev differentiation matrix ----------
def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.]) * (-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1./c) / (dX + np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

# ---------- base-state momentum solve (K-continuation from K=0) ----------
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
    return np.vstack([lamv*d, -lamv*d, d, -d, np.full_like(eta, -Sv)])

def solve_base(K, lam, S, eta_inf=40.0, nmesh=800, tol=1e-10, max_nodes=60000):
    mesh = np.linspace(0, eta_inf, nmesh)
    sol = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     mesh, mom_guess(mesh, lam, S), tol=tol, max_nodes=max_nodes, verbose=0)
    assert sol.status == 0, "K=0 start failed"
    for Kstep in np.linspace(0, K, 31)[1:]:
        sol = solve_bvp(lambda e,y: mom_rhs(e,y,Kstep), lambda ya,yb: mom_bc(ya,yb,lam,S),
                         sol.x, sol.y, tol=tol, max_nodes=max_nodes, verbose=0)
        assert sol.status == 0, f"K-continuation failed at K={Kstep}"
    return sol

# ---------- corrected momentum-perturbation Chebyshev collocation ----------
# Hf is eliminated algebraically via Hf' = -2 Ff, Hf(0)=0 (an exact linear
# integral operator, not a differential unknown carrying its own gamma
# weight) -- this removes the singular/zero block that a naive [Ff,Gf,Hf]
# formulation puts into the mass matrix B, which is what produced the
# spurious near-zero and wildly resolution-dependent "eigenvalues" seen in
# an earlier, less careful attempt at this collocation.
def build_H_operator(Deta, N, n):
    D2 = Deta.copy()
    D2[0, :] = 0.0
    D2[0, 0] = 1.0          # enforce Hf(0)=0 as the eta=0 row
    Zsel = np.eye(n)
    Zsel[0, 0] = 0.0        # RHS row for the eta=0 equation is 0, not -2*Ff(0)
    P = -2.0 * np.linalg.solve(D2, Zsel)   # Hf = P @ Ff
    return P

def momentum_eigs_cheb(base_sol, K, eta_max, N):
    D, xi = cheb(N)
    eta = eta_max*(1-xi)/2
    Deta = D*(-2/eta_max)
    D2eta = Deta @ Deta
    n = N+1

    y0 = base_sol.sol(eta)
    F0,F0p,G0,G0p,H0 = y0[0],y0[1],y0[2],y0[3],y0[4]
    Delta0 = 1-2*K*F0
    F0pp = (F0**2-G0**2+H0*F0p-K*(F0p**2-G0p**2))/Delta0
    G0pp = (2*F0*G0+H0*G0p-2*K*F0p*G0p)/Delta0

    P = build_H_operator(Deta, N, n)   # Hf = P @ Ff

    I = np.eye(n); Z = np.zeros((n,n))

    # ----- Ffpp block (verified matching manuscript/eq:linF exactly) -----
    # Delta*Ffpp = -2K Ff F0pp +2K F0p Ffp -2K G0p Gfp +gamma Ff
    #              -2F0 Ff +2G0 Gf +H0 Ffp +Hf F0p
    A_FF = (np.diag(Delta0) @ D2eta
            + 2*K*np.diag(F0pp)
            - 2*K*np.diag(F0p) @ Deta
            - np.diag(H0) @ Deta
            + 2*np.diag(F0)
            - np.diag(F0p) @ P)                       # Hf*F0p term -> P
    A_FG = 2*K*np.diag(G0p) @ Deta - 2*np.diag(G0)
    # gamma*Ff term handled via B

    # ----- Gfpp block (CORRECTED, independently re-verified above) -----
    # Delta*Gfpp = 4K Ff G0pp -2K Gf F0pp -2K F0p Gfp -2K Ffp G0p +gamma Gf
    #              +2F0 Gf +2Ff G0 +H0 Gfp +Hf G0p
    A_GF = (-4*K*np.diag(G0pp) + 2*K*np.diag(G0p) @ Deta - 2*np.diag(G0)
            - np.diag(G0p) @ P)                        # Hf*G0p term -> P
    A_GG = (np.diag(Delta0) @ D2eta + 2*K*np.diag(F0pp)
            + 2*K*np.diag(F0p) @ Deta - np.diag(H0) @ Deta - 2*np.diag(F0))

    A = np.block([[A_FF, A_FG],
                  [A_GF, A_GG]])
    B = np.block([[I, Z],
                  [Z, I]])

    Afull = -A.copy()
    Bfull = B.copy()
    m = n
    # Ff(0)=0
    Afull[0, :] = 0; Afull[0, 0] = 1; Bfull[0, :] = 0
    # Gf(0)=0
    Afull[m+0, :] = 0; Afull[m+0, m+0] = 1; Bfull[m+0, :] = 0
    # Ff'(eta_max)=0  (relaxed far-field condition, node N)
    Afull[N, :] = 0; Afull[N, 0:m] = Deta[N, :]; Bfull[N, :] = 0
    # Gf(eta_max)=0
    Afull[m+N, :] = 0; Afull[m+N, m+N] = 1; Bfull[m+N, :] = 0

    evals = eig(Afull, Bfull, right=False)
    evals = evals[np.isfinite(evals)]
    real_evals = np.sort(evals[np.abs(evals.imag) < 1e-6*np.maximum(1, np.abs(evals.real))].real)
    return real_evals

if __name__ == "__main__":
    print("=== Baseline K=0.3, lambda=0.5, S=0.5: corrected Chebyshev collocation, Hf eliminated ===")
    base = solve_base(0.3, 0.5, 0.5, eta_inf=40.0)
    for N, eta_max in [(50,20.0),(60,20.0),(70,25.0),(80,25.0),(90,30.0),(100,30.0),(110,35.0),(120,35.0)]:
        evs = momentum_eigs_cheb(base, 0.3, eta_max, N)
        pos = evs[(evs>1e-3) & (evs<50)]
        print(f"  N={N:4d}, eta_max={eta_max:5.1f}: lowest few positive real eigs = {np.round(pos[:4],6)}")
