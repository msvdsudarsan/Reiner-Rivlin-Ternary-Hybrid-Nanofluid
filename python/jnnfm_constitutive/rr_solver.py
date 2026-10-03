import numpy as np
from scipy.integrate import solve_bvp

def mom_rhs(eta, y, K):
    F, Fp, G, Gp, H = y
    D = 1 - 2*K*F
    D = np.where(np.abs(D) < 1e-10, np.sign(D)*1e-10 + 1e-14, D)
    Fpp = (F**2 - G**2 + H*Fp - K*(Fp**2 - Gp**2)) / D
    Gpp = (2*F*G + H*Gp - 2*K*Fp*Gp) / D
    return np.vstack([Fp, Fpp, Gp, Gpp, -2*F])

def mom_bc(ya, yb, lam, S):
    return np.array([ya[0]-lam, ya[2]-1.0, ya[4]+S, yb[0], yb[2]])

def guess(eta, lam, S):
    d = np.exp(-eta)
    return np.vstack([lam*d, -lam*d, d, -d, np.full_like(eta, -S)])

def solve(sol_prev, K, lam, S, tol=1e-9, mx=20000, einf=15.0):
    if sol_prev is None:
        eta = np.linspace(0, einf, 500); y0 = guess(eta, lam, S)
    else:
        eta, y0 = sol_prev.x, sol_prev.y
    return solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda a,b: mom_bc(a,b,lam,S),
                     eta, y0, tol=tol, max_nodes=mx, verbose=0)

def continue_lambda(K, S, dlam=-0.02, lam_min=-30, tol=1e-9):
    """Natural continuation in lambda from 0 downward. Returns (lam_c, minF, maxF, sol)."""
    sol = None; lam = 0.0; step = dlam; stalls = 0; last = None
    while lam > lam_min and stalls < 6:
        nxt = lam + step
        s = solve(sol, K, nxt, S, tol=tol)
        if s.status == 0:
            sol, lam = s, nxt
            e = np.linspace(0, 15, 400); Fv = s.sol(e)[0]
            last = (lam, Fv.min(), Fv.max())
            stalls = 0
        else:
            step /= 2; stalls += 1
            if abs(step) < 1e-5: break
    return last, sol
