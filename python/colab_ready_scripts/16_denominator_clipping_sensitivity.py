import numpy as np
from scipy.integrate import solve_bvp
import sys

K, S = 0.3, 1.5
tol_clip = float(sys.argv[1])

def mom_rhs_clipped(eta, y, Kv, tol):
    F, Fp, G, Gp, H = y
    Delta = 1 - 2*Kv*F
    Delta = np.where(np.abs(Delta) < tol, np.sign(Delta)*tol + tol*1e-4, Delta)
    Fpp = (F**2 - G**2 + H*Fp - Kv*(Fp**2 - Gp**2)) / Delta
    Gpp = (2*F*G + H*Gp - 2*Kv*Fp*Gp) / Delta
    return np.vstack([Fp, Fpp, Gp, Gpp, -2*F])

def mom_bc(ya, yb, lamv, Sv):
    return np.array([ya[0]-lamv, ya[2]-1.0, ya[4]+Sv, yb[0], yb[2]])
def mom_guess(eta, lamv, Sv):
    d = np.exp(-eta)
    return np.vstack([lamv*d, -lamv*d, d, -d, np.full_like(eta, -Sv)])

def solve_one(sol_prev, K, lam, S, tol_clip, max_nodes=15000):
    if sol_prev is None:
        eta = np.linspace(0, 35.0, 600); y0 = mom_guess(eta, lam, S)
    else:
        eta, y0 = sol_prev.x, sol_prev.y
    return solve_bvp(lambda e,y: mom_rhs_clipped(e,y,K,tol_clip), lambda ya,yb: mom_bc(ya,yb,lam,S),
                      eta, y0, tol=1e-8, max_nodes=max_nodes, verbose=0)

sol = None
for Kstep in np.linspace(0, K, 16)[1:]:
    sol = solve_one(sol, Kstep, 0.0, S, tol_clip)
    if sol.status != 0:
        print("K-continuation failed"); sys.exit()

lam = 0.0; dlam = -0.05; last_good_lam = lam; n = 0
while lam > -5.02 and n < 200:
    n += 1
    lam_next = max(lam+dlam, -5.02)
    s = solve_one(sol, K, lam_next, S, tol_clip)
    if s.status == 0:
        sol, lam = s, lam_next; last_good_lam = lam
        if abs(dlam) < 0.05: dlam = np.sign(dlam)*min(abs(dlam)*1.5, 0.05)
    else:
        dlam /= 2.0
        if abs(dlam) < 1e-4: break

maxF = np.max(sol.y[0])
print(f"clip tol={tol_clip:.0e}: reached lambda={last_good_lam:.6f} max F={maxF:.6f} (1/2K={1/(2*K):.6f})")
