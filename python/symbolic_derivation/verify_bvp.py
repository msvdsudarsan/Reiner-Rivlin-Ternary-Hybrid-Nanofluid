import numpy as np
from scipy.integrate import solve_bvp

def make_rhs(K, Pr, Sc, beta):
    def rhs(eta, y):
        F,Fp,G,Gp,H,th,thp,ph,php = y
        Delta = 1 - 2*K*F
        Delta = np.where(np.abs(Delta) < 1e-8, np.sign(Delta)*1e-8 + 1e-12, Delta)
        Fpp = (F**2 - G**2 + H*Fp - K*(Fp**2 - Gp**2)) / Delta
        Gpp = (2*F*G + H*Gp - 2*K*Fp*Gp) / Delta
        Hp = -2*F
        thpp = Pr*H*thp
        phpp = Sc*H*php + Sc*beta*ph
        return np.vstack([Fp,Fpp,Gp,Gpp,Hp,thp,thpp,php,phpp])
    return rhs

def make_bc(lam, S):
    def bc(ya, yb):
        return np.array([
            ya[0]-lam, ya[2]-1, ya[4]+S, ya[5]-1, ya[7]-1,
            yb[0], yb[2], yb[5], yb[7]
        ])
    return bc

def solve_case(K, lam, S, Pr=1.0, Sc=1.0, beta=0.0, eta_inf=15, n=400):
    eta = np.linspace(0, eta_inf, n)
    y0 = np.zeros((9, eta.size))
    decay = np.exp(-eta)
    y0[0] = lam*decay
    y0[1] = -lam*decay
    y0[2] = decay
    y0[3] = -decay
    y0[4] = -S + lam*(1-decay)
    y0[5] = decay
    y0[6] = -decay
    y0[7] = decay
    y0[8] = -decay

    sol = solve_bvp(make_rhs(K,Pr,Sc,beta), make_bc(lam,S), eta, y0, tol=1e-9, max_nodes=100000, verbose=0)
    return sol

print("=== Newtonian check: K=0, lambda=0, S=0 (classical von Karman) ===")
sol = solve_case(K=0.0, lam=0.0, S=0.0, eta_inf=20, n=600)
print("status:", sol.status, sol.message)
print("F'(0)  =", sol.y[1,0], "  (reference 0.5102)")
print("-G'(0) =", -sol.y[3,0], "  (reference 0.6159)")
print("-H(inf)=", -sol.y[4,-1], "  (reference 0.8845)")

print()
print("=== Reiner-Rivlin check: K=0.3, lambda=0.5 (stretching), S=0.5 ===")
sol2 = solve_case(K=0.3, lam=0.5, S=0.5, Pr=6.2, Sc=1.5, beta=0.2, eta_inf=15, n=400)
print("status:", sol2.status, sol2.message)
print("F'(0)  =", sol2.y[1,0])
print("-G'(0) =", -sol2.y[3,0])
print("-theta'(0) =", -sol2.y[6,0])
print("-phi'(0)   =", -sol2.y[8,0])
print("max F(eta) =", sol2.y[0].max(), " (must stay below 1/(2K)=", 1/(2*0.3), "for regularity)")

print()
print("=== Shrinking-disk sweep: K=0.3, S=1.5, scan lambda from 0 to -2 ===")
lam_grid = np.linspace(0, -2.0, 41)
Fp0_vals = []
prev_sol = None
for lam in lam_grid:
    if prev_sol is None:
        s = solve_case(K=0.3, lam=lam, S=1.5, eta_inf=15, n=400)
    else:
        eta = prev_sol.x
        s = solve_bvp(make_rhs(0.3,1.0,1.0,0.0), make_bc(lam,1.5), eta, prev_sol.y, tol=1e-8, max_nodes=100000)
    if s.status != 0:
        print(f"lambda={lam:.3f}: FAILED ({s.message})")
        break
    Fp0_vals.append(s.y[1,0])
    prev_sol = s
    print(f"lambda={lam:.4f}  F'(0)={s.y[1,0]:.6f}")

Fp0_vals = np.array(Fp0_vals)
d = np.diff(Fp0_vals)
sign_changes = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
print("Sign changes (possible fold) at lambda indices:", sign_changes, 
      [lam_grid[i+1] for i in sign_changes] if len(sign_changes) else "none in range")

print()
print("=== Extending shrinking sweep further, and trying different S ===")
for Stest in [0.5, 1.5, 3.0]:
    lam_grid = np.linspace(0, -5.0, 51)
    Fp0_vals = []
    prev_sol = None
    failed_at = None
    for lam in lam_grid:
        if prev_sol is None:
            s = solve_case(K=0.3, lam=lam, S=Stest, eta_inf=15, n=400)
        else:
            eta = prev_sol.x
            s = solve_bvp(make_rhs(0.3,1.0,1.0,0.0), make_bc(lam,Stest), eta, prev_sol.y, tol=1e-8, max_nodes=100000)
        if s.status != 0:
            failed_at = lam
            break
        Fp0_vals.append(s.y[1,0])
        prev_sol = s
    Fp0_vals = np.array(Fp0_vals)
    d = np.diff(Fp0_vals) if len(Fp0_vals)>1 else np.array([])
    sign_changes = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0] if len(d)>1 else []
    print(f"S={Stest}: reached lambda={lam_grid[len(Fp0_vals)-1] if len(Fp0_vals) else 'NA'}, "
          f"failed_at={failed_at}, sign_changes={list(sign_changes)}, "
          f"F'(0) range=[{Fp0_vals.min():.4f},{Fp0_vals.max():.4f}]" if len(Fp0_vals) else "no data")
