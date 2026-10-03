import numpy as np
from scipy.integrate import solve_bvp
import time

# Example 4: thermal and species wall gradients across shrinking
# parameter, nanoparticle loading, and reaction parameter (K=0.3, S=1.5).
# Reproduces the table in Main.tex. Takes about 20-30 seconds in Colab.

K, S = 0.3, 1.5
rho_f, cp_f, k_f = 997.1, 4179.0, 0.613
rho_s1, cp_s1, k_s1 = 8933.0, 385.0, 400.0
rho_s2, cp_s2, k_s2 = 3970.0, 765.0, 40.0
rho_s3, cp_s3, k_s3 = 4250.0, 686.2, 8.9538
mu_f, D_f = 8.9e-4, 1.0e-9

def props(phi):
    p1=p2=p3=phi
    rho_thnf=(1-p3)*((1-p2)*((1-p1)*rho_f+p1*rho_s1)+p2*rho_s2)+p3*rho_s3
    rhocp_thnf=(1-p3)*((1-p2)*((1-p1)*rho_f*cp_f+p1*rho_s1*cp_s1)+p2*rho_s2*cp_s2)+p3*rho_s3*cp_s3
    mu_thnf=mu_f*(1-p1)**-2.5*(1-p2)**-2.5*(1-p3)**-2.5
    k_nf=k_f*(k_s1+2*k_f-2*p1*(k_f-k_s1))/(k_s1+2*k_f+p1*(k_f-k_s1))
    k_hnf=k_nf*(k_s2+2*k_nf-2*p2*(k_nf-k_s2))/(k_s2+2*k_nf+p2*(k_nf-k_s2))
    k_thnf=k_hnf*(k_s3+2*k_hnf-2*p3*(k_hnf-k_s3))/(k_s3+2*k_hnf+p3*(k_hnf-k_s3))
    D_nf=D_f*2*(1-p1)/(2+p1); D_hnf=D_nf*2*(1-p2)/(2+p2); D_thnf=D_hnf*2*(1-p3)/(2+p3)
    nu_thnf=mu_thnf/rho_thnf
    Pr=(mu_thnf*rhocp_thnf/rho_thnf)/k_thnf
    Sc=nu_thnf/D_thnf
    return Pr, Sc

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

def solve_momentum_at(lamv):
    mesh = np.linspace(0, 8, 200)
    mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lamv,S),
                     mesh, mom_guess(mesh, lamv, S), tol=1e-10, max_nodes=300000, verbose=0)
    for eta_t in [15, 25, 40]:
        extra = np.linspace(mom.x[-1], eta_t, 100)[1:]
        new_x = np.concatenate([mom.x, extra])
        new_y = np.concatenate([mom.y, np.tile(mom.y[:,-1:], (1, extra.size))], axis=1)
        mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lamv,S),
                         new_x, new_y, tol=1e-10, max_nodes=300000, verbose=0)
    for Kstep in np.linspace(0, K, 31)[1:]:
        mom = solve_bvp(lambda e,y: mom_rhs(e,y,Kstep), lambda ya,yb: mom_bc(ya,yb,lamv,S),
                         mom.x, mom.y, tol=1e-10, max_nodes=300000, verbose=0)
    return mom

t0 = time.time()
rows = []
for lamv in [0.5, -1.0, -2.0]:
    mom = solve_momentum_at(lamv)
    H_interp = lambda e: mom.sol(e)[4]
    for phi in [0.0, 0.01, 0.03]:
        Pr, Sc = props(phi)
        en_mesh = np.linspace(mom.x[0], mom.x[-1], 300)
        en = solve_bvp(lambda e,y: np.vstack([y[1], Pr*H_interp(e)*y[1]]),
                        lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
                        en_mesh, np.vstack([np.exp(-en_mesh), -np.exp(-en_mesh)]), tol=1e-11, max_nodes=200000, verbose=0)
        for beta in [0.0, 0.5]:
            sp_mesh = np.unique(np.concatenate([np.linspace(mom.x[0], mom.x[-1], 200),
                                                 np.geomspace(1e-4, mom.x[-1], 400)]))
            sp = solve_bvp(lambda e,y: np.vstack([y[1], Sc*H_interp(e)*y[1] + Sc*beta*y[0]]),
                            lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
                            sp_mesh, np.vstack([np.exp(-sp_mesh), -np.exp(-sp_mesh)]), tol=1e-11, max_nodes=400000, verbose=0)
            rows.append((lamv, phi, beta, Pr, -en.y[1,0], -sp.y[1,0]))

print(f"Done in {time.time()-t0:.1f} seconds.")
print()
print("=== COMPARE AGAINST Main.tex Table (Example 4) ===")
print(f"{'lambda':>7} {'phi':>6} {'beta':>5} {'Pr_thnf':>9} {'-theta_p0':>11} {'-phi_p0':>11}")
for r in rows:
    print(f"{r[0]:7.1f} {r[1]:6.2f} {r[2]:5.1f} {r[3]:9.4f} {r[4]:11.4f} {r[5]:11.3f}")
print()
print("Compare each row above against Table 'Example 4' in Main.tex --")
print("they should match to 3-4 decimal places.")
