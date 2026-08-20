import numpy as np
from scipy.integrate import solve_bvp

# ---------- parameters (Example 1 baseline, matching the paper) ----------
K, lam, S, beta = 0.3, 0.5, 0.5, 0.2
phi1, phi2, phi3 = 0.01, 0.01, 0.01
rho_f, cp_f, k_f = 997.1, 4179.0, 0.613
rho_s1, cp_s1, k_s1 = 8933.0, 385.0, 400.0      # Cu
rho_s2, cp_s2, k_s2 = 3970.0, 765.0, 40.0       # Al2O3
rho_s3, cp_s3, k_s3 = 4250.0, 686.2, 8.9538     # TiO2
mu_f, D_f = 8.9e-4, 1.0e-9
eta_inf, Nmesh0 = 40.0, 800

# ---------- effective ternary hybrid nanofluid properties ----------
rho_thnf = (1-phi3)*((1-phi2)*((1-phi1)*rho_f + phi1*rho_s1) + phi2*rho_s2) + phi3*rho_s3
rhocp_f, rhocp_s1, rhocp_s2, rhocp_s3 = rho_f*cp_f, rho_s1*cp_s1, rho_s2*cp_s2, rho_s3*cp_s3
rhocp_thnf = (1-phi3)*((1-phi2)*((1-phi1)*rhocp_f + phi1*rhocp_s1) + phi2*rhocp_s2) + phi3*rhocp_s3
mu_thnf = mu_f*(1-phi1)**(-2.5)*(1-phi2)**(-2.5)*(1-phi3)**(-2.5)
k_nf  = k_f*(k_s1+2*k_f-2*phi1*(k_f-k_s1))/(k_s1+2*k_f+phi1*(k_f-k_s1))
k_hnf = k_nf*(k_s2+2*k_nf-2*phi2*(k_nf-k_s2))/(k_s2+2*k_nf+phi2*(k_nf-k_s2))
k_thnf= k_hnf*(k_s3+2*k_hnf-2*phi3*(k_hnf-k_s3))/(k_s3+2*k_hnf+phi3*(k_hnf-k_s3))
D_nf  = D_f*2*(1-phi1)/(2+phi1)
D_hnf = D_nf*2*(1-phi2)/(2+phi2)
D_thnf= D_hnf*2*(1-phi3)/(2+phi3)
nu_thnf = mu_thnf/rho_thnf
Pr_thnf = (mu_thnf*rhocp_thnf/rho_thnf)/k_thnf
Sc_thnf = nu_thnf/D_thnf
print(f"Pr_thnf = {Pr_thnf:.4f}   Sc_thnf = {Sc_thnf:.4f}")

# ---------- momentum sub-system (F, F', G, G', H) ----------
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

def solve_momentum(Kv, lamv, Sv, eta_mesh, y_guess, tol=1e-10):
    return solve_bvp(lambda e,y: mom_rhs(e,y,Kv), lambda ya,yb: mom_bc(ya,yb,lamv,Sv),
                      eta_mesh, y_guess, tol=tol, max_nodes=300000, verbose=0)

# ROBUST two-stage solve: (1) grow domain at K=0, (2) grow K to target
print("Stage 1: solving K=0 at small domain, then extending...")
mesh = np.linspace(0, 8, 200)
mom = solve_momentum(0.0, lam, S, mesh, mom_guess(mesh, lam, S))
for eta_target in [10, 15, 20, 25, 30, eta_inf]:
    if eta_target <= mom.x[-1]:
        continue
    extra = np.linspace(mom.x[-1], eta_target, 100)[1:]
    new_x = np.concatenate([mom.x, extra])
    new_y = np.concatenate([mom.y, np.tile(mom.y[:, -1:], (1, extra.size))], axis=1)
    mom = solve_momentum(0.0, lam, S, new_x, new_y)
print(f"  K=0 stable solution at eta_inf={eta_inf}: F'(0) = {mom.y[1,0]:.6f}")

print("Stage 2: stepping K from 0 up to target...")
for Kstep in np.linspace(0, K, 31)[1:]:
    mom = solve_momentum(Kstep, lam, S, mom.x, mom.y)
print(f"  Final momentum result: F'(0) = {mom.y[1,0]:.6f}   -G'(0) = {-mom.y[3,0]:.6f}")

# ---------- energy sub-system (theta, theta') ----------
H_interp = lambda e: mom.sol(e)[4]
def energy_rhs(eta, y):
    th, thp = y
    return np.vstack([thp, Pr_thnf*H_interp(eta)*thp])
en_mesh = np.linspace(mom.x[0], mom.x[-1], 300)
en = solve_bvp(energy_rhs, lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
               en_mesh, np.vstack([np.exp(-en_mesh), -np.exp(-en_mesh)]), tol=1e-11, max_nodes=200000)
print(f"  -theta'(0) = {-en.y[1,0]:.6f}")

# ---------- species sub-system (phi, phi') ----------
def species_rhs(eta, y):
    ph, php = y
    return np.vstack([php, Sc_thnf*H_interp(eta)*php + Sc_thnf*beta*ph])
sp_mesh = np.unique(np.concatenate([np.linspace(mom.x[0], mom.x[-1], 200),
                                     np.geomspace(1e-4, mom.x[-1], 400)]))
sp = solve_bvp(species_rhs, lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
               sp_mesh, np.vstack([np.exp(-sp_mesh), -np.exp(-sp_mesh)]), tol=1e-11, max_nodes=400000)
print(f"  -phi'(0)   = {-sp.y[1,0]:.6f}")

print()
print("=== COMPARE AGAINST Main.tex Example 1 (should match closely) ===")
print("Paper:  F'(0)=-0.442969   -G'(0)=1.368066   -theta'(0)=3.631009   -phi'(0)=443.932186")
print(f"Yours:  F'(0)={mom.y[1,0]:.6f}   -G'(0)={-mom.y[3,0]:.6f}   -theta'(0)={-en.y[1,0]:.6f}   -phi'(0)={-sp.y[1,0]:.6f}")
