import numpy as np
from scipy.integrate import solve_bvp

# Reproduces four new results added to the manuscript during peer-review
# incorporation: (1) monotonic lambda_c vs K trend across 5 K values,
# (2) extended beta sweep confirming near-insensitivity, (3) numerical
# illustration of Corollary 1 (wall-shear factor), (4) the H(eta)
# zero-crossing that explains theta'(0)->0 at lambda=-2.
# Takes about 2-3 minutes total in Colab.

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

def find_lambda_c(K, S, NMAX=6000):
    lam = 0.0
    mesh = np.linspace(0, 15, 400)
    sol = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     mesh, mom_guess(mesh, lam, S), tol=1e-9, max_nodes=NMAX, verbose=0)
    dlam = -0.05
    last_lam = lam
    bound = -30.0/K if K < 0.1 else -20.0
    while lam > bound:
        trial = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam+dlam,S),
                           sol.x, sol.y, tol=1e-9, max_nodes=NMAX, verbose=0)
        if trial.status != 0:
            dlam /= 2
            if abs(dlam) < 1e-5:
                break
            continue
        sol = trial
        lam += dlam
        last_lam = lam
    lo, hi = last_lam, last_lam + dlam
    mesh_lo, y_lo = sol.x, sol.y
    for _ in range(15):
        mid = 0.5*(lo+hi)
        trial = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,mid,S),
                           mesh_lo, y_lo, tol=1e-9, max_nodes=NMAX, verbose=0)
        if trial.status == 0:
            lo = mid; mesh_lo, y_lo = trial.x, trial.y
        else:
            hi = mid
        if abs(hi-lo) < 1e-4:
            break
    return lo

print("=== (1) Monotonic lambda_c vs K trend ===")
print("Paper: lambda_c = -17.01, -10.22, -7.83, -6.49, -5.03 at K = 0.05, 0.10, 0.15, 0.20, 0.30")
print("Yours:")
for K in [0.05, 0.1, 0.15, 0.2, 0.3]:
    lc = find_lambda_c(K, 1.5)
    print(f"  K={K:.2f}  lambda_c={lc:.2f}  (1/(2K)={1/(2*K):.2f})")

print()
print("=== (2) Extended beta sweep (Example 4, K=0.3, lambda=0.5, S=1.5, phi=0.01) ===")
K, S, lam = 0.3, 1.5, 0.5
phi1=phi2=phi3=0.01
rho_f, cp_f, k_f = 997.1, 4179.0, 0.613
rho_s1, cp_s1, k_s1 = 8933.0, 385.0, 400.0
rho_s2, cp_s2, k_s2 = 3970.0, 765.0, 40.0
rho_s3, cp_s3, k_s3 = 4250.0, 686.2, 8.9538
mu_f, D_f = 8.9e-4, 1.0e-9
rho_thnf=(1-phi3)*((1-phi2)*((1-phi1)*rho_f+phi1*rho_s1)+phi2*rho_s2)+phi3*rho_s3
mu_thnf=mu_f*(1-phi1)**-2.5*(1-phi2)**-2.5*(1-phi3)**-2.5
k_nf=k_f*(k_s1+2*k_f-2*phi1*(k_f-k_s1))/(k_s1+2*k_f+phi1*(k_f-k_s1))
k_hnf=k_nf*(k_s2+2*k_nf-2*phi2*(k_nf-k_s2))/(k_s2+2*k_nf+phi2*(k_nf-k_s2))
k_thnf=k_hnf*(k_s3+2*k_hnf-2*phi3*(k_hnf-k_s3))/(k_s3+2*k_hnf+phi3*(k_hnf-k_s3))
D_nf=D_f*2*(1-phi1)/(2+phi1); D_hnf=D_nf*2*(1-phi2)/(2+phi2); D_thnf=D_hnf*2*(1-phi3)/(2+phi3)
nu_thnf=mu_thnf/rho_thnf
Sc_thnf=nu_thnf/D_thnf

mesh = np.linspace(0, 8, 200)
mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lam,S),
                 mesh, mom_guess(mesh, lam, S), tol=1e-10, max_nodes=300000, verbose=0)
for eta_t in [15, 25, 40]:
    extra = np.linspace(mom.x[-1], eta_t, 100)[1:]
    new_x = np.concatenate([mom.x, extra])
    new_y = np.concatenate([mom.y, np.tile(mom.y[:,-1:], (1, extra.size))], axis=1)
    mom = solve_bvp(lambda e,y: mom_rhs(e,y,0.0), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     new_x, new_y, tol=1e-10, max_nodes=300000, verbose=0)
for Kstep in np.linspace(0, K, 31)[1:]:
    mom = solve_bvp(lambda e,y: mom_rhs(e,y,Kstep), lambda ya,yb: mom_bc(ya,yb,lam,S),
                     mom.x, mom.y, tol=1e-10, max_nodes=300000, verbose=0)
H_interp = lambda e: mom.sol(e)[4]
print("Paper: -phi'(0) = 1325.34, 1325.68, 1326.01, 1326.67, 1327.34 at beta = 0, 0.5, 1, 2, 3")
print("Yours:")
for beta in [0.0, 0.5, 1.0, 2.0, 3.0]:
    sp_mesh = np.unique(np.concatenate([np.linspace(mom.x[0], mom.x[-1], 200), np.geomspace(1e-4, mom.x[-1], 400)]))
    sp = solve_bvp(lambda e,y: np.vstack([y[1], Sc_thnf*H_interp(e)*y[1] + Sc_thnf*beta*y[0]]),
                    lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
                    sp_mesh, np.vstack([np.exp(-sp_mesh), -np.exp(-sp_mesh)]), tol=1e-11, max_nodes=400000, verbose=0)
    print(f"  beta={beta:.1f}   -phi'(0) = {-sp.y[1,0]:.2f}")

print()
print("=== (3) Corollary 1 numerical illustration (Example 1 baseline) ===")
Fp0, negGp0 = -0.442969, 1.368066
factor = 1 - 2*K*lam
print(f"Paper: (1-2*K*lambda)*F'(0) = -0.3101, -(1-2*K*lambda)*G'(0) = 0.9576")
print(f"Yours: (1-2*K*lambda)*F'(0) = {factor*Fp0:.4f}, -(1-2*K*lambda)*G'(0) = {factor*negGp0:.4f}")

print()
print("=== (4) H(eta) zero-crossing explaining theta'(0)->0 at lambda=-2 ===")
lam2 = -2.0
mesh2 = np.linspace(0, 15, 400)
sol2 = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam2,S),
                  mesh2, mom_guess(mesh2, lam2, S), tol=1e-10, max_nodes=200000, verbose=0)
from scipy.optimize import brentq
zero_cross = brentq(lambda e: sol2.sol(e)[4], 0.01, 2.0)
print(f"Paper: H(eta) crosses zero at eta ~= 0.50")
print(f"Yours: H(eta) crosses zero at eta = {zero_cross:.4f}")
