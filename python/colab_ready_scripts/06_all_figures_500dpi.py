import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time, os

# Generates the six baseline figures within this script's scope, at 500 DPI
# (plus vector PDF), and auto-downloads them (works in Google Colab; if
# not running in Colab, files are just saved in the current folder
# instead). Takes about 1-2 minutes total. The remaining two manuscript
# figures (the corrected eigenvalue spectrum and its trend along the
# branch, and the H(eta) mechanism comparison) are generated separately
# by scripts 14, 15, and a dedicated regeneration step, since they depend
# on the corrected stability equations fixed after this script was
# originally written.

DPI = 500
K, lam, S, beta = 0.3, 0.5, 0.5, 0.2
phi1 = phi2 = phi3 = 0.01
rho_f, cp_f, k_f = 997.1, 4179.0, 0.613
rho_s1, cp_s1, k_s1 = 8933.0, 385.0, 400.0
rho_s2, cp_s2, k_s2 = 3970.0, 765.0, 40.0
rho_s3, cp_s3, k_s3 = 4250.0, 686.2, 8.9538
mu_f, D_f = 8.9e-4, 1.0e-9

rho_thnf=(1-phi3)*((1-phi2)*((1-phi1)*rho_f+phi1*rho_s1)+phi2*rho_s2)+phi3*rho_s3
rhocp_thnf=(1-phi3)*((1-phi2)*((1-phi1)*rho_f*cp_f+phi1*rho_s1*cp_s1)+phi2*rho_s2*cp_s2)+phi3*rho_s3*cp_s3
mu_thnf=mu_f*(1-phi1)**-2.5*(1-phi2)**-2.5*(1-phi3)**-2.5
k_nf=k_f*(k_s1+2*k_f-2*phi1*(k_f-k_s1))/(k_s1+2*k_f+phi1*(k_f-k_s1))
k_hnf=k_nf*(k_s2+2*k_nf-2*phi2*(k_nf-k_s2))/(k_s2+2*k_nf+phi2*(k_nf-k_s2))
k_thnf=k_hnf*(k_s3+2*k_hnf-2*phi3*(k_hnf-k_s3))/(k_s3+2*k_hnf+phi3*(k_hnf-k_s3))
D_nf=D_f*2*(1-phi1)/(2+phi1); D_hnf=D_nf*2*(1-phi2)/(2+phi2); D_thnf=D_hnf*2*(1-phi3)/(2+phi3)
nu_thnf=mu_thnf/rho_thnf
Pr_thnf=(mu_thnf*rhocp_thnf/rho_thnf)/k_thnf
Sc_thnf=nu_thnf/D_thnf

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

print("[1/4] Solving baseline (robust two-stage)...")
t0=time.time()
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
en_mesh = np.linspace(mom.x[0], mom.x[-1], 300)
en = solve_bvp(lambda e,y: np.vstack([y[1], Pr_thnf*H_interp(e)*y[1]]),
               lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
               en_mesh, np.vstack([np.exp(-en_mesh), -np.exp(-en_mesh)]), tol=1e-11, max_nodes=200000, verbose=0)
sp_mesh = np.unique(np.concatenate([np.linspace(mom.x[0], mom.x[-1], 200), np.geomspace(1e-4, mom.x[-1], 400)]))
sp = solve_bvp(lambda e,y: np.vstack([y[1], Sc_thnf*H_interp(e)*y[1] + Sc_thnf*beta*y[0]]),
               lambda ya,yb: np.array([ya[0]-1.0, yb[0]]),
               sp_mesh, np.vstack([np.exp(-sp_mesh), -np.exp(-sp_mesh)]), tol=1e-11, max_nodes=400000, verbose=0)
print(f"      done ({time.time()-t0:.1f}s)")

print("[2/4] Running continuation search (Example 2)...")
t0=time.time()
NMAX = 6000
lam2, S2 = 0.0, 1.5
mesh2 = np.linspace(0, 15, 400)
sol2 = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam2,S2),
                  mesh2, mom_guess(mesh2, lam2, S2), tol=1e-9, max_nodes=NMAX, verbose=0)
dlam = -0.05
lams, Fp0s, maxFs = [lam2], [sol2.y[1,0]], [sol2.y[0].max()]
while lam2 > -6.0:
    trial = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam2+dlam,S2),
                       sol2.x, sol2.y, tol=1e-9, max_nodes=NMAX, verbose=0)
    if trial.status != 0:
        dlam /= 2
        if abs(dlam) < 1e-5:
            break
        continue
    sol2 = trial
    lam2 += dlam
    lams.append(lam2); Fp0s.append(sol2.y[1,0]); maxFs.append(sol2.y[0].max())
print(f"      done ({time.time()-t0:.1f}s), {len(lams)} points, stalled near lambda={lam2:.4f}")

print("[3/4] Searching for the momentum eigenvalue (Example 3)...")
t0=time.time()
def pert_rhs(eta, z, gamma, base_interp):
    y0 = base_interp(eta)
    F0,F0p,G0,G0p,H0 = y0
    Delta = 1 - 2*K*F0
    if abs(Delta) < 1e-10: Delta = 1e-10 if Delta>=0 else -1e-10
    F0pp = (F0**2-G0**2+H0*F0p-K*(F0p**2-G0p**2))/Delta
    G0pp = (2*F0*G0+H0*G0p-2*K*F0p*G0p)/Delta
    Ff,Ffp,Gf,Gfp,Hf = z
    Ffpp = (2*K*Ff*F0pp - 2*K*F0p*Ffp + 2*K*G0p*Gfp - gamma*Ff + 2*F0*Ff - 2*G0*Gf + H0*Ffp + Hf*F0p)/Delta
    Gfpp = (2*K*Ff*G0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p - gamma*Gf + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p)/Delta
    return [Ffp, Ffpp, Gfp, Gfpp, -2*Ff]
def far_field(gamma, alpha, base_interp, eta_max):
    s = solve_ivp(pert_rhs, [0, eta_max], [0,1,0,alpha,0], args=(gamma, base_interp),
                  method='RK45', rtol=1e-8, atol=1e-10, max_step=eta_max/50)
    if not s.success or np.any(np.abs(s.y) > 1e8): return 1e8, 1e8
    return s.y[0,-1], s.y[2,-1]
def alpha_root(gamma, base_interp, eta_max):
    alphas = np.linspace(-3, 3, 25)
    vals = [far_field(gamma, a, base_interp, eta_max)[0] for a in alphas]
    for i in range(len(alphas)-1):
        if np.sign(vals[i]) != np.sign(vals[i+1]):
            try: return brentq(lambda a: far_field(gamma, a, base_interp, eta_max)[0], alphas[i], alphas[i+1], xtol=1e-8)
            except Exception: pass
    return None
def gamma_residual(gamma, base_interp, eta_max):
    a = alpha_root(gamma, base_interp, eta_max)
    if a is None: return None
    return far_field(gamma, a, base_interp, eta_max)[1]

gammas_scan = np.linspace(-2, 2, 13)
resids = [gamma_residual(g, mom.sol, 8.0) for g in gammas_scan]
gamma_M = None
for i in range(len(gammas_scan)-1):
    if resids[i] is not None and resids[i+1] is not None and np.sign(resids[i]) != np.sign(resids[i+1]):
        gamma_M = brentq(lambda g: gamma_residual(g, mom.sol, 8.0) or 1e8, gammas_scan[i], gammas_scan[i+1], xtol=1e-6)
        break
print(f"      done ({time.time()-t0:.1f}s), gamma_1^(M) = {gamma_M:.4f}")

print("[4/4] Drawing and saving figures at 500 DPI...")
outdir = "paper_figures"
os.makedirs(outdir, exist_ok=True)
saved = []

def save(fig, name):
    path_png = os.path.join(outdir, name)
    fig.savefig(path_png, dpi=DPI, bbox_inches='tight')
    path_pdf = os.path.join(outdir, name.replace('.png', '.pdf'))
    fig.savefig(path_pdf, bbox_inches='tight')  # vector PDF for journal submission
    plt.close(fig)
    saved.append(path_png)
    saved.append(path_pdf)

eta_plot = np.linspace(0, 8, 300)
F,Fp,G,Gp,H = mom.sol(eta_plot)
fig, ax = plt.subplots(figsize=(6,4.2))
ax.plot(eta_plot, F, label=r'$F(\eta)$')
ax.plot(eta_plot, G, label=r'$G(\eta)$')
ax.plot(eta_plot, -H, label=r'$-H(\eta)$')
ax.set_xlabel(r'$\eta$'); ax.set_ylabel('Velocity components')
ax.set_title('Baseline velocity profiles (K=0.3, $\\lambda$=0.5, S=0.5)')
ax.legend(); ax.grid(alpha=0.3)
save(fig, "velocity_profiles.png")

fig, ax = plt.subplots(figsize=(6,4.2))
th_plot = en.sol(eta_plot)[0]
ax.plot(eta_plot, th_plot)
ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$\theta(\eta)$')
ax.set_title('Baseline temperature profile')
ax.grid(alpha=0.3)
save(fig, "temperature_profile.png")

eta_c = np.linspace(0, 0.3, 300)
fig, ax = plt.subplots(figsize=(6,4.2))
ax.plot(eta_c, sp.sol(eta_c)[0])
ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$\phi(\eta)$')
ax.set_title(f'Baseline concentration profile (Sc$\\approx${Sc_thnf:.0f})')
ax.grid(alpha=0.3)
save(fig, "concentration_profile.png")

fig, ax = plt.subplots(figsize=(6.5,4.5))
ax.plot(lams, Fp0s, '-', lw=1.5)
ax.axvline(lams[-1], color='red', ls='--', lw=1, label=r'$\lambda$ near singularity')
ax.set_xlabel(r'$\lambda$ (shrinking parameter)'); ax.set_ylabel(r"$F'(0)$")
ax.set_title('Momentum branch, K=0.3, S=1.5')
ax.legend(); ax.grid(alpha=0.3)
save(fig, "branch_diagram.png")

fig, ax = plt.subplots(figsize=(6.5,4.5))
ax.plot(lams, maxFs, '-', lw=1.5, color='tab:orange')
ax.axhline(1/(2*K), color='red', ls='--', lw=1, label=r'analytical singularity $F=1/(2K)$')
ax.set_xlabel(r'$\lambda$'); ax.set_ylabel(r'$\max_\eta F(\eta)$')
ax.set_title('Approach to the regularity boundary')
ax.legend(); ax.grid(alpha=0.3)
save(fig, "regularity_boundary.png")

fig, ax = plt.subplots(figsize=(6.5,4.5))
ax.plot(gammas_scan, resids, 'o-')
ax.axhline(0, color='k', lw=0.7)
if gamma_M is not None:
    ax.axvline(gamma_M, color='red', ls='--', lw=1, label=fr'root $\gamma_1\approx${gamma_M:.3f}')
ax.set_xlabel(r'$\gamma$ (trial eigenvalue)'); ax.set_ylabel('shooting residual')
ax.set_yscale('symlog')
ax.set_title('Momentum eigenvalue shooting residual, baseline case')
ax.legend(); ax.grid(alpha=0.3)
save(fig, "eigenvalue_spectrum.png")

print(f"      saved {len(saved)} figures to ./{outdir}/ at {DPI} DPI")

# ---- auto-download (Colab only; harmless no-op elsewhere) ----
try:
    from google.colab import files
    for p in saved:
        files.download(p)
    print()
    print("Download prompts triggered for all 6 figures -- check your browser's downloads.")
except ImportError:
    print()
    print(f"Not running in Colab -- figures are saved locally in ./{outdir}/")
