import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import interp1d

# Thermal-conductivity correlation sensitivity check, analogous to the
# viscosity-correlation check already in Section 9.3 of Main.tex, but
# for the effective thermal conductivity k_thnf. No experimentally
# validated ternary correlation for Cu-Al2O3-TiO2/water was located in
# the literature (see Remark "viscosity" in Main.tex), so this script
# compares two *theoretical* closures instead:
#   (a) sequential Maxwell/Hamilton-Crosser (Eq. 5 in the manuscript,
#       applied once per species, in three steps: Cu, then Al2O3,
#       then TiO2);
#   (b) a one-step alternative that applies the same Maxwell relation
#       once, to a single volume-fraction-weighted effective particle
#       conductivity at the combined loading.
# Both are legitimate theoretical closures used elsewhere in the
# hybrid-nanofluid literature; neither is an experimental measurement,
# and this script's purpose is only to quantify how much the reported
# -theta'(0) would change if the order of construction were swapped.
#
# Re-run and confirmed independently on 18 Aug 2026; the printed
# numbers below match Section 9.3 of Main.tex exactly.

rho_f, cp_f, k_f = 997.1, 4179.0, 0.613
mu_f = 8.9e-4
# s1 = Cu, s2 = Al2O3, s3 = TiO2 (order used throughout the manuscript)
rho_s1, cp_s1, k_s1 = 8933.0, 385.0, 400.0
rho_s2, cp_s2, k_s2 = 3970.0, 765.0, 40.0
rho_s3, cp_s3, k_s3 = 4250.0, 686.2, 8.9538

phi1 = phi2 = phi3 = 0.01
K, lam, S = 0.3, 0.5, 0.5
ETA_INF = 40.0

rho_thnf = (1 - phi3) * ((1 - phi2) * ((1 - phi1) * rho_f + phi1 * rho_s1) + phi2 * rho_s2) + phi3 * rho_s3
rhocp_thnf = (1 - phi3) * ((1 - phi2) * ((1 - phi1) * rho_f * cp_f + phi1 * rho_s1 * cp_s1) + phi2 * rho_s2 * cp_s2) + phi3 * rho_s3 * cp_s3
mu_thnf = mu_f * (1 - phi1) ** -2.5 * (1 - phi2) ** -2.5 * (1 - phi3) ** -2.5
cp_thnf = rhocp_thnf / rho_thnf

# (a) sequential Maxwell / Hamilton-Crosser, Eq. (5)
knf = k_f * (k_s1 + 2 * k_f - 2 * phi1 * (k_f - k_s1)) / (k_s1 + 2 * k_f + phi1 * (k_f - k_s1))
khnf = knf * (k_s2 + 2 * knf - 2 * phi2 * (knf - k_s2)) / (k_s2 + 2 * knf + phi2 * (knf - k_s2))
k_seq = khnf * (k_s3 + 2 * khnf - 2 * phi3 * (khnf - k_s3)) / (k_s3 + 2 * khnf + phi3 * (khnf - k_s3))

# (b) one-step Maxwell with volume-weighted effective particle conductivity
phitot = phi1 + phi2 + phi3
k_eff = (phi1 * k_s1 + phi2 * k_s2 + phi3 * k_s3) / phitot
k_alt = k_f * (k_eff + 2 * k_f - 2 * phitot * (k_f - k_eff)) / (k_eff + 2 * k_f + phitot * (k_f - k_eff))

Pr_seq = mu_thnf * cp_thnf / k_seq
Pr_alt = mu_thnf * cp_thnf / k_alt

print(f"k_thnf (sequential Maxwell/Hamilton-Crosser) = {k_seq:.6f} W/(m K)")
print(f"k_thnf (one-step, volume-weighted particle)  = {k_alt:.6f} W/(m K)")
print(f"pct change in k_thnf: {100*(k_alt-k_seq)/k_seq:+.3f}%")
print(f"Pr_thnf (sequential) = {Pr_seq:.6f}")
print(f"Pr_thnf (alt)        = {Pr_alt:.6f}")
print(f"pct change in Pr_thnf: {100*(Pr_alt-Pr_seq)/Pr_seq:+.3f}%")


def momentum_rhs(eta, y, Kp):
    F, Fp, G, Gp, H = y
    denom = 1 - 2 * Kp * F
    Fpp = (F**2 - G**2 + H * Fp - Kp * (Fp**2 - Gp**2)) / denom
    Gpp = (2 * F * G + H * Gp - 2 * Kp * Fp * Gp) / denom
    return np.vstack([Fp, Fpp, Gp, Gpp, -2 * F])


def momentum_bc(ya, yb, lam, S):
    return np.array([ya[0] - lam, ya[2] - 1.0, ya[4] + S, yb[0], yb[2]])


eta = np.linspace(0, ETA_INF, 400)
y0 = np.zeros((5, eta.size))
y0[0] = lam * np.exp(-eta)
y0[2] = np.exp(-eta)
y0[4] = -S
sol = solve_bvp(lambda e, y: momentum_rhs(e, y, K),
                 lambda ya, yb: momentum_bc(ya, yb, lam, S),
                 eta, y0, tol=1e-10, max_nodes=20000)
assert sol.status == 0, "momentum base state did not converge"
print(f"\nmomentum base state: F'(0) = {sol.y[1,0]:.6f}, -G'(0) = {-sol.y[3,0]:.6f} "
      "(cross-check against Example 1 baseline)")

H_of_eta = interp1d(sol.x, sol.y[4], kind="cubic", fill_value="extrapolate")


def energy_rhs(eta, y, Pr):
    theta, thetap = y
    return np.vstack([thetap, Pr * H_of_eta(eta) * thetap])


def energy_bc(ya, yb):
    return np.array([ya[0] - 1.0, yb[0]])


def solve_energy(Pr, n=800):
    eta = np.linspace(0, ETA_INF, n)
    y0 = np.zeros((2, eta.size))
    y0[0] = np.exp(-eta)
    y0[1] = -np.exp(-eta)
    return solve_bvp(lambda e, y: energy_rhs(e, y, Pr), energy_bc, eta, y0, tol=1e-10, max_nodes=20000)


sol_seq = solve_energy(Pr_seq)
sol_alt = solve_energy(Pr_alt)
tp_seq = -sol_seq.y[1, 0]
tp_alt = -sol_alt.y[1, 0]
print(f"\n-theta'(0) (sequential closure) = {tp_seq:.6f}")
print(f"-theta'(0) (one-step closure)   = {tp_alt:.6f}")
print(f"pct change in -theta'(0): {100*(tp_alt-tp_seq)/tp_seq:+.3f}%")
print()
print("Reported in Section 9.3 (Sensitivity to the ternary-mixing-rule")
print("correlations) of Main.tex: k_thnf changes by +0.53%, Pr_thnf by")
print("-0.52%, and -theta'(0) by -0.43% -- smaller than the +1.49% found")
print("for the viscosity correlation check in the same section.")
