import numpy as np
from scipy.integrate import solve_bvp

# Comparison with Puspanathan et al. (2024), Chinese Journal of Physics 88:198-211,
# Table 2, at K = 1.5 and lambda (their epsilon) = -2.4. Their suction parameter
# (vartheta) is half of the suction parameter S used here, so vartheta = 2.8
# corresponds to S = 5.6. Their Table 2 coefficients convert through
# (1 - 2*K*lambda) to F'(0) = 2.5190 (first) and 2.0836 (second).
#
# Part A: their own variables, multi-start shooting (vartheta = 2.8, S = 5.6).
# Part B: natural continuation in lambda from lambda = 0 in the present
#         variables, at S = 5.6 and, for reference, at S = 2.8, which is a
#         different operating point (their vartheta = 1.4).

K, LAM = 1.5, -2.4
fac = 1 - 2*K*LAM
print(f"Printed values converted: first {20.65572426/fac:.5f}, second {17.08537860/fac:.5f}")

def part_a(th):
    def rhs(x, y):
        f, f1, f2, g, g1 = y
        d = 1 - 2*K*f1
        return np.vstack([f1, f2, (f1**2-2*f*f2-g**2-K*f2**2+K*g1**2)/d, g1,
                          (2*f1*g-2*f*g1-2*K*f2*g1)/d])
    def bc(a, b):
        return np.array([a[0]-th, a[1]-LAM, a[3]-1, b[1], b[3]])
    found = set()
    for Z in (12, 20):
        for c in np.linspace(0.5, 6, 23):
            for e in (0.3, 1.0):
                x = np.linspace(0, Z, 500); y = np.zeros((5, x.size))
                y[0] = th+LAM*(1-np.exp(-x)); y[1] = LAM*np.exp(-x)
                y[2] = c*np.exp(-x); y[3] = np.exp(-e*x); y[4] = -e*np.exp(-e*x)
                try:
                    s = solve_bvp(rhs, bc, x, y, tol=1e-8, max_nodes=100000)
                except Exception:
                    continue
                if s.status == 0 and (1-2*K*s.y[1]).min() > 0:
                    found.add(round(s.y[2, 0]*fac, 4))
    return sorted(found)

print("Part A (vartheta = 2.8, S = 5.6): converged Cf values =", part_a(2.8))
print("  Cf/(1-2K lambda) gives F'(0); compare with the converted printed values above.")
print()

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


def part_b(S):
    lam = 0.0
    mesh = np.linspace(0, 15, 400)
    sol = solve_bvp(lambda e, y: mom_rhs(e, y, K), lambda ya, yb: mom_bc(ya, yb, lam, S),
                    mesh, mom_guess(mesh, lam, S), tol=1e-9, max_nodes=8000, verbose=0)
    dlam, n = -0.02, 0
    while lam > LAM and n < 400:
        n += 1
        trial = solve_bvp(lambda e, y: mom_rhs(e, y, K),
                          lambda ya, yb: mom_bc(ya, yb, max(lam+dlam, LAM), S),
                          sol.x, sol.y, tol=1e-9, max_nodes=8000, verbose=0)
        if trial.status != 0:
            dlam /= 2
            if abs(dlam) < 1e-6:
                break
            continue
        sol = trial
        lam = max(lam+dlam, LAM)
    return lam, sol.y[1, 0]

for S in (5.6, 2.8):
    lam, v = part_b(S)
    print(f"Part B, S = {S}: continuation reached lambda = {lam:.3f}, F'(0) = {v:.4f}")
print("S = 5.6 is the operating point of their Table 2. The S = 2.8 value is a different case.")
