import numpy as np
from scipy.integrate import solve_bvp

# Direct parameter-matched comparison against Puspanathan et al. (2024),
# Chinese Journal of Physics 88:198-211, Table 2, at their exact operating
# point K=1.5, S(their vartheta)=2.8, lambda(their epsilon)=-2.4.
# Their reported skin-friction coefficients convert (via their Eq. 24,
# identical in form to this manuscript's Corollary 1) to F'(0)=2.5190
# (first solution) and F'(0)=2.0836 (second solution). This script
# reports what natural continuation from lambda=0 gives at this same
# point -- reported honestly whether or not it matches either value.

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

K, S = 1.5, 2.8
target_lambda = -2.4
factor = 1 - 2*K*target_lambda

print(f"Puspanathan et al. (2024) Table 2, K={K}, S={S}, lambda={target_lambda}")
print(f"Conversion factor (1-2*K*lambda) = {factor}")
print(f"Their first solution:  F'(0) = 20.65572426/{factor} = {20.65572426/factor:.4f}")
print(f"Their second solution: F'(0) = 17.08537860/{factor} = {17.08537860/factor:.4f}")
print()

lam = 0.0
mesh = np.linspace(0, 15, 400)
sol = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam,S),
                 mesh, mom_guess(mesh, lam, S), tol=1e-9, max_nodes=8000, verbose=0)
dlam = -0.02
n = 0
while lam > target_lambda and n < 400:
    n += 1
    trial = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,max(lam+dlam,target_lambda),S),
                       sol.x, sol.y, tol=1e-9, max_nodes=8000, verbose=0)
    if trial.status != 0:
        dlam /= 2
        if abs(dlam) < 1e-6:
            break
        continue
    sol = trial
    lam = max(lam+dlam, target_lambda)

print("=== COMPARE ===")
print(f"Paper:  natural-continuation F'(0) at this point = 1.7190")
print(f"Yours:  natural-continuation F'(0) at this point = {sol.y[1,0]:.4f}")
print()
print("This value is not expected to closely match either of Puspanathan's")
print("branch values -- that mismatch is the point being documented, reported")
print("honestly in the manuscript's Limitations section as a genuine open item.")
