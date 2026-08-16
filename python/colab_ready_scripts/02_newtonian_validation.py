import numpy as np
from scipy.integrate import solve_bvp

# Newtonian check: K=0, lambda=0, S=0 -- must match the classical
# textbook von Karman rotating disk values.

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

K, lam, S = 0.0, 0.0, 0.0
mesh = np.linspace(0, 20, 600)
sol = solve_bvp(lambda e,y: mom_rhs(e,y,K), lambda ya,yb: mom_bc(ya,yb,lam,S),
                 mesh, mom_guess(mesh, lam, S), tol=1e-11, max_nodes=200000, verbose=0)

print("=== Newtonian limit check (K=0, lambda=0, S=0) ===")
print(f"F'(0)    = {sol.y[1,0]:.6f}   (classical tabulated value: 0.5102)")
print(f"-G'(0)   = {-sol.y[3,0]:.6f}   (classical tabulated value: 0.6159)")
print(f"-H(inf)  = {-sol.y[4,-1]:.6f}   (classical tabulated value: 0.8845)")
print()
print("These are the well-known textbook von Karman rotating-disk values")
print("(no stretching, no cross-viscosity, no suction) -- if the three lines")
print("above are each close to their reference, the equations and solver")
print("are behaving correctly.")
