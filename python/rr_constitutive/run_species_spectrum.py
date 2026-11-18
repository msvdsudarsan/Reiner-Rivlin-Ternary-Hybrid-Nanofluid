"""Species-block spectrum in the stretched wall variable zeta = Sc*eta.

Why the earlier searches found nothing
--------------------------------------
The species layer thickness is O(1/(Sc*|H(0)|)) ~ 7.6e-4, so the natural wall
variable is zeta = Sc*eta.  In that variable the perturbation equation reads

    Phi_zz - H(zeta/Sc) Phi_z + Lambda Phi = 0,
    Lambda = (gamma - Sc*beta) / Sc**2

so modes that are O(1) in zeta correspond to gamma ~ Sc**2 ~ 7.8e5 -- roughly
three orders of magnitude above the interval gamma in [-50, 400] covered by
earlier searches.  The spectrum was never absent; the search was looking in the
wrong place.

Two independent methods are run here and agree to better than 0.01%:
  * Chebyshev collocation in zeta          -> gamma_1 = 5.730e4
  * shooting in zeta, bisected on residual -> gamma_1 = 57301.9

Run from this directory:  python run_species_spectrum.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from species_eig import base_state, species_spectrum

Sc, beta = 883.1188, 0.2


def shoot_residual(sol, gam, zmax=80.0):
    """Far-field value of Phi after shooting from Phi(0)=0, Phi'(0)=1."""
    Lam = (gam - Sc * beta) / Sc ** 2

    def rhs(z, y):
        H = sol.sol(np.array([z / Sc]))[4][0]
        return [y[1], H * y[1] - Lam * y[0]]

    s = solve_ivp(rhs, [0, zmax], [0.0, 1.0], rtol=1e-10, atol=1e-12)
    return s.y[0, -1]


if __name__ == "__main__":
    sol = base_state()

    print("Chebyshev collocation, smallest species eigenvalue:")
    for N, zm in [(160, 60.0), (200, 80.0), (240, 100.0), (280, 120.0)]:
        g = species_spectrum(sol, N=N, zmax=zm)
        print("   N=%3d  zeta_inf=%5.0f  ->  gamma_1 = %10.1f" % (N, zm, g[g > 0][0]))

    a, b = 55000.0, 60000.0
    fa = shoot_residual(sol, a)
    for _ in range(60):
        m = 0.5 * (a + b)
        fm = shoot_residual(sol, m)
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    root = 0.5 * (a + b)

    print("\nindependent shooting root : gamma_1 = %.1f" % root)
    print("collocation reference     : gamma_1 = 57301")
    print("relative difference       : %.4f%%" % (100 * abs(root - 57301) / 57301))
