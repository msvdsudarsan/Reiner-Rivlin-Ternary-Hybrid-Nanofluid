"""Inviscid stationary crossflow (type I) diagnostic of Gregory, Stuart & Walker (1955).

In the frame rotating with the disk the mean velocity is (F, G-1). A stationary
inviscid disturbance with wave angle eps sees the effective profile

    U_eps = F cos(eps) + (G - 1) sin(eps),

and the neutral inviscid stationary mode requires U_eps and U_eps'' to vanish at the
same interior point eta_c. Eliminating eps:

    F(eta_c) G''(eta_c) + F''(eta_c) [1 - G(eta_c)] = 0,   tan(eps) = F / (1 - G).

At leading order in the local Reynolds number the Reiner-Rivlin stress (viscous and
cross-viscous) enters the disturbance equations at O(1/R), so the fluid affects this
mode only through the mean profile.

This is a NECESSARY condition for a neutral inviscid stationary mode, not a viscous
stability threshold. It does not give the critical Reynolds number.

Expected output:
    Newtonian disk (K=0, lambda=0, S=0): eta_c = 1.458, eps = 13.22 deg  (GSW: 1.46, 13.2)
    then the critical points along K=+0.3 and K=-0.3 at S=1.5, as in the manuscript table.
"""
import numpy as np
from rr_solver import solve


def base(K, lam, S, einf=20.0):
    sol = solve(None, K, 0.0, S, einf=einf); l = 0.0
    step = 0.02 if lam > 0 else -0.02
    while abs(l - lam) > 1e-9:
        nxt = lam if abs(lam - l) < abs(step) else l + step
        s = solve(sol, K, nxt, S, einf=einf)
        if s.status == 0: sol, l = s, nxt
        else:
            step /= 2
            if abs(step) < 1e-7: return None
    return sol


def gsw_points(sol, K, eta_max=8.0):
    e = np.linspace(1e-3, eta_max, 6000)
    F, Fp, G, Gp, H = sol.sol(e)
    D = 1 - 2 * K * F
    Fpp = (F**2 - G**2 + H * Fp - K * (Fp**2 - Gp**2)) / D
    Gpp = (2 * F * G + H * Gp - 2 * K * Fp * Gp) / D
    g = F * Gpp + Fpp * (1 - G)
    out = []
    for i in range(len(e) - 1):
        if g[i] * g[i + 1] < 0:
            ec = e[i] - g[i] * (e[i + 1] - e[i]) / (g[i + 1] - g[i])
            y = sol.sol(np.array([ec]))
            out.append((ec, np.degrees(np.arctan2(y[0][0], 1 - y[2][0]))))
    return out


if __name__ == "__main__":
    r = gsw_points(base(0.0, 0.0, 0.0), 0.0)
    print("Newtonian disk: " + "  ".join("eta_c=%.3f eps=%.2f deg" % x for x in r))
    for K, lams in [(0.3, (0.0, -0.4469, -1.0, -2.0, -3.0, -4.0)),
                    (-0.3, (0.0, -0.3, -0.6036, -1.0, -1.6))]:
        print("\nK=%+.1f, S=1.5" % K)
        for lam in lams:
            r = gsw_points(base(K, lam, 1.5), K)
            print("  lambda=%7.4f  %s" % (lam, "  ".join("(%.3f, %.2f deg)" % x for x in r) if r else "none"))
