"""Rheometer reading to K, wall threshold and torque, with Monte Carlo uncertainty (synthetic inputs)."""
import numpy as np
from compatible_model import solve

R, RPM, RHO = 0.025, 60.0, 1000.0
ETA0, PSI2 = 1.0, -0.02
S_ETA, S_PSI = 0.05, 0.15
DRAWS, SEED = 50000, 20261116


def torque_coefficient_table():
    ks = np.linspace(-0.4, 0.0, 41)
    g = [abs(solve(float(k), 0.0, 0.0, "compatible", None, 1e-9, L=30.0).y[3, 0]) for k in ks]
    return ks, np.array(g)


def main():
    om = RPM * 2 * np.pi / 60
    ks, g = torque_coefficient_table()
    k0 = PSI2 * om / ETA0
    torque = lambda e, gp: np.pi / 2 * e * R**4 * om * np.sqrt(om * RHO / e) * gp
    print("Omega %.4f K %.5f lambda_w %.4f G'(0) %.4f torque %.4e" % (om, k0, 1 / (2 * k0), np.interp(k0, ks, g), torque(ETA0, np.interp(k0, ks, g))))
    rng = np.random.default_rng(SEED)
    e = rng.normal(ETA0, S_ETA * ETA0, DRAWS)
    p = rng.normal(PSI2, S_PSI * abs(PSI2), DRAWS)
    k = p * om / e
    ok = (k > ks[0]) & (k < 0)
    m = torque(e[ok], np.interp(k[ok], ks, g))
    print("retained %.3f" % ok.mean())
    print("K 95%%: %s" % np.percentile(k[ok], [2.5, 97.5]))
    print("lambda_w 95%%: %s" % np.percentile(1 / (2 * k[ok]), [2.5, 97.5]))
    print("torque 95%%: %s" % np.percentile(m, [2.5, 97.5]))


if __name__ == "__main__":
    main()
