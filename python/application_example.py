"""Assumed laboratory disk: torque, entrainment and rim shear with the compatible closure."""
import numpy as np
from compatible_model import solve

R, OMEGA, RHO, ETA0, PSI2 = 0.02, 300 * 2 * np.pi / 60, 1000.0, 0.01, -1e-5
NU = ETA0 / RHO


def wall_values(K):
    s = solve(K, 0.0, 0.0, "compatible", None, 1e-10, L=30.0)
    return s.y[1, 0], s.y[3, 0], s.y[4, -1]


def main():
    K = PSI2 * OMEGA / ETA0
    out = {}
    for name, k in (("newtonian", 0.0), ("rr", K)):
        fp, gp, hinf = wall_values(k)
        torque = np.pi / 2 * ETA0 * R**4 * OMEGA * np.sqrt(OMEGA / NU) * abs(gp)
        flow = np.pi * R**2 * np.sqrt(NU * OMEGA) * abs(hinf)
        rate = R * OMEGA * np.sqrt(OMEGA / NU) * np.hypot(fp, gp)
        out[name] = (fp, gp, hinf, torque, flow, rate)
        print(name, "K=%.4f" % k, "F'(0)=%.4f G'(0)=%.4f H(inf)=%.4f" % (fp, gp, hinf),
              "M=%.4e Q=%.4e rate=%.1f" % (torque, flow, rate))
    print("torque change %.2f%%" % (100 * (out["rr"][3] / out["newtonian"][3] - 1)))
    print("entrainment change %.2f%%" % (100 * (out["rr"][4] / out["newtonian"][4] - 1)))
    rate = out["rr"][5]
    print("rim wall shear %.2f Pa, N2 %.2f Pa, local Wi %.2f" % (ETA0 * rate, PSI2 * rate**2, abs(PSI2) * rate / ETA0))


if __name__ == "__main__":
    main()
