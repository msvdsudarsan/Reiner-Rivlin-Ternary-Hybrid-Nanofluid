"""Illustrative shrinking-disk calculation with the full constitutive tensor (assumed inputs).

Fluid: Reiner-Rivlin, eta0 = 1 Pa s, rho = 970 kg/m^3, Psi2 = K*eta0/Omega with K = -0.1
(Omega = 10 rad/s, so Psi2 = -0.01 Pa s^2). Disk: R = 0.05 m, shrinking parameter lambda = -0.4,
suction S = 1.5, closure Pi = 0. None of the inputs is a measured value.
Run from this folder:  python3 shrinking_disk_example.py
"""
import os, sys
import numpy as np
from scipy.linalg import eig
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "rr_constitutive"))
from cheb_base import base_state_from_bvp
from stability_tracked import pencil

K, lam, S = -0.1, -0.4, 1.5
eta0, rho, Om, R = 1.0, 970.0, 10.0, 0.05
nu = eta0/rho; q = R*Om*np.sqrt(Om/nu); mu_c = K*eta0/Om

def state(Kv, N=100, L=25):
    return base_state_from_bvp(Kv, lam, S, N=N, einf=L)

eta, F, G, H, De, D2 = state(K); _, Fn, Gn, Hn, Den, _ = state(0.0)
Fp, Gp = (De@F)[0], (De@G)[0]; Gpn = (Den@Gn)[0]
print("min(1-2KF) = %.4f" % np.min(1 - 2*K*F))
print("wall gradients: F'(0)=%.6f  G'(0)=%.6f  H(inf)=%.6f (Newtonian G'(0)=%.6f)" % (Fp, Gp, H[-1], Gpn))
print("torque ratio RR/Newtonian - 1 = %.4f %%" % (100*(abs(Gp)*abs(1-2*K*lam)/abs(Gpn) - 1)))
# wall velocity-gradient tensor in (r, theta, z); F(0)=lambda, G(0)=1
Lm = np.array([[Om*lam, -Om, q*Fp], [Om, Om*lam, q*Gp], [0, 0, -2*Om*lam]])
A1 = Lm + Lm.T; T = eta0*A1 + mu_c*A1@A1
print("tau_thz = %.5f Pa (formula eta0*q*G'(0)*(1-2K*lam) = %.5f)" % (T[1, 2], eta0*q*Gp*(1-2*K*lam)))
print("tau_rr-tau_thth = %.4f Pa, tau_thth-tau_zz = %.4f Pa, tau_rr-tau_zz = %.4f Pa" % (T[0,0]-T[1,1], T[1,1]-T[2,2], T[0,0]-T[2,2]))
A, B = pencil(K, eta, F, G, H, De, D2)
w = eig(A, B, right=False); w = w[np.isfinite(w)]; w = np.sort(w[abs(w.imag) < 1e-8].real)
print("smallest real eigenvalues (|gamma|>1e-3):", w[abs(w) > 1e-3][:3])

# ---- component breakdown: Newtonian part and second-normal-stress part -------------------
Tn = eta0*A1; Tm = mu_c*A1@A1
print("A1 (Pa-free, 1/s) =\n", np.round(A1, 4))
print("A1^2 =\n", np.round(A1@A1, 2))
for nm, i, j in (("rr-thth", 0, 1), ("thth-zz", 1, 2), ("rr-zz", 0, 2)):
    print("%-8s total %.4f = Newtonian %.4f + mu_c*A1^2 part %.4f" % (nm, T[i,i]-T[j,j], Tn[i,i]-Tn[j,j], Tm[i,i]-Tm[j,j]))
print("closed form tau_thth-tau_zz = 6*eta0*Om*lam + mu_c*(-12*Om^2*lam^2 - q^2*F'(0)^2) = %.4f" %
      (6*eta0*Om*lam + mu_c*(-12*Om**2*lam**2 - q**2*Fp**2)))
# dimensional kinematics (u = Om r F, v = Om r G, w = sqrt(nu Om) H)
w0 = np.sqrt(nu*Om)
print("wall velocities at r=R: u_r = %.4f m/s, u_theta = %.4f m/s" % (Om*R*lam, Om*R))
print("axial velocity at the wall: %.5f m/s (suction S=%.1f); far field (entrainment) %.5f m/s" % (-S*w0, S, w0*H[-1]))
print("check: F(0)=%.4f G(0)=%.4f H(0)=%.4f" % (F[0], G[0], H[0]))
