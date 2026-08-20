import sympy as sp, pickle

r = sp.symbols('r', positive=True)
Omega, nu, mu, muc, rho = sp.symbols('Omega nu mu mu_c rho', positive=True)
eta = sp.symbols('eta', real=True)
zeta = sp.sqrt(Omega/nu)

with open('/home/claude/A1_E.pkl','rb') as f:
    d = pickle.load(f)
u, v, w = d['u'], d['v'], d['w']
A1sq = d['A1sq']
A1 = d['A1']

def ddr(expr):
    return sp.diff(expr, r)
def ddz(expr):
    return zeta*sp.diff(expr, eta)

# extra-stress tensor S = mu*A1 + mu_c*A1^2  (physical components r,th,z)
S = mu*A1 + muc*A1sq
S_rr, S_rth, S_rz = S[0,0], S[0,1], S[0,2]
S_thth, S_thz = S[1,1], S[1,2]
S_zz = S[2,2]

# divergence of symmetric tensor S in cylindrical coords (axisymmetric, d/dtheta=0)
divS_r  = ddr(S_rr) + ddz(S_rz) + (S_rr - S_thth)/r
divS_th = ddr(S_rth) + ddz(S_thz) + 2*S_rth/r
divS_z  = ddr(S_rz) + ddz(S_zz) + S_rz/r

divS_r  = sp.simplify(divS_r)
divS_th = sp.simplify(divS_th)
divS_z  = sp.simplify(divS_z)

# convective inertia
a_r  = u*ddr(u) + w*ddz(u) - v**2/r
a_th = u*ddr(v) + w*ddz(v) + u*v/r
a_z  = u*ddr(w) + w*ddz(w)

a_r  = sp.simplify(a_r)
a_th = sp.simplify(a_th)
a_z  = sp.simplify(a_z)

print("=== inertia terms ===")
print("a_r  =", a_r)
print("a_th =", a_th)
print("a_z  =", a_z)
print()
print("=== div S terms ===")
print("divS_r  =", divS_r)
print("divS_th =", divS_th)
print("divS_z  =", divS_z)

# r-momentum: rho*a_r = -dp/dr + divS_r   (assume p = p(z) only, i.e. dp/dr = 0; check consistency)
rmom = sp.simplify(rho*a_r - divS_r)
thmom = sp.simplify(rho*a_th - divS_th)
print()
print("=== rho*a_r - divS_r  (should reduce cleanly, proportional to r) ===")
print(sp.simplify(rmom/r))
print()
print("=== rho*a_th - divS_th (should reduce cleanly, proportional to r) ===")
print(sp.simplify(thmom/r))
print()
print("=== z-momentum: rho*a_z - divS_z =  -dp/dz  (defines P(eta)) ===")
zmom = sp.simplify(rho*a_z - divS_z)
print(zmom)
