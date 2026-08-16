import sympy as sp, pickle

r = sp.symbols('r', positive=True)
Omega, nu, mu, rho, K = sp.symbols('Omega nu mu rho K', positive=True)
eta = sp.symbols('eta', real=True)
zeta = sp.sqrt(Omega/nu)
muc = K*mu/Omega   # dimensionless cross-viscosity parameter K = mu_c*Omega/mu

with open('/home/claude/A1_E.pkl','rb') as f:
    d = pickle.load(f)
u, v, w = d['u'], d['v'], d['w']

def ddr(expr): return sp.diff(expr, r)
def ddz(expr): return zeta*sp.diff(expr, eta)

# ---- continuity check: (1/r) d(ru)/dr + dw/dz = 0 ?
cont = sp.simplify(sp.diff(r*u, r)/r + ddz(w))
print("continuity residual (want it in form const*(2F+H')):", cont)

F = sp.Function('F')(eta); G = sp.Function('G')(eta); H = sp.Function('H')(eta)
Fp, Fpp = sp.Derivative(F,eta), sp.Derivative(F,eta,2)
Gp, Gpp = sp.Derivative(G,eta), sp.Derivative(G,eta,2)
Hp, Hpp = sp.Derivative(H,eta), sp.Derivative(H,eta,2)

E = sp.Matrix([[u.diff(r), sp.Rational(1,2)*(v.diff(r)-v/r), sp.Rational(1,2)*(u.diff(0*r+1,evaluate=False) if False else 0)]]) # placeholder, unused

# rebuild strain tensor cleanly (as before)
e_rr   = ddr(u)
e_thth = u/r
e_zz   = ddz(w)
e_rz   = sp.Rational(1,2)*(ddz(u) + ddr(w))
e_rth  = sp.Rational(1,2)*(ddr(v) - v/r)
e_thz  = sp.Rational(1,2)*ddz(v)

Emat = sp.Matrix([[e_rr, e_rth, e_rz],[e_rth, e_thth, e_thz],[e_rz, e_thz, e_zz]])
A1 = 2*Emat
A1sq = sp.expand(A1*A1)

S = mu*A1 + muc*A1sq
S_rr, S_rth, S_rz = S[0,0], S[0,1], S[0,2]
S_thth, S_thz = S[1,1], S[1,2]
S_zz = S[2,2]

divS_r  = ddr(S_rr) + ddz(S_rz) + (S_rr - S_thth)/r
divS_th = ddr(S_rth) + ddz(S_thz) + 2*S_rth/r

a_r  = u*ddr(u) + w*ddz(u) - v**2/r
a_th = u*ddr(v) + w*ddz(v) + u*v/r

rmom = sp.simplify((rho*a_r - divS_r)/(r*rho*Omega**2))
thmom = sp.simplify((rho*a_th - divS_th)/(r*rho*Omega**2))

print()
print("=== r-momentum / (r*rho*Omega^2) = 0  ===")
print(sp.nsimplify(sp.expand(rmom)))
print()
print("=== theta-momentum / (r*rho*Omega^2) = 0 ===")
print(sp.nsimplify(sp.expand(thmom)))
