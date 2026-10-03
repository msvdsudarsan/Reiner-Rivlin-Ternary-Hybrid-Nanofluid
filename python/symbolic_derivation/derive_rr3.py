import sympy as sp

r = sp.symbols('r', positive=True)
Omega, nu, mu, muc, rho = sp.symbols('Omega nu mu mu_c rho', positive=True)
eta = sp.symbols('eta', real=True)
zeta = sp.sqrt(Omega/nu)

F = sp.Function('F')(eta)
G = sp.Function('G')(eta)
H = sp.Function('H')(eta)

u = r*Omega*F
v = r*Omega*G
w = sp.sqrt(nu*Omega)*H

def ddr(expr):
    return sp.diff(expr, r)

def ddz(expr):
    return zeta*sp.diff(expr, eta)

# --- rate-of-strain tensor, physical components (axisymmetric, d/dtheta = 0) ---
e_rr   = sp.simplify(ddr(u))
e_thth = sp.simplify(u/r)
e_zz   = sp.simplify(ddz(w))
e_rz   = sp.simplify(sp.Rational(1,2)*(ddz(u) + ddr(w)))
e_rth  = sp.simplify(sp.Rational(1,2)*(ddr(v) - v/r))
e_thz  = sp.simplify(sp.Rational(1,2)*ddz(v))

print("e_rr   =", e_rr)
print("e_thth =", e_thth)
print("e_zz   =", e_zz)
print("e_rz   =", e_rz)
print("e_rth  =", e_rth, "  (should be 0)")
print("e_thz  =", e_thz)

E = sp.Matrix([[e_rr, e_rth, e_rz],
               [e_rth, e_thth, e_thz],
               [e_rz, e_thz, e_zz]])
A1 = 2*E
A1sq = sp.expand(A1*A1)
print()
print("A1^2 components:")
for i,ii in enumerate(['r','th','z']):
    for j,jj in enumerate(['r','th','z']):
        if j>=i:
            print(f"A1sq_{ii}{jj} =", sp.simplify(A1sq[i,j]))

with open('/home/claude/A1_E.pkl','wb') as f:
    import pickle
    pickle.dump({'E':E,'A1':A1,'A1sq':A1sq,'u':u,'v':v,'w':w}, f)
