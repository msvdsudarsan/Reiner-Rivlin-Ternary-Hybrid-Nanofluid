import sympy as sp

r, z = sp.symbols('r z', positive=True)
Omega, nu = sp.symbols('Omega nu', positive=True)
zeta = sp.sqrt(Omega/nu)
eta = z*zeta

F = sp.Function('F')
G = sp.Function('G')
H = sp.Function('H')

u = r*Omega*F(eta)
v = r*Omega*G(eta)
w = sp.sqrt(nu*Omega)*H(eta)

# --- rate of strain tensor (physical components), axisymmetric ---
e_rr   = sp.diff(u, r)
e_thth = u/r
e_zz   = sp.diff(w, z)
e_rz   = sp.Rational(1,2)*(sp.diff(u, z) + sp.diff(w, r))
e_rth  = sp.Rational(1,2)*(sp.diff(v, r) - v/r)
e_thz  = sp.Rational(1,2)*sp.diff(v, z)

E = sp.Matrix([[e_rr, e_rth, e_rz],
               [e_rth, e_thth, e_thz],
               [e_rz, e_thz, e_zz]])

A1 = 2*E
A1sq = sp.simplify(A1*A1)

print("A1^2 =")
sp.pprint(A1sq)
