import sympy as sp

r, z, t = sp.symbols('r z t', positive=True)
Omega, nu, mu, muc, rho = sp.symbols('Omega nu mu mu_c rho', positive=True)
eta = sp.symbols('eta')

F = sp.Function('F')
G = sp.Function('G')
H = sp.Function('H')

eta_expr = z*sp.sqrt(Omega/nu)

u = r*Omega*F(eta_expr)
v = r*Omega*G(eta_expr)
w = sp.sqrt(nu*Omega)*H(eta_expr)

# rate of strain (physical/orthonormal components), axisymmetric (no theta dependence)
e_rr = sp.diff(u, r)
e_thth = u/r
e_zz = sp.diff(w, z)
e_rz = sp.Rational(1,2)*(sp.diff(u,z) + sp.diff(w,r))
e_rth = sp.Rational(1,2)*(r*sp.diff(v/r, r))   # = (1/2)(dv/dr - v/r) since v depends on r explicitly and via eta(z) only -> dv/dr at fixed z
e_thz = sp.Rational(1,2)*sp.diff(v, z)

e_rr = sp.simplify(e_rr)
e_thth = sp.simplify(e_thth)
e_zz = sp.simplify(e_zz)
e_rz = sp.simplify(e_rz)
e_rth = sp.simplify(e_rth)
e_thz = sp.simplify(e_thz)

print("e_rr =", e_rr)
print("e_thth =", e_thth)
print("e_zz =", e_zz)
print("e_rz =", e_rz)
print("e_rth =", e_rth)
print("e_thz =", e_thz)
