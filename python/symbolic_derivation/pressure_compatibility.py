import sympy as sp
x,y,z,K=sp.symbols('x y z K',real=True)
F=sp.Function('F')(z);G=sp.Function('G')(z);H=sp.Function('H')(z)
# units: Omega=1, nu=mu/rho=1, rho=1, mu=1, mu_c=K  (K = mu_c*Omega/mu). H' = -2F
u=sp.Matrix([x*F-y*G, y*F+x*G, H])
X=sp.Matrix([x,y,z])
L=u.jacobian(X)                      # grad u
A1=L+L.T                             # Rivlin-Ericksen A1
S=A1+K*A1*A1                         # extra stress (mu=1, mu_c=K)
divS=sp.Matrix([sum(sp.diff(S[i,j],X[j]) for j in range(3)) for i in range(3)])
conv=L*u                              # (u.grad)u
subsH={sp.Derivative(H,z):-2*F}
def red(expr):
    e=expr.subs({y:0}).doit()
    return sp.simplify(e.subs(sp.Derivative(H,z,z),-2*sp.diff(F,z)).subs(sp.Derivative(H,z),-2*F))
mom=(divS-conv)       # = grad p
r=sp.symbols('r',positive=True)
Fz,Fzz,Gz=sp.diff(F,z),sp.diff(F,z,2),sp.diff(G,z)
# radial component at y=0 is x-component, divided by x=r
mx=red(mom[0].subs(x,r)); my=red(mom[1].subs(x,r)); mz=red(mom[2].subs(x,r))
print("radial/r  :",sp.expand(sp.simplify(mx/r)))
print("azimuth/r :",sp.expand(sp.simplify(my/r)))
print("z-comp    :",sp.expand(mz))
# Our pressure: p = P0(z) + r^2 K (F'^2+G'^2)  (rho=Omega=1);  d p/dr / r = 2K(F'^2+G'^2)
extra=2*K*(Fz**2+Gz**2)
print("radial/r minus dp/dr/r (ours):",sp.expand(sp.simplify(mx/r-extra)))
