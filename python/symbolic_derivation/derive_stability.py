import sympy as sp

eta, tau, K, gamma, eps, Pr, Sc, beta = sp.symbols('eta tau K gamma epsilon Pr Sc beta')

F0 = sp.Function('F0')(eta); G0 = sp.Function('G0')(eta); H0 = sp.Function('H0')(eta)
th0 = sp.Function('theta0')(eta); ph0 = sp.Function('phi0')(eta)
Ff = sp.Function('Ff')(eta); Gf = sp.Function('Gf')(eta); Hf = sp.Function('Hf')(eta)
Th = sp.Function('Theta')(eta); Ph = sp.Function('Phi')(eta)

# perturbed fields
F = F0 + eps*sp.exp(-gamma*tau)*Ff
G = G0 + eps*sp.exp(-gamma*tau)*Gf
H = H0 + eps*sp.exp(-gamma*tau)*Hf
th = th0 + eps*sp.exp(-gamma*tau)*Th
ph = ph0 + eps*sp.exp(-gamma*tau)*Ph

d = lambda expr: sp.diff(expr, eta)

def RHS_F(F,G,H):
    return sp.diff(F,eta,2) + K*(2*F*sp.diff(F,eta,2) + 5*sp.diff(F,eta)**2
            + 2*sp.diff(F,eta)*sp.diff(H,eta,2) + 2*sp.diff(F,eta,2)*sp.diff(H,eta)
            - sp.diff(G,eta)**2) - F**2 + G**2 - H*sp.diff(F,eta)

def RHS_G(F,G,H):
    return sp.diff(G,eta,2) + K*(2*F*sp.diff(G,eta,2) + 6*sp.diff(F,eta)*sp.diff(G,eta)
            + 2*sp.diff(G,eta)*sp.diff(H,eta,2) + 2*sp.diff(G,eta,2)*sp.diff(H,eta)) \
            - 2*F*G - H*sp.diff(G,eta)

def RHS_th(H, th):
    return sp.diff(th,eta,2) - Pr*H*sp.diff(th,eta)

def RHS_ph(H, ph):
    return sp.diff(ph,eta,2) - Sc*H*sp.diff(ph,eta) - Sc*beta*ph

eqF = sp.Eq(sp.diff(F,tau), RHS_F(F,G,H))
eqG = sp.Eq(sp.diff(G,tau), RHS_G(F,G,H))
eqTh = sp.Eq(sp.diff(th,tau), RHS_th(H, th))
eqPh = sp.Eq(sp.diff(ph,tau), RHS_ph(H, ph))

def order_eps1(eq):
    lhs_minus_rhs = eq.lhs - eq.rhs
    series = sp.series(lhs_minus_rhs, eps, 0, 2).removeO()
    coeff1 = sp.diff(series, eps).subs(eps,0)
    return sp.simplify(coeff1 / sp.exp(-gamma*tau))

lin_F = order_eps1(eqF)
lin_G = order_eps1(eqG)
lin_Th = order_eps1(eqTh)
lin_Ph = order_eps1(eqPh)

print("=== Linearized F-perturbation equation (=0) ===")
sp.pprint(sp.expand(lin_F))
print()
print("=== Linearized G-perturbation equation (=0) ===")
sp.pprint(sp.expand(lin_G))
print()
print("=== Linearized Theta-perturbation equation (=0) ===")
sp.pprint(sp.expand(lin_Th))
print()
print("=== Linearized Phi-perturbation equation (=0) ===")
sp.pprint(sp.expand(lin_Ph))
