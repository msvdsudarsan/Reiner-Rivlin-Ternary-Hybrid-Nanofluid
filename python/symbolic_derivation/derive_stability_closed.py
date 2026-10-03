import sympy as sp

eta, K, gamma = sp.symbols('eta K gamma')
F0=sp.Function('F0')(eta); G0=sp.Function('G0')(eta); H0=sp.Function('H0')(eta)
Ff=sp.Function('Ff')(eta); Gf=sp.Function('Gf')(eta); Hf=sp.Function('Hf')(eta)

F0p,F0pp = sp.diff(F0,eta), sp.diff(F0,eta,2)
G0p,G0pp = sp.diff(G0,eta), sp.diff(G0,eta,2)
H0p = sp.diff(H0,eta)
Ffp,Ffpp = sp.diff(Ff,eta), sp.diff(Ff,eta,2)
Gfp,Gfpp = sp.diff(Gf,eta), sp.diff(Gf,eta,2)
Hfp = sp.diff(Hf,eta)

# linearized eqs from previous verified derivation (equations = 0)
linF = ( -2*K*F0*Ffpp - 2*K*Ff*F0pp - 10*K*F0p*Ffp - 2*K*F0p*sp.diff(Hf,eta,2)
         - 2*K*sp.diff(H0,eta,2)*Ffp - 2*K*Ffp*sp.diff(H0,eta,2)  # placeholder fix below
       )
# Instead of re-typing by hand (error prone), just re-derive with sympy series expansion again, but now
# immediately substitute H0''=-2F0', Hf''=-2Ff', H0'=-2F0, Hf'=-2Ff, then solve for Ffpp,Gfpp.

tau, eps = sp.symbols('tau epsilon')
F0f=sp.Function('F0')(eta); G0f=sp.Function('G0')(eta); H0f=sp.Function('H0')(eta)
Fff=sp.Function('Ff')(eta); Gff=sp.Function('Gf')(eta); Hff=sp.Function('Hf')(eta)

F = F0f + eps*sp.exp(-gamma*tau)*Fff
G = G0f + eps*sp.exp(-gamma*tau)*Gff
H = H0f + eps*sp.exp(-gamma*tau)*Hff

def RHS_F(F,G,H):
    return sp.diff(F,eta,2) + K*(2*F*sp.diff(F,eta,2) + 5*sp.diff(F,eta)**2
            + 2*sp.diff(F,eta)*sp.diff(H,eta,2) + 2*sp.diff(F,eta,2)*sp.diff(H,eta)
            - sp.diff(G,eta)**2) - F**2 + G**2 - H*sp.diff(F,eta)

def RHS_G(F,G,H):
    return sp.diff(G,eta,2) + K*(2*F*sp.diff(G,eta,2) + 6*sp.diff(F,eta)*sp.diff(G,eta)
            + 2*sp.diff(G,eta)*sp.diff(H,eta,2) + 2*sp.diff(G,eta,2)*sp.diff(H,eta)) \
            - 2*F*G - H*sp.diff(G,eta)

eqF = sp.diff(F,tau) - RHS_F(F,G,H)
eqG = sp.diff(G,tau) - RHS_G(F,G,H)

def order_eps1(expr):
    s = sp.series(expr, eps, 0, 2).removeO()
    c1 = sp.diff(s, eps).subs(eps,0)
    return sp.simplify(c1/sp.exp(-gamma*tau))

linF = order_eps1(eqF)
linG = order_eps1(eqG)

# now substitute continuity identities: H0'=-2F0, H0''=-2F0', Hf'=-2Ff, Hf''=-2Ff'
subs_dict = {
    sp.diff(H0f,eta,2): -2*sp.diff(F0f,eta),
    sp.diff(H0f,eta): -2*F0f,
    sp.diff(Hff,eta,2): -2*sp.diff(Fff,eta),
    sp.diff(Hff,eta): -2*Fff,
}
linF2 = linF.subs(subs_dict)
linG2 = linG.subs(subs_dict)

Ffpp, Gfpp = sp.symbols('Ffpp Gfpp')
linF3 = linF2.subs(sp.diff(Fff,eta,2), Ffpp).subs(sp.diff(Gff,eta,2), Gfpp)
linG3 = linG2.subs(sp.diff(Fff,eta,2), Ffpp).subs(sp.diff(Gff,eta,2), Gfpp)

sol = sp.solve([linF3, linG3], [Ffpp, Gfpp])
print("Ff'' =")
sp.pprint(sp.simplify(sol[Ffpp]))
print()
print("Gf'' =")
sp.pprint(sp.simplify(sol[Gfpp]))
