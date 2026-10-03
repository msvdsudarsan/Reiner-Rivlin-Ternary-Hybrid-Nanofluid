"""Is the regularity boundary an artefact of the Reiner-Rivlin N1 = 0 closure?

Compares, symbolically, the radial momentum equation produced by

    Reiner-Rivlin   tau = -p I + mu A1 + mu_c A1^2          (N1 = 0)
    second-order    tau = -p I + mu A1 + alpha_1 A2 + mu_c A1^2   (N1 != 0)

under the von Karman similarity ansatz, and reports the coefficients of the two
highest derivatives of F.

Result (reproduced by running this file):

    Reiner-Rivlin   coeff F''' = 0
                    coeff F''  = mu Omega^2 r (1 - 2 K F) / nu      [K = mu_c Omega / mu]

    second-order    coeff F''' = alpha_1 H Omega^3 r / nu   (nonzero -> third order)
                    coeff F''  = identical to Reiner-Rivlin; the alpha_1 terms
                                 cancel against the continuity identity H' = -2F

So the singular factor is produced by the A1^2 term alone, i.e. by the second
normal-stress response, and is untouched by the term that generates N1.
"""
import sympy as sp

r, Om, nu, mu, muc, a1 = sp.symbols('r Omega nu mu mu_c alpha_1', positive=True)
e = sp.symbols('eta'); k = sp.sqrt(Om / nu)
F, G, H = sp.Function('F'), sp.Function('G'), sp.Function('H')

u = r * Om * F(e); v = r * Om * G(e); w = sp.sqrt(nu * Om) * H(e)
ddz = lambda x: sp.diff(x, e) * k

L = sp.Matrix([[sp.diff(u, r), -v / r, ddz(u)],
               [sp.diff(v, r),  u / r, ddz(v)],
               [sp.diff(w, r),      0, ddz(w)]])
A1 = sp.simplify(L + L.T)
A2 = sp.simplify(u * sp.diff(A1, r) + w * sp.diff(A1, e) * k + A1 * L + L.T * A1)

divr = lambda T: sp.diff(T[0, 0], r) + ddz(T[0, 2]) + (T[0, 0] - T[1, 1]) / r

F0, F1, F2, F3 = sp.symbols('F0 F1 F2 F3')
H0, H1 = sp.symbols('H0 H1')
sub = {sp.Derivative(F(e), (e, 3)): F3, sp.Derivative(F(e), (e, 2)): F2,
       sp.Derivative(F(e), e): F1, F(e): F0,
       sp.Derivative(H(e), e): H1, H(e): H0}

if __name__ == "__main__":
    out = {}
    for name, T in [("Reiner-Rivlin", mu * A1 + muc * (A1 * A1)),
                    ("second-order ", mu * A1 + a1 * A2 + muc * (A1 * A1))]:
        d = sp.expand(divr(T).doit()).subs(sub)
        c3 = sp.simplify(sp.expand(d).coeff(F3))
        c2 = sp.factor(sp.simplify(sp.expand(d).coeff(F2)).subs({H1: -2 * F0}))
        out[name] = c2
        print("%s  coeff F''' = %s" % (name, c3))
        print("%s  coeff F''  = %s   (after H' = -2F)" % (name, c2))
    diff = sp.simplify(out["second-order "] - out["Reiner-Rivlin"])
    print("\ndifference in the F'' coefficient:", diff)
    assert diff == 0, "alpha_1 did not cancel"
    print("the alpha_1 contributions cancel exactly: the critical factor is set by mu_c alone")
