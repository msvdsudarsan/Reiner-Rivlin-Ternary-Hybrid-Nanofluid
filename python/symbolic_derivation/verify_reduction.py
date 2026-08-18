import sympy as sp

eta = sp.symbols('eta')
K = sp.symbols('K')
F0=sp.Function('F0')(eta); G0=sp.Function('G0')(eta); H0=sp.Function('H0')(eta)
Fp,Fpp = sp.diff(F0,eta), sp.diff(F0,eta,2)
Gp,Gpp = sp.diff(G0,eta), sp.diff(G0,eta,2)
Hp = sp.diff(H0,eta)

# momentum residuals (verified earlier)
eqF = Fpp + K*(2*F0*Fpp + 5*Fp**2 + 2*Fp*sp.diff(H0,eta,2) + 2*Fpp*Hp - Gp**2) - F0**2 + G0**2 - H0*Fp
eqG = Gpp + K*(2*F0*Gpp + 6*Fp*Gp + 2*Gp*sp.diff(H0,eta,2) + 2*Gpp*Hp) - 2*F0*G0 - H0*Gp

# continuity: H' = -2F0  =>  H'' = -2F'
eqF2 = eqF.subs(sp.diff(H0,eta,2), -2*Fp).subs(Hp, -2*F0)
eqG2 = eqG.subs(sp.diff(H0,eta,2), -2*Fp).subs(Hp, -2*F0)

Fpp_sol = sp.solve(eqF2, Fpp)[0]
print("F'' =", sp.simplify(Fpp_sol))
Gpp_sol = sp.solve(eqG2.subs(Fpp, Fpp_sol), Gpp)
print("G'' =", sp.simplify(Gpp_sol[0]) if Gpp_sol else "check manually")

# direct solve of the 2x2 linear system in (Fpp,Gpp) simultaneously (more reliable)
Fpp_s, Gpp_s = sp.symbols('Fpp_s Gpp_s')
eqF3 = eqF2.subs(Fpp, Fpp_s).subs(Gpp, Gpp_s)
eqG3 = eqG2.subs(Fpp, Fpp_s).subs(Gpp, Gpp_s)
sol = sp.solve([eqF3, eqG3], [Fpp_s, Gpp_s])
print()
print("Simultaneous solve:")
print("F'' =", sp.simplify(sol[Fpp_s]))
print("G'' =", sp.simplify(sol[Gpp_s]))
