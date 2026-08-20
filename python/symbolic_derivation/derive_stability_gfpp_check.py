"""
derive_stability_gfpp_check.py

Independent symbolic (sympy) re-derivation of manuscript Eq. (linG) (the
azimuthal perturbation momentum equation), added in the V26 revision pass
while investigating the momentum-eigenvalue-near-boundary open item.

Method: build Eq. (linG) symbolically with ALL EIGHT of its K-bracket
terms exactly as printed in Main.tex, eliminate H0''=-2F0' and
Hf''=-2Ff' (the manuscript's own stated elimination, used identically
for Eq. linF), solve for Gf'', and symbolically subtract the formula
previously used in stability.py and stability_shooting.py. A nonzero,
non-simplifying-to-zero remainder is conclusive evidence of a coding
error, independent of any numerical run.

Result: the remainder is 2*K*(Ff*G0'' - Gf*F0''), a term that was
missing from the previously-used formula. The corrected formula
(now used throughout the codebase) is:

    Delta*Gfpp = 4*K*Ff*G0pp - 2*K*Gf*F0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p
                 - gamma*Gf + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p
"""
import sympy as sp

eta = sp.symbols('eta')
F0 = sp.Function('F0')(eta); G0 = sp.Function('G0')(eta); H0 = sp.Function('H0')(eta)
Ff = sp.Function('Ff')(eta); Gf = sp.Function('Gf')(eta); Hf = sp.Function('Hf')(eta)
K, gamma = sp.symbols('K gamma')

F0p, F0pp = F0.diff(eta), F0.diff(eta, 2)
G0p, G0pp = G0.diff(eta), G0.diff(eta, 2)
H0p, H0pp = H0.diff(eta), H0.diff(eta, 2)
Ffp, Ffpp = Ff.diff(eta), Ff.diff(eta, 2)
Gfp, Gfpp = Gf.diff(eta), Gf.diff(eta, 2)
Hfp, Hfpp = Hf.diff(eta), Hf.diff(eta, 2)

# Eq. (linF), exactly as printed (included for completeness/symmetry;
# independently confirmed to already match the coded Ffpp exactly).
eqF = (Ffpp + K * (2*F0*Ffpp + 2*F0pp*Ff + 10*F0p*Ffp + 2*F0p*Hfpp
                    + 2*H0pp*Ffp + 2*F0pp*Hfp + 2*H0p*Ffpp - 2*G0p*Gfp)
       - 2*F0*Ff + 2*G0*Gf - H0*Ffp - F0p*Hf + gamma*Ff)

# Eq. (linG), exactly as printed (all 8 K-bracket terms):
# K[ 2F0G'' + 2F0''G + 6F0'G' + 6G0'F' + 2G0'H'' + 2H0''G' + 2G0''H' + 2H0'G'' ]
eqG = (Gfpp + K * (2*F0*Gfpp + 2*F0pp*Gf + 6*F0p*Gfp + 6*G0p*Ffp + 2*G0p*Hfpp
                    + 2*H0pp*Gfp + 2*G0pp*Hfp + 2*H0p*Gfpp)
       - 2*F0*Gf - 2*G0*Ff - H0*Gfp - G0p*Hf + gamma*Gf)

subs = {H0pp: -2*F0p, Hfpp: -2*Ffp, H0p: -2*F0, Hfp: -2*Ff}
solF = sp.solve(sp.expand(eqF.subs(subs)), Ffpp)[0]
solG = sp.solve(sp.expand(eqG.subs(subs)), Gfpp)[0]
Delta = 1 - 2*K*F0

code_Ffpp = (2*K*Ff*F0pp - 2*K*F0p*Ffp + 2*K*G0p*Gfp - gamma*Ff
             + 2*F0*Ff - 2*G0*Gf + H0*Ffp + Hf*F0p)
code_Gfpp_OLD = (2*K*Ff*G0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p - gamma*Gf
                 + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p)
code_Gfpp_NEW = (4*K*Ff*G0pp - 2*K*Gf*F0pp - 2*K*F0p*Gfp - 2*K*Ffp*G0p - gamma*Gf
                 + 2*F0*Gf + 2*Ff*G0 + H0*Gfp + Hf*G0p)

print("Ffpp check (0 expected -- confirms this formula was always correct):")
print("  ", sp.simplify(solF*Delta - code_Ffpp))
print("Gfpp check against pre-V26 code (nonzero => bug confirmed):")
print("  ", sp.simplify(solG*Delta - code_Gfpp_OLD))
print("Gfpp check against V26-corrected code (0 expected):")
print("  ", sp.simplify(solG*Delta - code_Gfpp_NEW))
