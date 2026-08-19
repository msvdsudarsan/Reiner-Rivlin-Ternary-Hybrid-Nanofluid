"""
derive_species_sign_check.py

Independent symbolic (sympy) re-derivation of the species perturbation
equation (linPhi in Main.tex), added in the V33 revision pass after an
independent AI peer review flagged a sign inconsistency between the
printed equation and the standard linearization procedure.

Method: the SAME substitution-and-linearize procedure used for the
thermal equation is applied first to the (already-published, unchanged,
and here re-confirmed correct) thermal case, to validate the method
itself; then to the species equation. If the method reproduces the
correct thermal equation, the species result it produces is trusted.

Result: the correct linearized species equation is

    Phi'' - Sc*H0*Phi' - Sc*phi0'*Hf - Sc*beta*Phi + gamma*Phi = 0

i.e. the reaction term enters as -Sc*beta*Phi, not +Sc*beta*Phi as an
earlier revision of the manuscript printed. This also corrects the
far-field threshold from gamma_thr = -Sc*beta to gamma_thr = +Sc*beta.
The Python implementation in stability.py was found, on inspection, to
already use the correct sign; only the printed equation in Main.tex
and the analytical argument built on it needed correcting. The
Chebyshev spectral script (10_chebyshev_species_search.py) DID have
the same sign bug as the old printed equation, and has been corrected
and re-run alongside this check.
"""
import sympy as sp

eta, tau, gamma = sp.symbols('eta tau gamma')
Sc, Pr, beta = sp.symbols('Sc Pr beta')

theta0 = sp.Function('theta0')(eta)
Theta = sp.Function('Theta')(eta)
H0 = sp.Function('H0')(eta)
Hf = sp.Function('Hf')(eta)
phi0 = sp.Function('phi0')(eta)
Phi = sp.Function('Phi')(eta)

eps = sp.exp(-gamma * tau)
e = sp.symbols('e')  # bookkeeping symbol standing in for eps, to extract the O(eps^1) term


def linearize(base, forcing_rhs_builder):
    """base + eps*perturbation substituted into base_tau = RHS(base+eps*pert);
    returns the coefficient of e^1 after expanding and substituting eps->e."""
    full = base
    resid = sp.diff(full, tau) - forcing_rhs_builder(full)
    resid = sp.expand(resid.subs(eps, e))
    poly = sp.Poly(resid, e)
    return sp.simplify(poly.coeff_monomial(e))


# ---- Step 1: validate the method against the thermal equation (known correct) ----
theta_full = theta0 + eps * Theta
H_full = H0 + eps * Hf


def thermal_rhs(field):
    # steady: theta'' = Pr*H*theta'  =>  unsteady: theta_tau = theta'' - Pr*H*theta'
    return sp.diff(field, eta, 2) - Pr * H_full * sp.diff(field, eta)


coeff_thermal = linearize(theta_full, thermal_rhs)
print("Thermal linearization (method-validation step):")
print(" ", coeff_thermal, "= 0")
print("  Matches manuscript Eq. (linTheta):",
      sp.simplify(coeff_thermal - (Pr * H0 * sp.diff(Theta, eta) + Pr * Hf * sp.diff(theta0, eta)
                                    - gamma * Theta - sp.diff(Theta, eta, 2))) == 0)
print()

# ---- Step 2: apply the SAME validated method to the species equation ----
phi_full = phi0 + eps * Phi


def species_rhs(field):
    # steady: phi'' = Sc*H*phi' + Sc*beta*phi
    #  =>  unsteady: phi_tau = phi'' - Sc*H*phi' - Sc*beta*phi
    return sp.diff(field, eta, 2) - Sc * H_full * sp.diff(field, eta) - Sc * beta * field


coeff_species = linearize(phi_full, species_rhs)
print("Species linearization (from the same, now-validated, method):")
print(" ", coeff_species, "= 0")

correct_form = (Sc * H0 * sp.diff(Phi, eta) + Sc * Hf * sp.diff(phi0, eta)
                 + Sc * beta * Phi - gamma * Phi - sp.diff(Phi, eta, 2))
old_wrong_form = (Sc * H0 * sp.diff(Phi, eta) + Sc * Hf * sp.diff(phi0, eta)
                   - Sc * beta * Phi - gamma * Phi - sp.diff(Phi, eta, 2))

print("Matches CORRECTED manuscript Eq. (linPhi) [-Sc*beta*Phi term, i.e. +Sc*beta*Phi here after sign flip]:",
      sp.simplify(coeff_species - correct_form) == 0)
print("Matches OLD (pre-V33) printed equation [+Sc*beta*Phi term, i.e. -Sc*beta*Phi here after sign flip]:",
      sp.simplify(coeff_species - old_wrong_form) == 0)
print()
print("(Note: coeff_species is the linearization written as 'stuff - Phi'' = 0';")
print(" the corrected manuscript equation Phi''-Sc*H0*Phi'-Sc*Hf*phi0'-Sc*beta*Phi+gamma*Phi=0")
print(" is the same equation multiplied by -1, i.e. Sc*beta*Phi+Sc*H0*Phi'+Sc*Hf*phi0'-gamma*Phi-Phi''=0,")
print(" which is exactly what 'correct_form' encodes above.)")
