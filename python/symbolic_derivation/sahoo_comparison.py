"""Exact cylindrical-coordinate audit of Sahoo (2009), Eqs. (9)-(11), (16a).

Constitutive law in Sahoo:
    tau = -p I + 2 mu e + 2 mu_c e.e, e = sym(grad v).
Set mu=Omega=nu=1 and L=mu_c*Omega/mu. Since A1=2e, the paper's
main-text law mu A1 + mu_c(main) A1^2 maps to L=2K.

The script derives all three momentum components from the cylindrical
stress divergence, imposes H'=-2F, and compares them with the printed
equations. It does not infer that Sahoo's numerical results are wrong.
"""
import sympy as sp
r, z, L = sp.symbols("r z L", real=True)
F, G, H, P = [sp.Function(n)(z) for n in "FGHP"]
u, v, w = r*F, r*G, H

# Cylindrical components of e=sym(grad v), axisymmetric velocity.
e = sp.zeros(3, 3)
e[0, 0] = sp.diff(u, r)
e[1, 1] = u/r
e[2, 2] = sp.diff(w, z)
e[0, 2] = e[2, 0] = sp.diff(u, z)/2
e[1, 2] = e[2, 1] = sp.diff(v, z)/2
e[0, 1] = e[1, 0] = (sp.diff(v, r)-v/r)/2

# Deviatoric stress divided by mu*Omega:
T = 2*e + 2*L*(e*e)
div_r = sp.diff(T[0, 0], r) + sp.diff(T[0, 2], z) + (T[0, 0]-T[1, 1])/r
div_phi = sp.diff(T[0, 1], r) + sp.diff(T[1, 2], z) + 2*T[0, 1]/r
div_z = sp.diff(T[0, 2], r) + T[0, 2]/r + sp.diff(T[2, 2], z)

inert_r = u*sp.diff(u, r) + w*sp.diff(u, z) - v**2/r
inert_phi = u*sp.diff(v, r) + w*sp.diff(v, z) + u*v/r
inert_z = u*sp.diff(w, r) + w*sp.diff(w, z)

cont = {sp.diff(H, z): -2*F, sp.diff(H, z, 2): -2*sp.diff(F, z)}
# Pressure p=-P(z); momentum residual is inertia-div(stress)+grad(p).
Rr = sp.factor(((inert_r-div_r)/r).subs(cont))
Rphi = sp.factor(((inert_phi-div_phi)/r).subs(cont))
Rz = sp.factor((inert_z-div_z-sp.diff(P, z)).subs(cont))

Fp, Fpp = sp.diff(F,z), sp.diff(F,z,2)
Gp, Gpp = sp.diff(G,z), sp.diff(G,z,2)

# Normalize radial/azimuthal equations to coefficient +1 in Newtonian limit.
our_radial = sp.expand(-Rr)
our_azimuthal = sp.expand(-Rphi)
our_axial = sp.expand(-Rz)

# Printed Sahoo (2009), Eqs. (9)-(11), with Mn=0.
sahoo_radial = Fpp-H*Fp-F**2+G**2-sp.Rational(1,2)*L*(Fp**2-3*Gp**2-2*F*Fpp)
sahoo_azimuthal = Gpp-H*Gp-2*F*G+L*(Fp*Gp+F*Gpp)
sahoo_axial = sp.diff(H,z,2)-H*sp.diff(H,z)-sp.Rational(7,2)*L*sp.diff(H,z)*sp.diff(H,z,2)+sp.diff(P,z)

# Exact normalized equations from this derivation.
expected_radial = (1-L*F)*Fpp-H*Fp-F**2+G**2+sp.Rational(1,2)*L*(Fp**2-Gp**2)
expected_azimuthal = (1-L*F)*Gpp+(L*Fp-H)*Gp-2*F*G
expected_axial = sp.diff(H,z,2)-H*sp.diff(H,z)+sp.Rational(7,2)*L*sp.diff(H,z)*sp.diff(H,z,2) \
                 + L*r**2*(Fp*Fpp+Gp*Gpp)+sp.diff(P,z)

assert sp.simplify(our_radial-expected_radial) == 0
assert sp.simplify(our_azimuthal-expected_azimuthal) == 0
assert sp.simplify(our_axial-expected_axial.subs(cont)) == 0

# Wall tractions in the same Sahoo scaling.
wall_r = sp.factor((T[0,2]/r).subs(cont))
wall_phi = sp.factor((T[1,2]/r).subs(cont))

print("Sahoo parameter mapping: L = mu_c(Sahoo)*Omega/mu = 2*K(main).")
print("DERIVED radial equation:")
print(sp.collect(expected_radial, [Fpp, Fp, Gp]))
print("PRINTED Sahoo Eq. (9):")
print(sp.collect(sahoo_radial, [Fpp, Fp, Gp]))
print("DERIVED azimuthal equation:")
print(sp.collect(expected_azimuthal, [Gpp, Gp]))
print("PRINTED Sahoo Eq. (10):")
print(sp.collect(sahoo_azimuthal, [Gpp, Gp]))
print("DERIVED full axial equation:")
print(sp.collect(expected_axial.subs(cont), [r**2, sp.diff(P,z)]))
print("PRINTED Sahoo Eq. (11), after continuity:")
print(sp.expand(sahoo_axial.subs(cont)))
print("DERIVED wall shear factors:")
print("radial:", wall_r)
print("azimuthal:", wall_phi)
print("PRINTED Sahoo slip factors: radial (1-LF), azimuthal (1-2LF).")
print("All three exact symbolic identities: PASS.")
