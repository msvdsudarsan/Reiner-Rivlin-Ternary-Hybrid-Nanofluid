"""
effective_properties.py

Ternary hybrid nanofluid effective-property correlations (viscosity,
thermal conductivity, mass diffusivity) implementing Eqs. 7-10 of the
manuscript.
"""

from dataclasses import dataclass
from parameters import Parameters


@dataclass
class EffectiveProperties:
    rho_thnf: float
    rhocp_thnf: float
    mu_thnf: float
    k_thnf: float
    D_thnf: float
    nu_thnf: float
    Pr_thnf: float
    Sc_thnf: float


def effective_properties(p: Parameters) -> EffectiveProperties:
    phi1, phi2, phi3 = p.phi1, p.phi2, p.phi3

    # density, Eq. (7)
    rho_thnf = (1 - phi3) * ((1 - phi2) * ((1 - phi1) * p.rho_f + phi1 * p.rho_s1) + phi2 * p.rho_s2) + phi3 * p.rho_s3

    # heat capacity, Eq. (8)
    rhocp_f = p.rho_f * p.cp_f
    rhocp_s1 = p.rho_s1 * p.cp_s1
    rhocp_s2 = p.rho_s2 * p.cp_s2
    rhocp_s3 = p.rho_s3 * p.cp_s3
    rhocp_thnf = (1 - phi3) * ((1 - phi2) * ((1 - phi1) * rhocp_f + phi1 * rhocp_s1) + phi2 * rhocp_s2) + phi3 * rhocp_s3

    # viscosity, Eq. (9): extended Brinkman form
    mu_thnf = p.mu_f * (1 - phi1) ** (-2.5) * (1 - phi2) ** (-2.5) * (1 - phi3) ** (-2.5)

    # thermal conductivity, Eq. (10): sequential Maxwell/Hamilton-Crosser
    k_nf = p.k_f * (p.k_s1 + 2 * p.k_f - 2 * phi1 * (p.k_f - p.k_s1)) / (p.k_s1 + 2 * p.k_f + phi1 * (p.k_f - p.k_s1))
    k_hnf = k_nf * (p.k_s2 + 2 * k_nf - 2 * phi2 * (k_nf - p.k_s2)) / (p.k_s2 + 2 * k_nf + phi2 * (k_nf - p.k_s2))
    k_thnf = k_hnf * (p.k_s3 + 2 * k_hnf - 2 * phi3 * (k_hnf - p.k_s3)) / (p.k_s3 + 2 * k_hnf + phi3 * (k_hnf - p.k_s3))

    # mass diffusivity: Maxwell/Fricke obstruction limit for impermeable
    # solid inclusions, applied sequentially (see effective_properties.m
    # for the physical justification and the caveat that this is a
    # modelling choice, not an experimentally validated ternary result).
    D_nf = p.D_f * 2 * (1 - phi1) / (2 + phi1)
    D_hnf = D_nf * 2 * (1 - phi2) / (2 + phi2)
    D_thnf = D_hnf * 2 * (1 - phi3) / (2 + phi3)

    nu_thnf = mu_thnf / rho_thnf
    Pr_thnf = (mu_thnf * rhocp_thnf / rho_thnf) / k_thnf
    Sc_thnf = nu_thnf / D_thnf

    return EffectiveProperties(rho_thnf, rhocp_thnf, mu_thnf, k_thnf, D_thnf, nu_thnf, Pr_thnf, Sc_thnf)
