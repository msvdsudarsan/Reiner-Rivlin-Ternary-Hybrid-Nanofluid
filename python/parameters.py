"""
parameters.py

Python mirror of matlab/parameters.m. Values are kept identical so that
results are directly comparable to the MATLAB R2026a run the user will
perform independently. This is a PYTHON (numpy/scipy) implementation,
NOT MATLAB R2026a -- every result derived from this module must be
labelled as such in the manuscript.
"""

from dataclasses import dataclass, field


@dataclass
class Parameters:
    # Reiner-Rivlin cross-viscosity parameter
    K: float = 0.3

    # disk kinematics
    lam: float = 0.5      # F(0) = lambda ; corresponds to MATLAB p.lambda
    S: float = 0.5        # H(0) = -S

    # chemical reaction
    beta: float = 0.2

    # nanoparticle volume fractions (Cu, Al2O3, TiO2)
    phi1: float = 0.01
    phi2: float = 0.01
    phi3: float = 0.01

    # base fluid (water) and nanoparticle thermophysical properties (25 C)
    rho_f: float = 997.1
    cp_f: float = 4179.0
    k_f: float = 0.613

    rho_s1: float = 8933.0   # Cu
    cp_s1: float = 385.0
    k_s1: float = 400.0

    rho_s2: float = 3970.0   # Al2O3
    cp_s2: float = 765.0
    k_s2: float = 40.0

    rho_s3: float = 4250.0   # TiO2
    cp_s3: float = 686.2
    k_s3: float = 8.9538

    mu_f: float = 8.9e-4      # Pa.s, water at 25 C
    D_f: float = 1.0e-9       # m^2/s, representative placeholder (species-specific)

    # numerical / domain settings
    eta_inf: float = 40.0     # see base_state.py module docstring: eta_inf=15
                              # was found insufficient for K=0.3 cold-start
                              # solves (converges to a different, wrong-
                              # domain root); 35-50 confirmed stable/correct
    RelTol: float = 1e-8
    AbsTol: float = 1e-10
    Nmesh0: int = 800         # scaled up to match the larger eta_inf above


def default_parameters() -> Parameters:
    return Parameters()
