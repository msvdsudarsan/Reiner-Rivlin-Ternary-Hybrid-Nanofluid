# Colab-ready scripts

Seventeen self-contained scripts, each independently tested and each
independently re-run by the corresponding author with matching
results (see PYTHON_OUTPUT_VALUES... output archived with this
repository). No installation needed -- open Google Colab, paste one
script per cell, run.

1. `01_example1_baseline.py` -- baseline solution (Example 1)
2. `02_newtonian_validation.py` -- classical limit check
3. `03_example2_lambda_c.py` -- critical shrinking parameter (Example 2)
4. `04_example3_eigenvalues.py` -- momentum + thermal stability eigenvalues (Example 3); illustrative shooting cross-check only for the momentum eigenvalue -- see script 14 for the authoritative value
5. `05_example4_table.py` -- thermal/species response table (Example 4)
6. `06_all_figures_500dpi.py` -- generates and downloads all manuscript figures at 500 DPI
7. `07_peer_review_additions.py` -- reproduces four new results added in the extended analysis round: the monotonic lambda_c vs K trend, the extended beta sweep, the Corollary 1 numerical illustration, and the H(eta) zero-crossing explanation for theta'(0)->0 at lambda=-2
8. `08_broad_search_and_species_threshold.py` -- reproduces two major results added in the extended analysis round: a systematic 14-point (K,S) search confirming every one terminates at the regularity boundary once adequately resolved, directly responding to reviewer requests for a broader multiplicity search; and the corrected species-block far-field threshold after a sign error was found. The discrete species eigenvalue itself was later located at gamma_1 = 5.730e4 by rescaling to the wall variable; see python/rr_constitutive/run_species_spectrum.py. Takes several minutes to run in full -- this is expected, not a bug.
9. `09_puspanathan_direct_comparison.py` -- comparison with Puspanathan et al. (2024) at K=1.5, lambda=-2.4. Their vartheta=2.8 equals S=5.6 here: Part A solves in their variables and finds the two published branches; Part B continues in lambda from 0 (S=5.6 gives a third solution, F'(0)=4.0383; S=2.8 is a different operating point, F'(0)=1.7190)
9b. `09b_puspanathan_kcontinuation.py` -- continuation in K at fixed lambda=-2.4 and S=5.6; reaches the same third solution (4.0383) as continuation in lambda; writes `puspanathan_kcontinuation_check.csv`
10. `10_chebyshev_species_search.py` -- Chebyshev-Tau spectral collocation search for higher discrete species-block modes: validates the method against a known exact spectrum first, then applies it to the real problem, finding no resolution-stable eigenvalue -- a second, independent line of evidence supporting the analytical stability threshold in Remark on the species block
11. `11_convergence_study.py` -- genuine mesh/tolerance-refinement convergence study for the momentum block (tol = 1e-6 .. 1e-11, all cleanly converged and matching the manuscript's Table "tolsweep" exactly), including an honest floating-point over-refinement finding at tol=1e-12
12. `12_thermal_conductivity_sensitivity.py` -- sensitivity of the effective thermal conductivity, and hence -theta'(0), to the choice between the manuscript's sequential Maxwell/Hamilton-Crosser closure and a one-step alternative (volume-weighted effective particle conductivity); both are theoretical closures, since no experimentally validated ternary thermal-conductivity correlation for this combination was located in the literature
13. `13_H_eta_comparison.py` -- generates the H(eta) comparison across the shrinking cases used in Example 4, giving direct visual confirmation of the near-wall sign-reversal mechanism behind the theta'(0)->0 finding at lambda=-2
14. `14_momentum_eigenvalue_chebyshev.py` -- finite-domain Chebyshev collocation at the Example 1 baseline, with the base state solved by K-continuation: lists the lowest positive real eigenvalues (0.41857) for several N and domain lengths. Positive eigenvalues only; the stability crossings are handled by script 15.
15. `15_stability_crossings.py` -- the zero crossings of the smallest real eigenvalue at K = +0.3 and K = -0.3 (S = 1.5), the eigenpair-residual and localisation evidence, and the stability figure. Numerical reproducibility only, not a half-line spectral certificate.
16. `16_denominator_clipping_sensitivity.py` -- tests whether the reported approach to the regularity boundary (lambda_c) depends on the artificial floor placed on the singular denominator Delta=1-2KF near F=1/(2K); re-runs the K=0.3, S=1.5 continuation to lambda=-5.02 at three floor tolerances spanning four orders of magnitude (1e-8, 1e-10, 1e-12) and confirms an identical result at each.
