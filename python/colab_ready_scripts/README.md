# Colab-ready scripts

Twelve self-contained scripts, each independently tested and each
independently re-run by the corresponding author with matching
results (see PYTHON_OUTPUT_VALUES... output archived with this
repository). No installation needed -- open Google Colab, paste one
script per cell, run.

1. `01_example1_baseline.py` -- baseline solution (Example 1)
2. `02_newtonian_validation.py` -- classical limit check
3. `03_example2_lambda_c.py` -- critical shrinking parameter (Example 2)
4. `04_example3_eigenvalues.py` -- momentum + thermal stability eigenvalues (Example 3)
5. `05_example4_table.py` -- thermal/species response table (Example 4)
6. `06_all_figures_500dpi.py` -- generates and downloads all 6 manuscript figures at 500 DPI
7. `07_peer_review_additions.py` -- reproduces four new results added during the first peer-review incorporation round: the monotonic lambda_c vs K trend, the extended beta sweep, the Corollary 1 numerical illustration, and the H(eta) zero-crossing explanation for theta'(0)->0 at lambda=-2
8. `08_broad_search_and_species_threshold.py` -- reproduces two major results added during the second peer-review incorporation round: a systematic 3x5 (K,S) grid search (14 combinations) confirming every one terminates at the regularity boundary once adequately resolved, directly responding to reviewer requests for a broader multiplicity search; and the analytical resolution of the species-block stability question (an exact stability threshold, not a discrete eigenvalue search). Takes several minutes to run in full -- this is expected, not a bug.
9. `09_puspanathan_direct_comparison.py` -- direct parameter-matched comparison against Puspanathan et al. (2024)'s exact reported operating point (K=1.5, S=2.8, lambda=-2.4), reporting the honest (not fully resolved) discrepancy documented in the manuscript's Limitations section
10. `10_chebyshev_species_search.py` -- Chebyshev-Tau spectral collocation search for higher discrete species-block modes: validates the method against a known exact spectrum first, then applies it to the real problem, finding no resolution-stable eigenvalue -- a second, independent line of evidence supporting the analytical stability threshold in Remark on the species block
11. `11_convergence_study.py` -- genuine mesh/tolerance-refinement convergence study for the momentum block (tol = 1e-6 .. 1e-11, all cleanly converged and matching the manuscript's Table "tolsweep" exactly), including an honest floating-point over-refinement finding at tol=1e-12
12. `12_thermal_conductivity_sensitivity.py` -- sensitivity of the effective thermal conductivity, and hence -theta'(0), to the choice between the manuscript's sequential Maxwell/Hamilton-Crosser closure and a one-step alternative (volume-weighted effective particle conductivity); both are theoretical closures, since no experimentally validated ternary thermal-conductivity correlation for this combination was located in the literature
