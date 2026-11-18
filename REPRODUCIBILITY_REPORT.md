# Reproducibility report

Raw logs are in `logs/script_runs/` (this pass) and `logs/full_run/` (earlier complete run; kept for reference, and superseded by this pass wherever the same script appears). "Executed, exit 0" means the script ran to completion; where a comparison with the manuscript was made, it is stated.

| Script | Status | Comparison with the manuscript |
|---|---|---|
| `python/run_all.py` | executed, exit 0 | baseline values agree (F'(0) = -0.442969, -G'(0) = 1.368066, -theta'(0) = 3.631009, -phi'(0) = 443.932186) |
| `colab_ready_scripts/01` | executed, exit 0 | baseline values agree |
| `02` | executed, exit 0 (earlier pass) | Newtonian limit check; output read, not tabulated |
| `03` | executed, exit 0 | reaches the regularity boundary near -5.03 with a smaller mesh budget; the 4th-5th decimal of lambda_c is not reproduced by this script (it says so) |
| `04` | executed, exit 0 | thermal eigenvalue 10.9623 agrees; the shooting search for the momentum eigenvalue returns 0.7095, not 0.41857, so it does not confirm the collocation value (documented in the script and in the Supplementary Material) |
| `05` | executed, exit 0 | table printed; row-by-row comparison with the main-text table not repeated in this pass |
| `06` | executed, exit 0 | five figures regenerated; the eigenvalue figure is now made by `rr_constitutive/make_eigenvalue_spectrum_figure.py` |
| `07` | executed, exit 0 | H(eta) zero crossing 0.5040 against 0.50; other outputs not tabulated |
| `08` | executed, exit 0 | (K, S) search and species threshold printed; row-by-row comparison not repeated |
| `09`, `09b` | executed, exit 0 | operating-point statements as printed in the script |
| `10` | executed, exit 0 | unscaled species collocation shows no resolution-stable eigenvalue, as stated |
| `11` | executed, exit 0 (earlier pass) | convergence study printed |
| `12` | executed, exit 0 | k_thnf +0.53 %, Pr_thnf -0.52 %, -theta'(0) -0.43 % |
| `13` | executed, exit 0 | H(eta) comparison figure regenerated |
| `14` | executed, exit 0 | 0.41857 for N = 70-120 at eta_max = 25-35 (finite-domain value) |
| `15_stability_crossings` | executed, exit 0 | crossings, evidence table and figure |
| `16` | executed, exit 0 (needs a tolerance argument; default 1e-8) | the three tolerances give the same maximum F = 1.663002 |
| `rr_constitutive/run_stability_check.py`, `mode_evidence.py`, `make_stability_figure.py`, `make_eigenvalue_spectrum_figure.py`, `run_species_spectrum.py` | executed, exit 0 | stability values and the species eigenvalue agree with the Supplementary Material |
| `python/shrinking_disk_example.py` | executed, exit 0 | agrees with Supplementary Section S18 |
| `symbolic_derivation/*.py` (14 scripts) | executed, exit 0 | exit status only; their printed checks were not re-read line by line in this pass |
| `rr_constitutive/run_conditioning_check.py`, `run_crossflow_diagnostic.py` | executed, exit 0 | output not compared line by line with the supplement tables |
| `rr_constitutive/run_negK_sweep.py` | PASS (all 11 (K,S) cases completed, run as two parallel parts, `python run_negK_sweep.py 0,3,4,10 part_a` and `1,2,5,6,7,8,9 part_b`, then merged; about 8 min each) | `logs/script_runs/run_negK_sweep.log`; max deviation of lambda_c from 1/(2K) = 1.54% |
| LAMMPS | BLOCKED: no executable, inputs, logs or trajectories exist | no molecular-dynamics claim is made |

This is not a certificate of full reproducibility: several scripts were compared with the manuscript only for selected values, and the symbolic scripts were checked by exit status.


## LAMMPS (exploratory)
`lammps_exploratory/`: generic Kremer-Grest melt, four shear rates, one seed; run with `run_all.sh`, analysed with `analyze.py` (results in `kg_results.json`). Not a validation of the Reiner-Rivlin closure; see its README.
