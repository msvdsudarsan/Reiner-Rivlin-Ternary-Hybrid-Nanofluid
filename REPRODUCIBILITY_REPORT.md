# Reproducibility report

Everything below was run from the repository root with Python 3.13, NumPy, SciPy, Matplotlib and SymPy. Raw output of the full run is in `logs/run_all.log`.

| Item | Command | Result |
|---|---|---|
| Full workflow | `python python/run_all.py` | completed without error; regenerates `data/compatible_*.csv`, the figures, the disk example and the rheometry workflow |
| Tests | `python -m pytest tests` | 4 passed (Newtonian limit, coincidence of closures at K = 0, closures differ at K > 0, minimum of 1 - 2KF at the wall for K < 0) |
| Symbolic check | `python python/symbolic_derivation/pressure_compatibility.py` | radial numerator F^2 - G^2 + HF' + K(F'^2 + 3G'^2); axial component fixes P_0 |
| Newtonian limit | tests | F'(0) = 0.5102, -G'(0) = 0.6159, -H(inf) = 0.8845 within 5e-5 |

Comparison with the manuscript: the wall derivatives, the eleven critical rates, the stability rows, the disk example (torque +0.17 %, entrainment +1.78 %, rim shear rate 905 1/s) and the rheometry intervals printed in `logs/run_all.log` are the values quoted in the manuscript and the Supplementary Material.

Limits. A recheck at tolerance 1e-10 converged for only one of eleven compatible cases near termination, so those values are quoted at 1e-8; the gap to 1/(2K) reflects the continuation step. Wall attainment and the absence of a fold are numerical evidence, not proofs. No published solution of the compatible problem exists to compare with; the zero-gradient limit reproduces the published first branch at K = 1.5, lambda = -2.4, S = 5.6 (F'(0) = 2.5190). The environment was not pinned to a lock file. The Chebyshev stability computation uses strict decaying conditions; earlier thresholds obtained with a relaxed far-field condition are withdrawn.
