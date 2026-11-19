# Reiner-Rivlin rotating-disk computational archive

Code and data for the manuscript "Regularity termination of the Reiner-Rivlin rotating-disk similarity reduction governed by the second normal-stress coefficient".

## What is computed

The steady von Karman flow of a Reiner-Rivlin fluid over a rotating, permeable, shrinking disk, with the pressure field that satisfies all three momentum components (compatible closure) and with the common zero-gradient closure for comparison. The dimensionless cross-viscosity is K = Psi_2 Omega / eta_0.

| Quantity | Script | Output |
|---|---|---|
| Wall shear, both closures | `python/primary_results.py` | `data/compatible_closure_comparison.csv`, `data/compatible_baseline_and_positive_K.csv` |
| Critical rate 1/(2K), eleven cases | `python/primary_results.py` | `data/compatible_critical_rates.csv` |
| Stability to similarity-preserving perturbations | `python/primary_results.py`, `python/make_primary_figures.py` | `data/compatible_stability.csv`, `data/compatible_growth_rate_curves.csv` |
| Figures | `python/make_primary_figures.py` | `figures/*.pdf` |
| Assumed laboratory disk | `python/application_example.py` | printed |
| Rheometry workflow (synthetic inputs) | `python/rheometry_workflow.py` | printed |
| Symbolic derivation | `python/symbolic_derivation/pressure_compatibility.py`, `second_order_check.py` | printed |

`python/compatible_model.py` holds the solver (`solve_bvp` with continuation) and the Chebyshev stability computation.

## Run

```
pip install -r requirements.txt
python python/run_all.py
python -m pytest tests
```

## Notes

`REPRODUCIBILITY_REPORT.md` lists what was run and what agreed. `REFERENCE_VERIFICATION.md` and `data/reference_audit_ledger.csv` record how the citations were checked. The code is released under the MIT licence.
