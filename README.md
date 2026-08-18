# Reiner-Rivlin Ternary Hybrid Nanofluid — Computational Reproducibility Repository

Python code, data, and figures supporting the manuscript "Structural
Decoupling and Regularity-Boundary Analysis of Reiner-Rivlin Ternary
Hybrid Nanofluid Flow over a Shrinking Rotating Disk with Chemical
Reaction Effects" (S. V. D. S. Madhyannapu, K. Subbarao; manuscript
under preparation for submission).

This repository is maintained as a **computational reproducibility
archive** — it contains the code, data, and figures needed to
reproduce every numerical result in the manuscript, but not the
manuscript source itself, which remains unpublished pending journal
submission.

## Repository contents

- `python/` — the complete computational implementation (NumPy/SciPy).
  Every numerical value and figure in the manuscript was produced by
  this code and independently re-executed by the corresponding author
  in Google Colab, with console output confirming exact agreement.
  See `python/run_all.py` to reproduce the full workflow, or run the
  individual scripts for Examples 1-4, validation, and figure
  generation separately. `python/README_python_notes.md` documents
  real numerical issues found and fixed during development.
- `python/colab_ready_scripts/` — six self-contained, independently
  tested scripts, each runnable directly in Google Colab with no
  installation needed.
- `data/` — CSV/JSON numerical output underlying the manuscript's
  tables.
- `figures/` — the 6 figures used in the manuscript, generated at
  500 DPI and independently reproduced by the corresponding author.

## Reproducing the results

Open any script in `python/` (or, more easily, `python/colab_ready_scripts/`)
in Google Colab, or run locally with `numpy`, `scipy`, and
`matplotlib` installed. Each script prints its output alongside the
corresponding value from the manuscript for direct comparison.

**Package versions used for the reported results:** Python 3, NumPy
2.4.4, SciPy 1.17.1 (`scipy.integrate.solve_bvp`). All values were
re-verified against these exact versions on 18 Aug 2026; earlier
verification passes used the current Google Colab default versions
at the time of that run. `solve_bvp`'s adaptive mesh refinement is
deterministic given identical solver inputs (initial mesh, initial
guess, tolerance, `max_nodes`), but converged mesh node counts are
not guaranteed to be identical across different SciPy versions or
different initial guesses; the physical quantities reported in the
manuscript (e.g. $F'(0)$) were confirmed stable across the versions
and initial guesses tested.

## Citation

If you use this code, please cite the associated manuscript
(citation details to be added upon publication).

## Revision history

This repository reflects a full peer-review incorporation pass
(five independent reviews). Key updates: corrected six bibliography
entries with verified real authors, corrected one mis-cited DOI,
grounded the chemical species with a real literature diffusivity
value, and added new computational results -- a monotonic lambda_c
vs. K trend across 5 values of K, an extended reaction-rate parameter
sweep, and a verified physical explanation for a thermal-gradient
finding at lambda=-2 (all confirmed by direct computation, not
asserted). See the manuscript's Nomenclature table and Results
section for details.
