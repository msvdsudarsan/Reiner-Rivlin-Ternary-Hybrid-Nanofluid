# Nonlinear Stability and Bifurcation Analysis of Reiner–Rivlin Ternary Hybrid Nanofluid Flow over a Shrinking Rotating Disk with Chemical Reaction Effects

Python code, data, and figures for the manuscript `Main.pdf`
(LaTeX source: `Main.tex`, references: `References.bib`).

**Authors:** Sri Venkata Durga Sudarsan Madhyannapu (corresponding author),
Kankipati Subbarao — Department of Mathematics, Dr. RVR NRI Institute
of Technology (Deemed to be University), Andhra Pradesh, India.

## Repository contents

- `Main.pdf` / `Main.tex` / `References.bib` — the manuscript.
- `python/` — the complete computational implementation (NumPy/SciPy).
  Every numerical value and figure in the manuscript was produced by
  this code and independently re-executed by the corresponding author
  in Google Colab, with console output confirming exact agreement.
  See `python/run_all.py` to reproduce the full workflow, or run the
  individual scripts for Examples 1–4, validation, and figure
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

## Citation

If you use this code, please cite the associated manuscript
(citation details to be added upon publication).
