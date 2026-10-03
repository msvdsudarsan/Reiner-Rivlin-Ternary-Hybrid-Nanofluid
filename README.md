# Reiner-Rivlin rotating-disk flow — computational reproducibility archive

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22009268.svg)](https://doi.org/10.5281/zenodo.22009268)

Python code, data and figures supporting the manuscript (release 2.2.0)

> **The second normal-stress coefficient sets the regularity boundary of the
> Reiner-Rivlin rotating-disk similarity reduction**
>
> S. V. D. S. Madhyannapu, K. Subbarao, G. Arunagiri, P. K. Cintaginjala,
> A. Kiran Kumar, D. R. Krishna Thippisetti
>
> Submitted to the Journal of Non-Newtonian Fluid Mechanics.

This repository is a computational reproducibility archive: it holds the code,
data and figures needed to reproduce every numerical result in the manuscript,
not the manuscript source. The release corresponding to the submission is
archived at Zenodo under the DOI above.

**A note on the history below.** The baseline solver, validation record and
development notes in the sections that follow date from an earlier formulation
of this work as a ternary hybrid nanofluid heat-transfer study, and they are kept
unchanged, so they still refer to that earlier title. What is new to the present
manuscript — the rheological identification of the cross-viscosity, the
negative-K sweep, the species spectrum and the Weissenberg-number collapse — is
collected under "JNNFM version" further down.

## Repository contents

- `python/` — the complete computational implementation (NumPy/SciPy).
  Every numerical value and figure in the manuscript was produced by
  this code and independently re-executed by the corresponding author
  in Google Colab, with console output confirming exact agreement.
  See `python/run_all.py` to reproduce the core numerical workflow
  (baseline, validation, an illustrative continuation sweep, stability
  eigenvalues, and the Example 4 sweep); the extended analyses --
  the systematic 14-point (K,S) search and its bisection refinement
  of lambda_c, the Puspanathan et al. comparison, the corrected
  Chebyshev momentum-eigenvalue solve and its trend along the branch,
  the property-closure sensitivity checks, the H(eta) mechanism
  comparison, and the reference-verification pass -- are each
  reproduced by their own dedicated script in
  `python/colab_ready_scripts/` (see that directory's README for the
  full list). `python/README_python_notes.md` documents real numerical
  issues found and fixed during development.
- `python/colab_ready_scripts/` — seventeen self-contained, independently
  tested scripts, each runnable directly in Google Colab with no
  installation needed.
- `data/` — CSV/JSON numerical output underlying the manuscript's
  tables and figures, including `broad_search_KS.csv`
  (the 14-point (K,S) regularity-boundary search, Table "broad-search"),
  `momentum_eigenvalue_trend.csv` (the corrected gamma_1^(M)(lambda)
  trend, Figure "eigen-trend", flagged by resolution status),
  `eigenvalue_baseline_convergence.csv` (the Chebyshev convergence
  table underlying gamma_1^(M)=0.41857), `species_chebyshev_spurious_modes.csv`
  (the resolution-dependent species-block spurious eigenvalues),
  `closure_sensitivity.csv` (viscosity/thermal-conductivity/cross-viscosity
  closure sensitivity checks), and `H_eta_comparison.csv` (the H(eta)
  profiles underlying the theta'(0)->0 mechanism figure, independently
  confirming the reported zero-crossing near eta=0.5 at lambda=-2).
- `figures/` — the 12 figures used in the manuscript, as vector PDF with
  PNG previews, each regenerated from the scripts and data in this archive.

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

## Development history

The full record of how the analysis evolved, including each numerical error that was found and corrected, is in `DEVELOPMENT_HISTORY.md`. It is not needed to reproduce the current results.

# Results reproduced by this archive

The results below are the ones the manuscript rests on. Their code and data
are in `python/jnnfm_constitutive/` and the correspondingly named files under
`data/` and `figures/`; the rest of the archive is the baseline solver and its
validation record.

## 1. The cross-viscosity is the second normal-stress coefficient

Evaluating the Reiner-Rivlin law in steady simple shear gives `N1 = 0` and
`N2 = mu_c * gammadot^2`, so `mu_c = Psi_2` exactly and the governing group is
`K = Psi_2 * Omega / eta_0`. Because measured second normal-stress
coefficients are consistently negative (Maklad and Poole, J. Non-Newtonian
Fluid Mech. 292:104522, 2021), the surveyed systems have negative `Psi_2`, so
`K < 0` is the indicated branch for fluids to which this closure applies. The
rotating-disk literature located for this work had not examined that branch.

## 2. On that branch the boundary is wall-attained, and the critical rate is closed-form

For `K > 0` the singular set `2KF = 1` is reached at an interior maximum of
`F`, and the critical shrinking rate has no closed form. For `K < 0` it is
reached at the wall, where `F(0) = lambda` exactly, giving

    lambda_c = 1 / (2K)

`python/jnnfm_constitutive/run_negK_sweep.py` verifies this across eleven
`(K,S)` combinations spanning `K` from -0.05 to -0.40 and `S` from 0.5 to 3.0
(`data/negK_sweep.csv`). Agreement runs from 0.02% to 1.54%, and `argmin F`
sits at `eta = 0` in every case, confirming the mechanism rather than only the
number. No second solution branch was found anywhere on the sweep.

## 3. The species spectrum — found, and the earlier null result explained

Earlier revisions reported the species block as having no discrete eigenvalue,
after searches over `gamma` in `[-50, 400]`. That conclusion was an artefact of
the search range. The species layer has thickness
`O(1/(Sc*|H(0)|)) ~ 7.6e-4`, so the natural wall variable is `zeta = Sc*eta`,
in which the perturbation equation becomes

    Phi_zz - H(zeta/Sc) Phi_z + Lambda Phi = 0,
    Lambda = (gamma - Sc*beta) / Sc^2

Modes that are `O(1)` in `zeta` therefore correspond to `gamma ~ Sc^2 ~ 7.8e5`
— about three orders of magnitude above the interval searched. The spectrum was
never absent; the search was looking in the wrong place.

`python/jnnfm_constitutive/run_species_spectrum.py` returns
`gamma_1 = 5.730e4` by Chebyshev collocation, stable to 0.06% across
`N = 200-280` and `zeta_inf = 80-150`, and `gamma_1 = 57301.9` by independent
shooting bisected on the far-field residual — agreement better than 0.01%.
The eigenvalue is positive, so the species block is linearly stable, and all
three spectra of the block-triangular structure are now computed rather than
two computed and one open.

## Validation before any new result was computed

`rr_solver.py` was first checked against the published baseline: at
`K = +0.3`, `S = 1.5` it returns `lambda_c = -5.0262` with `max F -> 1.6665`
against the analytical bound `1/(2K) = 1.6667`, matching the value reported in
the original study. The Chebyshev eigenvalue routine reproduces the published
`gamma_1^(M) = 0.41857` at the reference state to within one grid step.

## Honest numerical notes

* Near-boundary momentum eigenvalues are ill-conditioned.
  `data/momentum_eigenvalue_negK.csv` reports values only where three
  resolutions agree, and marks the final point (`lambda = -1.62`, ~9% spread)
  as marginally resolved with its value left blank rather than recording a
  number the data do not support.
* The parameter-matched discrepancy against the closest prior shrinking-disk
  study remains open. Twelve convention combinations were tested (group
  halved, unchanged or doubled, crossed with both sign choices for the suction
  parameter and for lambda). Every reversed-sign case and every doubled-group
  case fails to converge; only the halved group converges, returning 2.3344,
  which lies between the two published branch values without reproducing
  either. This narrows the explanation without closing it.

## Running the new scripts

    pip install -r requirements.txt
    cd python/jnnfm_constitutive
    python run_negK_sweep.py          # ~5 min
    python run_species_spectrum.py    # ~2 min


---


## Weissenberg-number collapse and measurement uncertainty

`figures/regularity_collapse.pdf` shows the eleven negative-K cases against
the prediction `2*Wi2*|lambda_c| = 1`, where `Wi2 = |Psi_2|*Omega/eta_0 = |K|`
is a Weissenberg number built on the second normal-stress coefficient rather
than on a relaxation time. The computed product has mean 0.9963 with a maximum
departure of 1.54%, the largest residual belonging to the weakest cross-viscosity
(K = -0.05), where the boundary lies furthest out. The same panel set shows how
relative uncertainty in Psi_2 and eta_0 propagates into the predicted boundary:
because lambda_c depends on them only through their ratio, the propagation is
direct and unamplified.

Reproduce with `python/jnnfm_constitutive/run_negK_sweep.py`; the collapse is a
post-processing step on `data/negK_sweep.csv`.

## Is the regularity boundary an artefact of N1 = 0?

`python/symbolic_derivation/second_order_check.py` answers this symbolically.
It builds the radial momentum equation under von Karman similarity for the
Reiner-Rivlin law and for the second-order fluid, which carries both normal-stress
differences, and compares the coefficients of the two highest derivatives of F.

    Reiner-Rivlin   coeff F''' = 0
                    coeff F''  = mu Omega^2 r (1 - 2 K F) / nu

    second-order    coeff F''' = alpha_1 H Omega^3 r / nu    (nonzero -> third order)
                    coeff F''  = identical; the alpha_1 terms cancel exactly
                                 against the continuity identity H' = -2F

The singular factor is therefore produced by the A1^2 term alone - by the second
normal-stress response - and is untouched by the term generating N1. This is
invariance at the level of the coefficient only: the alpha_1 term raises the order
of the problem, the complete higher-order boundary-value problem is distinct, and
its termination point is not computed here.

The script asserts the cancellation, so it fails loudly if it ever stops holding.

## Correction: the shrinking branch loses linear stability (October 2026)

An earlier version of this archive, and of the manuscript, reported the smallest
momentum eigenvalue as positive along the shrinking branch, with a non-monotonic
dip. **That conclusion was wrong.** The eigenvalue had been located by scanning
the smallest singular value of the collocation pencil over gamma in [0, 1.2]. A
scan restricted to non-negative gamma cannot see a negative eigenvalue. The
values it returned are correct as the smallest *positive* eigenvalues, but they
are not the smallest eigenvalues.

`python/jnnfm_constitutive/run_stability_check.py` redoes the analysis without
that blind spot:

* the base state is solved by Chebyshev collocation on the same grid as the
  eigenvalue problem, reproducing F'(0) = -0.442969 and -G'(0) = 1.368066 at the
  baseline exactly, and the spectrum is found by a direct generalised eigensolve;
* independently, the linearised equations are integrated forward in time and the
  growth rate measured, with no eigensolver involved.

Results at S = 1.5: the smallest eigenvalue crosses zero at lambda* = -0.4469
for K = +0.3 and lambda* = -0.6036 for K = -0.3, both identical to four decimals
at two resolutions. At K = +0.3, lambda = -1.0 the eigensolve gives -0.5886 and
time integration a growth rate of +0.5974, the two agreeing exactly once the
damping of the implicit time step is accounted for. The negative eigenvalue is
unchanged across N = 80-150, survives imposing the far-field condition on the
perturbation itself, has a properly decaying interior eigenfunction, and is
absent at the baseline and in the Newtonian limit.

The branch therefore loses linear stability well before the regularity boundary.
The constitutive results - the identification mu_c = Psi_2, the closed form
lambda_c = 1/(2K), the Weissenberg collapse and the robustness to N1 - concern the
existence of the steady solution and are unaffected.

## Mode identification and two-branch validation

`python/jnnfm_constitutive/run_mode_validation.py` answers the two questions an
eigenvalue result has to answer before it can be trusted.

**Which eigenvalue crosses?** A routine returning the most negative eigenvalue
could in principle jump between modes, so the script prints the eigenvalues
nearest zero on either side of each crossing. In both cases a single simple real
eigenvalue passes through zero, and the next is separated by a clear gap: 0.26 at
K = +0.3 and 0.32 at K = -0.3. Following the branch further, K = +0.3 acquires a
second unstable real mode at lambda = -2.289, while K = -0.3 carries exactly one
unstable mode up to the regularity boundary. No complex eigenvalue with negative
real part appears on either branch, so the instability is non-oscillatory.

**Is it real on both branches?** The linearised equations are integrated forward
in time with no eigensolver involved. On the negative branch the measured growth
rates are 0.4595 at lambda = -0.8 and 1.0858 at lambda = -1.0, equal to the
implicit-Euler predictions from the eigensolve. Together with the K = +0.3 check
(0.5974), both crossings are confirmed by an independent method.

## Checking the manuscript against this code

Every numerical value in the manuscript was checked against the output of the
scripts in this archive before release: the baseline wall gradients, the
momentum eigenvalues, both critical shrink rates, both stability crossings, the
second crossing, the local spectrum table, all three time-integration growth
rates, the species eigenvalue, all eleven rows of the negative-K table, the
collapse statistics and the symbolic N1 check. All agree to the precision
reported in the paper.

    cd python/jnnfm_constitutive
    python run_stability_check.py       # crossings, -0.5886, +0.5974
    python run_mode_validation.py       # local spectrum, K=-0.3 time integration

## Three-dimensional disturbances: inviscid crossflow diagnostic

`python/jnnfm_constitutive/run_crossflow_diagnostic.py` applies the inviscid
stationary-crossflow condition of Gregory, Stuart & Walker (1955) to the mean
profiles. It is validated first on the Newtonian von Karman disk, where it returns
eta_c = 1.458 and a wave angle of 13.22 deg, against the classical 1.46 and 13.2 deg.

On the negative branch (K = -0.3) a critical point exists at every shrink rate
tested, including at and beyond lambda*; on the positive branch (K = +0.3) none
exists near lambda* or at lambda = -3 and -4. The critical points are unchanged to
three decimals for eta_inf = 15-30.

This is a necessary condition for a neutral inviscid stationary mode. It does not
give the viscous critical Reynolds number, and it says nothing about travelling or
streamline-curvature modes; the full three-dimensional threshold remains open.

## Stability threshold against regularity boundary

`python/jnnfm_constitutive/run_threshold_map.py` locates the stability crossing
lambda* for ten values of K between 0 and -1 at S = 1.5 and compares it with the
regularity boundary lambda_c = 1/(2K) (`data/lamstar_map.csv`,
`figures/threshold_map.pdf`). The Newtonian disk already loses stability, at
lambda* = -0.5332, so the instability is hydrodynamic in origin. lambda* stays
between -0.40 and -0.61, while lambda_c runs from -10 to -0.5; the ratio
R* = |lambda*|/|lambda_c| rises from 0.055 to 0.808 but stays below one, so the
branch becomes unstable before the boundary throughout. Crossings at K = -0.05
and K = -1.0 agree to four decimals at N = 130, eta_inf = 36.

## Suction and the extended range of the closed form

`suction_sweep()` in `run_threshold_map.py` gives the stability crossing at K = -1.0,
where the gap to the boundary is narrowest, for S from 0.5 to 6
(`data/suction_sweep_K-1.csv`). Suction delays the instability, but R* levels off
near 0.91 and the instability still precedes lambda_c = -0.5; the S = 6 value moves
by 0.3% at N = 140. `data/negK_extension.csv` records continuation checks that
extend the verified range of lambda_c = 1/(2K) to K in [-1.0, -0.05] and S in [0.5, 6].

## Conditioning near the regularity boundary

`python/jnnfm_constitutive/run_conditioning_check.py` asks whether the end of the
negative branch is a property of the solution or of the numerics. It reports three
things and writes them to `data/conditioning_check.csv` (K = -0.3, S = 1.5).

1. **Bordered Newton Jacobian.** Its condition number stays between about 1.9e7 and
   1.7e7 while 1 - 2K*lambda falls from 0.52 to 0.004, so no progressive deterioration
   of this matrix is observed. This test is limited: the wall equation is replaced by
   the boundary condition F(0) = lambda, so the singular coefficient does not enter the
   matrix at the wall, and the result does not by itself certify the singular point.
2. **Coefficient-level diagnostic.** The matrix multiplying (F'', G'') in the closed
   momentum system is (1 - 2KF) times the 2x2 identity, so its smallest singular value
   is min |1 - 2KF|. It falls from 0.52 to 0.004 and is attained at the wall (eta = 0)
   in every row.
3. **Wall curvature.** F''(0) grows from -4.1 (lambda = -0.8) to -40.4 (lambda = -1.6),
   resolved to 0.2% across N = 160, 200, 240. At lambda = -1.66 it lies between about
   -74 and -77 and is no longer converged in N (spread 3.9%), so that value is an order
   of magnitude only. Its growth slows near the boundary, so no divergence law is claimed.

Taken together with the denominator-floor test and the continuation results this is
consistent with termination being a property of the solution; none of these checks is a proof.

