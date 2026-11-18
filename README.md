# Reiner-Rivlin rotating-disk flow — computational reproducibility archive

Python code, data and figures supporting the manuscript

> **The second normal-stress coefficient sets the regularity boundary of the
> Reiner-Rivlin rotating-disk similarity reduction**
>
> S. V. D. S. Madhyannapu, K. Subbarao, G. Arunagiri, P. K. Cintaginjala,
> A. Kiran Kumar, D. R. Krishna Thippisetti
>

This repository is a computational reproducibility archive: it holds the code,
data and figures needed to reproduce every numerical result in the manuscript,
not the manuscript source.

The baseline solver and validation record date from an earlier formulation of this work as a ternary hybrid nanofluid heat-transfer study; the rheological identification of the cross-viscosity, the negative-K sweep, the species spectrum, the Weissenberg-number collapse and the pressure-compatibility and application checks are collected in `python/rr_constitutive/` and `python/symbolic_derivation/`.

## Repository contents

- `python/` — the complete computational implementation (NumPy/SciPy).
  The numerical values and figures in the manuscript were produced by
  this code. The authors re-ran the scripts in Google Colab earlier; the
  runs recorded in `logs/` (November 2026, six scripts listed in
  `logs/SUMMARY.txt`) are the only ones whose output is archived here.
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
  `momentum_eigenvalue_trend.csv` (a historical positive-only eigenvalue
  diagnostic, retained to document the one-sided scan over gamma >= 0; it lists
  the smallest positive eigenvalues, not the smallest eigenvalues, and the
  the stability results use `stability_branch.csv` and the script described below),
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

**Execution status.** All repository scripts were executed without reported runtime errors, including the documented argument-dependent scripts (logs in `logs/` and `logs/full_run/`). The conditioning table (relative difference below 1e-11), the negative-K sweep (up to rounding of printed digits) and the rheometry workflow output (all 36 numerical entries, relative tolerance 1e-9) were regenerated and checked against the shipped data. A complete output-by-output comparison of all generated results against the archived data has not been performed; full repository-wide reproducibility is therefore not claimed.

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


## Rheometry workflow and independent Sahoo equation audit

`python/rheometry_workflow.py` demonstrates the mapping
`K = Psi2 * Omega / eta0`, a preliminary wall-factor screen, a full-profile
denominator check for the computed `lambda=S=0` state, and uncertainty
propagation. The current input dictionary is explicitly **synthetic**.
It is not a measured fluid dataset. The Monte Carlo torque calculation
recomputes the sample-specific kinematic viscosity `nu_i=eta_i/rho` for
each draw; only negative-K samples within the computed continuation grid
are interpolated. The regime figure labels the exact wall-singularity
locus, not a proved global existence boundary. Run from this directory:

```bash
cd python
python rheometry_workflow.py
python symbolic_derivation/sahoo_comparison.py
```

`python/symbolic_derivation/sahoo_comparison.py` independently constructs
the cylindrical rate-of-strain tensor and all three components of the
stress divergence for Sahoo (2009)'s printed constitutive law. It checks
the radial, azimuthal and axial reductions symbolically. The coefficient
mapping is `L_Sahoo = 2 K_main`, because Sahoo writes
`2 mu e + 2 mu_c e.e` while the main text writes
`mu A1 + mu_c(main) A1^2` with `A1=2e`. The audit records equation-level
differences without declaring the published numerical results incorrect;
those results have not been recomputed. See the archived stdout logs.

## 1. The cross-viscosity is the second normal-stress coefficient

Evaluating the Reiner-Rivlin law in steady simple shear gives `N1 = 0` and
`N2 = mu_c * gammadot^2`, so `mu_c = Psi_2` exactly and the governing group is
`K = Psi_2 * Omega / eta_0`. The sign of
the second normal-stress coefficient depends on the fluid class (Maklad and
Poole, J. Non-Newtonian Fluid Mech. 292:104522, 2021). It is negative for
polymeric solutions, micellar systems and dense suspensions, but positive
values are recorded for some liquid-crystalline and magnetorheological
samples. The negative-`K` branch is rheologically motivated by these fluid classes;
this does not establish that the Reiner-Rivlin constitutive model
quantitatively describes any particular fluid in those classes. The
rotating-disk literature located for this work had not examined that branch.

## 2. On that branch the boundary is wall-attained, and the critical rate is closed-form

For `K > 0` the singular set `2KF = 1` is reached at an interior maximum of
`F`, and the critical shrinking rate has no closed form. For `K < 0` it is
reached at the wall, where `F(0) = lambda` exactly, giving

    lambda_c = 1 / (2K)

`python/rr_constitutive/run_negK_sweep.py` verifies this across eleven
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

`python/rr_constitutive/run_species_spectrum.py` returns
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
`gamma_1^(M) = 0.41857` at the reference state to within one grid step (finite-domain value, see mode_evidence.py).

## Honest numerical notes

* Near-boundary momentum eigenvalues are ill-conditioned; the historical positive-only scan files were removed and the stability results rely on `data/stability_branch.csv` only.
* Comparison with the closest prior shrinking-disk study (Puspanathan et al., 2024).
  Their suction parameter is half of ours (S = 2 x their value), so their
  operating point corresponds to S = 5.6, K = 1.5, lambda = -2.4. With the
  no-radial-pressure-gradient closure used throughout, the first branch is
  reproduced (F'(0) = 2.51900), a second solution lies near their second
  branch (F'(0) about 2.083-2.086, matched only to about 1e-3), and natural
  continuation from lambda = 0 reaches a third solution (F'(0) = 4.0383) that is
  not a published branch. The value 1.7190 belongs to S = 2.8, a different
  operating point. See `python/colab_ready_scripts/09_puspanathan_direct_comparison.py`.

## Application example and pressure compatibility

`python/application_example.py` computes the idealized single rotating-disk
scenario of the manuscript (all inputs assumed, no measured fluid) and writes
its output to `data/application_example_log.txt`.
`python/symbolic_derivation/pressure_compatibility.py` checks symbolically that
the closure with no radial pressure gradient is an assumption for K different
from zero, not an exact three-dimensional reduction.

## Running the new scripts

    pip install -r requirements.txt
    cd python/rr_constitutive
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

Reproduce with `python/rr_constitutive/run_negK_sweep.py`; the collapse is a
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

`python/rr_constitutive/run_stability_check.py` redoes the analysis without
that blind spot:

* the base state is solved by Chebyshev collocation on the same grid as the
  eigenvalue problem, reproducing F'(0) = -0.442969 and -G'(0) = 1.368066 at the
  baseline exactly, and the spectrum is found by a direct generalised eigensolve;
* independently, the linearised equations are integrated forward in time and the
  growth rate measured, with no eigensolver involved.

Results at S = 1.5: the smallest reliable eigenvalue crosses zero at
lambda* = -0.4469 for K = +0.3 and lambda* = -0.6036 for K = -0.3, identical to
five decimals for (N, eta_inf) = (100, 30), (160, 50) and (200, 60). The
perturbation boundary conditions are f(0) = g(0) = 0, g(inf) = 0, f'(inf) = 0, with
the nonlocal pressure term kept. The far-field operator has its spectrum edge at
the origin, so eigenvalues very close to zero on a truncated domain can belong to
the continuum; only the isolated interior mode that crosses zero and the clearly
negative eigenvalues beyond it are used. The time integration uses the same spatial
operator, so it checks the eigenvalue extraction, not the discretisation.

Not claimed: a stability threshold for any other K, any suction other than S = 1.5,
or the Newtonian disk. The earlier ten-K map and the suction sweep were produced with
a routine that imposed the wrong far-field condition and dropped small eigenvalues;
they were withdrawn and their files removed.

## Mode identification

`mode_evidence.py` prints eigenpair residual, outer-half amplitude share and N / eta_inf convergence for the baseline mode (0.41857, a finite-domain value; erratic for eta_inf >= 60 at that baseline) and for both crossing modes. `shrinking_disk_example.py` reproduces the full-tensor shrinking-disk calculation (assumed inputs). `run_stability_check.py` locates the crossings; the
higher eigenvalues belong to a dense group whose values shift in the third decimal
with the resolution. Only real eigenvalues are examined; nothing is claimed about
oscillatory modes or about further unstable modes along the branch. The time
marching (implicit Euler, dt = 0.01) gives growth rates 0.5903 (K = +0.3,
lambda = -1.0), 0.4541 (K = -0.3, lambda = -0.8) and 1.0625 (K = -0.3, lambda = -1.0),
against 0.5903, 0.4553 and 1.0625 predicted from the eigenvalues.

## Checking the manuscript against this code

The numerical values in the manuscript were checked against the output of the
scripts in this archive before release (the archived November 2026 logs cover
six scripts only, see `logs/SUMMARY.txt`): the baseline wall gradients, the
momentum eigenvalues, both critical shrink rates, both stability crossings, the
second crossing, all three time-integration growth
rates, the species eigenvalue, all eleven rows of the negative-K table, the
collapse statistics and the symbolic N1 check. All agree to the precision
reported in the paper.

    cd python/rr_constitutive
    python run_stability_check.py       # crossings and time-marching check
    python make_stability_figure.py     # data/stability_branch.csv, figures/stability_corrected.pdf

## Three-dimensional disturbances: inviscid crossflow diagnostic

`python/rr_constitutive/run_crossflow_diagnostic.py` applies the inviscid
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

## Extended range of the closed form

`data/negK_extension.csv` records continuation checks that
extend the verified range of lambda_c = 1/(2K) to K in [-1.0, -0.05] and S in [0.5, 6].

## Conditioning near the regularity boundary

`python/rr_constitutive/run_conditioning_check.py` asks whether the end of the
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

