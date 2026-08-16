# Python implementation notes

These notes document what actually happened while running the full
computational workflow in Python/SciPy, including two real bugs found
and fixed along the way (also fixed in the corresponding MATLAB files).
Full details and numerical values are in `main.tex`, but this file
collects the debugging narrative in one place.

## Bug 1: monolithic 9-state base-state solve was not reliably converged

A first attempt solved F, G, H, theta, phi together as one 9-state
`scipy.integrate.solve_bvp` problem. Symptom: `F'(0)` shifted from
`-0.4430` at `tol=1e-9` to `-0.4454` at `tol=1e-10`, and tightening
further exhausted memory (`max_nodes` exceeded).

Root cause: the Schmidt number for these parameters is
`Sc_thnf ~ 883` (an ordinary base-fluid mass diffusivity, `D_f=1e-9
m^2/s`, against `nu_thnf ~ 8e-7 m^2/s`), which makes the concentration
boundary layer roughly 30x thinner than the momentum layer. Because
`solve_bvp` refines its mesh globally, this dragged mesh refinement
across all nine states even though only `phi` needed it.

Fix: exploit Proposition 1 of `main.tex` directly. The momentum
sub-system (F,G,H) does not depend on theta or phi at all, so it is
solved on its own (5 states); theta and phi are then each solved as
their own separate 2-state problems using the already-known `H(eta)`.
See `base_state.py`. Converged to 6 decimal places across
`tol = 1e-9` to `1e-11` after this fix.

## Bug 2 (methodological, not a code bug per se): naive fresh solves at
## different eta_inf gave spurious non-monotonic "mesh independence" results

A first mesh/domain-independence check solved each `eta_inf` value
independently from a fresh, generic exponential-decay initial guess.
Result: `F'(0)` at `eta_inf=15` and `eta_inf=30` agreed
(`-0.442969`), but `eta_inf=20` and `25` gave different, seemingly
non-converged values (`-0.444237`, `-0.444476`) -- non-monotonic and
inconsistent with genuine domain-independence.

Fix: warm-start each larger domain from the previous domain's
converged solution (continuation in `eta_inf`) instead of a fresh
guess each time. After this fix, all `eta_inf` in {15,20,25,30,35}
agree to 6 decimal places. Interpretation: the naive fresh guess was
landing near a nearby-but-distinct fixed point of Newton's method for
this moderately nonlinear system at some domain sizes -- a genuine,
reportable numerical-methodology finding (documented in `main.tex` as
part of Section "Numerical Methodology"), not evidence of multiple
physical solutions at those `eta_inf` values.

## Bug 3: eigenvalue solve initial-guess sign error

The perturbation eigenvalue problems (momentum, thermal, species
blocks) require a normalization boundary condition, e.g.
`Ff'(0) = 1`. An early implementation used an initial guess with
`Ff'(0) = -1` (the derivative of a plain `exp(-eta)` decay guess),
which is the WRONG SIGN relative to the boundary condition. This
caused every eigenvalue `solve_bvp` call to fail with "maximum number
of mesh nodes exceeded" regardless of tolerance or initial eigenvalue
guess, and was at first misdiagnosed as numerical stiffness near a
fold.

Fix: use `eta*exp(-eta)` as the guess shape, which satisfies
`[field](0)=0` and `[field]'(0)=1` exactly. This immediately fixed the
thermal-block eigenvalue solve. The momentum and species blocks
additionally needed a switch from `solve_bvp`'s unknown-parameter
mechanism to a direct two-point shooting method (`stability_shooting.py`),
which proved more robust for this specific linear homogeneous
eigenvalue problem. **This exact bug also existed in the delivered
MATLAB `stability_eigenvalues.m`** and has been fixed there too.

## Correction to Proposition 2 (block-triangular -> block-diagonal, 3 blocks)

Re-reading the linearized equations during implementation showed that
Theta and Phi perturbations do not couple to EACH OTHER at all (each
is forced only by the momentum perturbation Hf). The original
`main.tex` and `matlab/stability_eigenvalues.m` treated them as one
combined "scalar-transport" block sharing a single eigenvalue -- this
is wrong. The corrected structure is three independent diagonal
blocks: momentum {F,G,H}, thermal {Theta}, species {Phi}, each with
its own eigenvalue spectrum. Both `main.tex` (Remark on finer
block-diagonal structure) and `matlab/stability_eigenvalues.m` have
been updated to reflect this.

## What was NOT successfully resolved in this run

- The species-block eigenvalue could not be reliably converged (large
  Sc_thnf makes the perturbation layer very thin) by either
  `solve_bvp` or the shooting method within the available time
  budget. Open item -- see `main.tex` "Results and Discussion" and
  "Limitations".
- A momentum-eigenvalue trend approaching the regularity boundary
  (lambda -> -5.026) was attempted but not completed reliably.
- A genuine pseudo-arclength crossing onto a second branch was not
  achieved (the reduced-plane scheme converges to the same boundary
  point as bisection but does not cross it). Given the regularity-
  boundary finding (see `main.tex` Example 2), it is not clear there
  IS a second branch to cross onto for the parameter combinations
  tested.

## Key finding worth restating

For K=0.3 and S in {0.5, 1.0, 1.5, 2.0, 3.0}, no genuine saddle-node
(dual-solution) bifurcation was found. Apparent "folds" (natural
continuation failing to converge) coincide, in every case checked,
with `max(F(eta))` approaching the analytical singular value
`1/(2K)` from `main.tex` Eq. (19) -- confirmed by reducing `K` to 0.05
(pushing the singularity out to `F=10`), at which point continuation
runs cleanly with no failure at all from lambda=0 to lambda=-15. This
is reported as an honest negative/clarifying result, not reframed to
look like a successful bifurcation study.

## Bug 4 (found via MATLAB cross-check, most consequential so far): eta_inf=15 was insufficient for K=0.3 cold-start solves

When the user ran the corrected MATLAB code, `main_rotating_disk.m`
reported a fully-converged (residual < requested tolerance, no
warning) momentum baseline of `F'(0) = -0.445403` at K=0.3, lambda=0.5,
S=0.5, eta_inf=15 -- disagreeing with the Python reference value of
`-0.442969` used throughout `main.tex`, despite both being genuinely
converged solutions (not a tolerance/mesh-cap issue).

Reproduced directly in Python: at eta_inf=15, a coarser initial mesh
(100 points) also converges to `-0.445403`, matching MATLAB, while a
finer one (200+ points) converges to `-0.442969`, BOTH to residuals
~1e-10. This first looked like genuine solution multiplicity of the
Reiner-Rivlin momentum BVP, and a "robust" K-continuation solver
(`solve_momentum_robust`, starting from the unique K=0 classical
solution and stepping K up to the target) was added, which reproducibly
gave `-0.445403` regardless of starting mesh density -- suggesting
`-0.445403` was the "true", continuously-connected branch and
`-0.442969` was the spurious one.

**This conclusion was itself wrong and had to be corrected again.**
Testing the SAME K-continuation procedure at LARGER eta_inf revealed
that `Fp0` drifts with eta_inf up to eta_inf~35, at which point it
snaps to, and then stays stable at, `-0.442969` for eta_inf in
{35, 40, 50} (confirmed to 8 decimal places, residuals ~1e-10). So:
**eta_inf=15 is simply too short a computational domain for K=0.3 at
these parameters** -- both `-0.445403` (found by MATLAB and by a
coarse-mesh or short-domain Python solve) and the various intermediate
values obtained by K-continuation at eta_inf between 15 and 30 are
converged solutions of the TRUNCATED (eta_inf too small) problem, not
good approximations to the true semi-infinite-domain boundary-value
problem that the model is meant to represent. `-0.442969` -- the
original value used throughout `main.tex` -- is confirmed to be the
correct one.

Separately checked and NOT affected: the Example 2 continuation search
(K=0.3, S=1.5, tracing lambda from 0 downward with proper warm-starting
at every step) gives an IDENTICAL stall point (lambda=-5.0000,
F'(0)=1.7025, max(F)=1.6513) at eta_inf=15 and eta_inf=40 -- so the
regularity-boundary finding of Example 2 does not need revision. The
sensitivity found here appears specific to COLD-START solves at a
single lambda value, not to properly warm-started continuation.

**Fix applied:** `eta_inf` default increased from 15 to 40 in both
`parameters.py` and `matlab/parameters.m`; `Nmesh0` scaled up to match.
`solve_momentum_robust` (K-continuation from K=0) is retained as the
default solving method for extra insurance against any further
multiplicity, now combined with a domain size confirmed sufficient.

**Lesson:** a "converged" (low-residual, no-warning) solution from
`bvp4c`/`solve_bvp` is not by itself evidence that the truncated-domain
approximation is adequate -- domain independence must be checked by
explicitly pushing eta_inf well beyond the initially chosen value, not
assumed from a single converged run, and not fully validated by
continuation-in-eta_inf alone if the continuation is only run over a
narrow eta_inf range (the original mesh-independence check in this
project only tested eta_inf up to 35 starting FROM a wrong eta_inf=15
solution, which is why it did not catch this earlier).

## Bug 5 (found via a THIRD round of MATLAB cross-checking): even K=0 was not reliably obtained by a cold solve at large eta_inf

After Bug 4 was "fixed" (eta_inf raised to 40, K-continuation from K=0
added as a safeguard), the user ran the corrected MATLAB code and got a
NEW failure: the K=0 starting solve for the K-continuation itself
failed to converge, with `bvp4c` reporting residuals around 0.03 even
after growing to 20,000+ mesh points -- not just slow, genuinely stuck.

Root cause, found by reproducing directly in Python: **the K=0
(classical Newtonian) problem at lambda=0.5, S=0.5 is ALSO not reliably
obtained by a direct cold solve at a large eta_inf.** Solving fresh at
each of eta_inf in {10,15,20,25,30,35,40,50,60}, all converging to low
residual (~1e-10, no warning), gave DIFFERENT answers: -0.287, -0.321,
-0.313, -0.309, -0.306, -0.287, -0.287, -0.287, then FAILED OUTRIGHT at
eta_inf=60. So the premise of Bug 4's fix -- "K=0 is the guaranteed-
unique classical case, safe to cold-start" -- was itself wrong for this
lambda=0.5 (non-zero stretching), S=0.5 (suction) combination. Only
pure rotation (lambda=0) is textbook-unique; stretching combined with
suction reintroduces exactly the kind of solver-path sensitivity found
for the K=0.3 case in Bug 4.

**The fix that actually works**, confirmed stable from eta_inf=8 to 80
(agreeing to 6 decimal places): a TWO-STAGE continuation.
1. Solve at a SMALL, easily-converged eta_inf (8) first.
2. Incrementally extend eta_inf (warm-starting each step from the
   previous, constant-tail-extended solution) up to the target.
3. Only once the domain-correct K=0 solution is obtained, do the
   K-continuation from stage 2 of Bug 4's fix (step K from 0 to the
   target, warm-starting each step) at the FIXED final eta_inf.

Jumping straight to a large eta_inf at K=0 (what Bug 4's fix did) is a
poorly-conditioned starting point regardless of solver -- confirmed by
reproducing MATLAB's exact failure mode in scipy (large residual,
requires far more mesh points than the well-posed small-domain problem
needs) and then confirming the two-stage approach resolves it cleanly
in both languages.

**Also fixed**: `pseudo_arclength_continuation.m`'s (and
`continuation.py`'s) OWN starting point (lambda=0, at the target K)
had exactly the same cold-start-at-large-domain flaw and has been
updated to use the same robust two-stage solver.

**Final confirmed value** (unchanged from before -- this was a
methodology fix, not a numerical correction): K=0.3, lambda=0.5, S=0.5,
eta_inf=40 gives F'(0)=-0.442969, matching every value reported
throughout main.tex.

**Lesson, restated more strongly than after Bug 4**: for boundary-layer
similarity solutions with stretching/suction, do not trust ANY
cold-start solve at a large truncated domain, even at parameter values
that look "obviously unique" (like the Newtonian limit) -- always
build up to the target domain size and parameter values incrementally
from a small, well-behaved starting point.
