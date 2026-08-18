# Python implementation notes

These notes document real numerical issues found and fixed during
development, kept here for transparency and reproducibility. Full
details and numerical values are in `Main.tex`; this file collects
the debugging narrative in one place.

## Finding 1: a monolithic base-state solve was not reliably converged

A first attempt solved F, G, H, theta, phi together as one 9-state
`scipy.integrate.solve_bvp` problem. Symptom: `F'(0)` shifted from
`-0.4430` at `tol=1e-9` to `-0.4454` at `tol=1e-10`, and tightening
further exhausted memory (`max_nodes` exceeded).

Root cause: the Schmidt number for these parameters is `Sc_thnf ~
883` (an ordinary base-fluid mass diffusivity, `D_f = 1e-9 m^2/s`,
against `nu_thnf ~ 8e-7 m^2/s`), which makes the concentration
boundary layer roughly 30x thinner than the momentum layer. Because
`solve_bvp` refines its mesh globally, this dragged mesh refinement
across all nine states even though only `phi` needed it.

Fix: exploit Proposition 1 of `Main.tex` directly. The momentum
sub-system (F, G, H) does not depend on theta or phi at all, so it is
solved on its own (5 states); theta and phi are then each solved as
their own separate 2-state problems using the already-known `H(eta)`.
See `base_state.py`. Converged to 6 decimal places once solved this
way.

## Finding 2: naive fresh solves at different eta_inf gave inconsistent results

A first mesh/domain-independence check solved each `eta_inf` value
independently from a fresh, generic initial guess and got
non-monotonic, inconsistent-looking results across eta_inf.

Fix: warm-start each larger domain from the previous domain's
converged solution (continuation in `eta_inf`) instead of a fresh
guess each time. This resolved the inconsistency completely -- the
naive fresh guess had been landing near a nearby-but-distinct fixed
point of Newton's method for this moderately nonlinear system at some
domain sizes.

## Finding 3: eigenvalue solve initial-guess sign error

The perturbation eigenvalue problems (momentum, thermal, species
blocks) require a normalization boundary condition, e.g. `Ff'(0) =
1`. An early implementation used an initial guess with `Ff'(0) =
-1`, the wrong sign relative to the boundary condition. This caused
every eigenvalue `solve_bvp` call to fail regardless of tolerance or
initial eigenvalue guess, and was at first misdiagnosed as numerical
stiffness. Fixed by using `eta*exp(-eta)` as the guess shape, which
satisfies the boundary condition exactly. The momentum and species
blocks additionally needed a switch to a direct shooting method
(`stability_shooting.py`), which proved more robust for this specific
linear homogeneous eigenvalue problem.

## Finding 4: eta_inf=15 was insufficient at K=0.3

A cold-start momentum solve at `K=0.3, lambda=0.5, S=0.5, eta_inf=15`
converges (low residual, no warning) to a DIFFERENT value than the
one reported throughout `Main.tex`. Both are genuine converged
solutions of the discretized, domain-truncated problem; testing at
larger `eta_inf` resolved the ambiguity -- `F'(0) = -0.442969` is
stable to 8 decimal places for `eta_inf` in {35, 40, 50}, while the
other value is an artifact of a domain too short to represent the
semi-infinite boundary layer at this K. Fixed by raising the default
`eta_inf` to 40.

A follow-up check found that even the K=0 starting point for
K-continuation was not reliably obtained by a direct cold solve at a
large `eta_inf` (different `eta_inf` gave different low-residual
answers, and one setting failed outright). Only pure rotation
(lambda=0) is textbook-unique; stretching combined with suction
reintroduces this sensitivity. Resolved with a genuine two-stage
continuation: first extend `eta_inf` incrementally from a small,
well-behaved starting domain at K=0, only then step K up to the
target at the fixed final domain. Confirmed stable from `eta_inf=8`
to `80`, reproducing `F'(0) = -0.442969` exactly. See
`solve_momentum_robust` in `base_state.py`.

## Finding 5: species-block eigenvalue is not a discrete-mode search problem

Repeated shooting attempts at a wide range of trial eigenvalues
produced far-field values that varied smoothly and never crossed
zero -- unlike the momentum and thermal blocks, which both show
clean sign changes. Tracing the far-field (eta -> infinity) limit of
the species perturbation equation analytically shows why: the two
characteristic roots of that limiting constant-coefficient equation
have product `Sc_thnf*beta + gamma`, independent of the base-state
value at infinity. Both roots are negative (decaying) for any
`gamma > -Sc_thnf*beta`, so the far-field decay condition is
satisfied by a continuum of gamma, not a discrete set -- exactly the
smooth non-crossing behaviour observed numerically. The genuine
transition is a single analytical threshold, `gamma_thr =
-Sc_thnf*beta`, not a discrete eigenvalue ladder. See
`Main.tex`'s Remark on the species stability threshold, and
`colab_ready_scripts/08_broad_search_and_species_threshold.py` and
`10_chebyshev_species_search.py` for the numerical confirmation
(including a Chebyshev spectral cross-check that independently finds
no consistent discrete eigenvalue either).

## Key finding worth restating

For K=0.3 and S in {0.5, 1.0, 1.5, 2.0, 3.0}, and more broadly across
a systematic 14-point (K,S) grid (script 08), no genuine saddle-node
(dual-solution) bifurcation was found. Apparent "folds" (natural
continuation failing to converge) coincide, in every case checked
at adequate resolution, with `max(F(eta))` approaching the analytical
singular value `1/(2K)` from `Main.tex` Eq. (19). This is reported as
an honest negative/clarifying result, not reframed to look like a
successful bifurcation study.
