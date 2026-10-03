# Development history

This file records how the code and analysis evolved, including every numerical error that was found and corrected along the way. It is kept for transparency and is not needed to reproduce the current results; for that, see README.md. The pass numbers below are internal working numbers used while the analysis was being developed; they are not release numbers and do not correspond to the version of this archive.

## Revision log

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

A subsequent pass (18 Aug 2026) closed several smaller items
flagged across seven review reports: an explicit sensitivity
estimate for the continuation search under an alternative
cross-viscosity closure (Remark 3 in the manuscript), obtained by
log-linear interpolation against the already-tabulated lambda_c(K)
data and reported at that level of confidence -- this repository
does not yet contain a dedicated re-solve script under the
alternative closure, and the manuscript is explicit that the
Remark 3 estimate is an interpolation, not an independent
continuation run; a citation-order correction so
the two 2019 Asma et al. references render in ascending order; an
H(eta) comparison figure across the shrinking cases used in Example
4, added to give direct visual confirmation of the near-wall
sign-reversal mechanism already described in the text rather than
leaving it asserted from H(0) and H(infinity) alone; and an appendix
section that runs the repository's own symbolic-verification script
and reports its literal output against the manuscript's closed-form
momentum equations, rather than only describing the script's
existence. The manuscript does not name a specific reacting species
or fix a measured reaction-rate constant; this is stated as a
deliberate scope limitation in a new appendix note, not resolved
by an invented example.

A further pass (pass 25) addressed the two remaining MUST items from an
independent cross-check against these seven reports: (1) the
Table 6 caption previously described the broadened search as a
"3x5 grid" while the table itself lists 14 combinations (one cell
of the nominal 15-cell grid was never run); the caption and the one
remaining in-text reference to "3x5" now both say "14 systematically
tested (K,S) combinations", matching the table exactly; (2) the
Puspanathan et al. (2024) discrepancy (Section 14 of the
manuscript) is now also tested by a warm-started continuation-in-K
run (`python/colab_ready_scripts/09_puspanathan_direct_comparison.py`
extended with a K-continuation stage), stepping K from the validated
0.3 baseline up to their reported K=1.5 at fixed S=2.8, lambda=-2.4,
rather than jumping there directly. This gives F'(0)=1.718979,
matching the natural-lambda-continuation value (1.7190) to four
decimal places and still matching neither of their two reported
branch values -- reported in the manuscript as evidence that the
discrepancy is not an artefact of the continuation path taken.

**A further pass (pass 26) fixed a genuine algebra bug in the momentum
perturbation equation** (`stability.py`, `stability_shooting.py`),
found while investigating the previously-open momentum-eigenvalue-
near-boundary item. An independent symbolic (sympy) re-derivation of
the manuscript's own Eq. (linG), carried out from scratch in
`python/symbolic_derivation/derive_stability_gfpp_check.py`, found
that the coded formula for `Gfpp` was missing the term
`2*K*(Ff*G0pp - Gf*F0pp)` relative to the equation as printed in the
manuscript -- confirmed by symbolic subtraction, which is identically
zero only once the term is restored. This means the previously
reported baseline momentum eigenvalue, gamma_1^(M) ~ 0.712, was
**incorrect**. With the corrected equation, a new, properly-posed
Chebyshev collocation script
(`python/colab_ready_scripts/14_momentum_eigenvalue_chebyshev.py`)
was added: it eliminates the perturbation axial velocity Hf
algebraically (via Hf'=-2*Ff, Hf(0)=0, an exact linear integral
operator) before forming the generalized eigenvalue problem, which
avoids the singular/zero mass-matrix block that a naive
[Ff,Gf,Hf]-collocation would otherwise carry, and that was found (in
an earlier, less careful attempt) to produce spurious,
resolution-dependent near-zero "eigenvalues". The corrected,
properly-conditioned solve finds the smallest positive real
eigenvalue at the Example 1 baseline (K=0.3, lambda=0.5, S=0.5) to be
**gamma_1^(M) = 0.41857**, stable to 5 significant figures across
N=50-140 Chebyshev nodes and eta_max=20-40 -- see the script's
printed convergence table. The corrected system also reveals several
further eigenvalues closely spaced just above this one (~0.47-0.55),
a genuinely different spectral structure from the single, cleanly
isolated root the buggy equation had produced; this closer spacing is
also why naive shooting-based branch-tracking, which worked cleanly
for the (unaffected) thermal block, is markedly less reliable for the
corrected momentum block, and why the Chebyshev collocation result is
adopted as the manuscript's primary value. Both the old and new
values are positive, so the qualitative stability conclusion (the
baseline state is linearly stable) is unchanged; only the specific
number is corrected. The momentum-eigenvalue *trend* immediately
approaching the regularity boundary (as distinct from this baseline
value) was re-attempted with the corrected equations during this same
pass; base-state continuation itself became numerically fragile close
to lambda_c under the step sizes and computational budget available
in this pass, and no reliable trend was obtained -- this item remains
genuinely open, and is reported as such in the manuscript rather than
filled in with an unreliable number.

**A further pass (pass 27) resolved the momentum-eigenvalue trend over
roughly the first two-thirds of the shrinking branch**, using an
adaptive-step warm-started continuation in lambda (K=0.3, S=1.5)
cross-checked at three independent Chebyshev resolutions per point
(`python/colab_ready_scripts/15_momentum_eigenvalue_trend.py`). Over
lambda in [0,-3], the three resolutions agree to within a few percent
at every point, revealing a non-monotonic trend not observed in the
range resolved: gamma_1^(M) is non-monotonic, falling from 0.567 at lambda=0
to a local minimum of 0.0493 at lambda=-2 (still positive/stable, but
the closest approach to marginal stability found anywhere on the
branch) before rising again toward lambda=-3. Beyond lambda~-3, the
same three-resolution check shows the calculation has not converged
(10% to qualitative disagreement between resolutions), traced to
cubic-spline base-state interpolation accuracy degrading as the
base-state curvature grows near the regularity boundary; this final
segment remains honestly reported as open, exactly as before, but the
resolved two-thirds of the branch is new.

**A further pass (pass 29) eliminated all LaTeX compiler warnings** (not
just errors) from the Overleaf source -- two genuine Overfull hboxes
were fixed (a table column was too narrow; one hyphenated-compound
phrase was reworded) and cosmetic Underfull hbox noise was suppressed
via tolerance/hbadness/hfuzz settings, without reintroducing
hyphenation. Verified: zero hyphenated word-breaks anywhere in the
compiled PDF.

**A further pass (pass 30) corrected a genuine internal inconsistency**
flagged by an independent AI peer review: Proposition 2 in the
manuscript states the perturbation operator is block
*lower-triangular* (thermal and species perturbations are forced by,
but do not feed back into, the momentum block through the axial
velocity perturbation H), but an earlier revision's explicit operator
matrix and an accompanying remark had drifted into describing this as
block-*diagonal*, and the displayed matrix incorrectly showed zero
forcing entries. The operator matrix has been corrected to show H as
its own row (rather than eliminated), making the nonzero forcing of
the thermal and species rows through H explicit and the matrix
genuinely lower-triangular, matching the text. This is a
manuscript-presentation correction only -- the accompanying Python
code already solves the three eigenvalue blocks independently
(exactly as the corrected mathematics requires) and required no
change. pass 30 also reduced the manuscript from 43 to 34 pages (smaller
font/margins, resized figures, tightened list/caption spacing) with
no reduction in scientific content, at the request of the
corresponding author ahead of journal submission.

**A further pass (pass 32) fixed five reproducibility inconsistencies**
found by an independent AI peer review of pass 31, and reformatted the
manuscript for submission to International Communications in Heat and
Mass Transfer (Elsevier/elsarticle class, author-year citation style):
(1) `04_example3_eigenvalues.py` carried the pre-pass 26 buggy Gfpp formula
and still printed the old, incorrect gamma_1^(M)~0.712 comparison line;
the formula is now corrected and the script's output is explicitly
labelled as an illustrative shooting cross-check (not the manuscript's
authoritative value, which comes from script 14) with a note
explaining why a coarse shooting scan can land on the wrong one of the
corrected system's closely-spaced modes; (2) `run_all.py` described
itself as reproducing "the full workflow" when it in fact reproduces
only the core numerical results, with all extended analyses (the
14-point search, the Puspanathan comparison, the corrected Chebyshev
eigenvalue and its trend, the sensitivity checks, etc.) implemented in
separate dedicated scripts; the docstring and in-line comments now say
this explicitly, and the stale "species eigenvalue: ... open item"
print statement (the species question is resolved analytically) has
been corrected; (3) the manuscript's Table 6 caption and summary bullet
said "gap <= 1e-3" when the table's own largest entries are 1.4e-3;
both now say "gap <= 1.4e-3"; (4) the manuscript's Table 7 now states
explicitly that it shows a compact 12-row subset of the complete
18-row dataset, which is provided in `data/`; (5) both READMEs'
stated script and figure counts (previously "six scripts" / "6
figures") are corrected to the actual 16 scripts and 8 figures.

**A further pass (pass 33) fixed a genuine sign error in the species
stability perturbation equation**, found by an independent AI peer
review and confirmed by an independent symbolic re-derivation
(`python/symbolic_derivation/derive_species_sign_check.py`, validated
first against the already-correct thermal equation before being
applied to the species case): the correct linearization of the
species equation has the reaction term entering as `-Sc*beta*Phi`, not
`+Sc*beta*Phi` as an earlier revision of the manuscript printed. The
Python implementation in `stability.py` was found, on inspection, to
already use the correct sign; only the printed equation in Main.tex
and the analytical far-field argument built on it needed correcting.
`10_chebyshev_species_search.py` DID have the same sign bug as the old
printed equation and has been corrected and re-run; the corrected
far-field threshold is `gamma_thr = +Sc_thnf*beta = 176.62` (the
reverse sign from the previously reported `-176.62`). Both the
corrected shooting search (across gamma in [-50, 400]) and the
corrected Chebyshev spectral search still find no resolution-stable
discrete eigenvalue, so the manuscript's practical conclusion is
unchanged, but the analytical threshold, its derivation, and the
manuscript's narrative around it have all been corrected and are now
reported with appropriate honesty about what is and is not fully
understood, rather than asserted as a clean resolution. pass 33 also: (a)
corrected the Example 2 procedure description, which had overstated
pseudo-arclength continuation's role in the 14-point broadened search
(natural continuation with step-halving did that search; pseudo-arclength
was used as a separate diagnostic near the baseline regularity boundary
only); (b) added a genuine denominator-clipping sensitivity check
(`clip_test2.py`-style computation, see the corresponding manuscript
remark), re-running the K=0.3,S=1.5 continuation to lambda=-5.02 at
three floor tolerances spanning four orders of magnitude and obtaining
an identical result at each; (c) removed `__pycache__` directories,
added `.gitignore` and `requirements.txt` with the exact tested
package versions; (d) reformatted the manuscript for submission to
International Communications in Heat and Mass Transfer (Elsevier
`elsarticle` class, author-year citation style): trimmed the abstract
to 244 words, reduced keywords from 8 to 7, rewrote the five highlights
to fit the journal's 85-character limit, moved the appendices ahead of
the CRediT/funding/competing-interests/AI-declaration/reference-list
block to match Elsevier's own template ordering, removed the internal
audit-blue text coloring from the submission copy, and added a
Declaration of generative AI use section as the journal's guide
requires.

**Still open after pass 33, noted honestly rather than left implicit:**
the manuscript's total word count (title through the appendices,
excluding the bibliography) is approximately 18,700 words, well above
ICHMT's suggested 7,000-8,000-word length for a "Communications"
article; a dedicated compression pass was not completed in this
revision round and remains a genuine follow-up item before
submission, alongside depositing the numerical-data archive at a
citable repository such as Zenodo (GitHub alone does not satisfy a
journal's research-data-availability requirement in the way a
versioned, DOI-bearing archive does).

**A further pass (pass 36) fixed a real internal inconsistency** found by
an independent AI peer review: this repository's
`08_broad_search_and_species_threshold.py`, `run_all.py`, and
`README_python_notes.md` all still stated the species-block far-field
threshold with the OLD (pre-pass 33) sign, `gamma_thr = -Sc_thnf*beta`,
even though the manuscript itself, `stability.py`, and
`10_chebyshev_species_search.py` had already been corrected to
`gamma_thr = +Sc_thnf*beta`. All three files are now corrected and
internally consistent with the rest of the repository. Also fixed:
stale "3x5 grid" wording in script 8 and its README description
(the search is, and has been reported as, a 14-point search); the
manuscript's Table 7 previously showed only 12 of the 18 computed
Example-4 combinations, now expanded to the complete dataset that
`05_example4_table.py` already computes; the manuscript had one
leftover duplicate section heading and several instances of
elsarticle's automatic "Appendix" label being duplicated by a manual
"Appendix~" prefix in the text (both LaTeX-source bugs, now fixed);
and the eight main-text figures are now also exported as vector PDF
(`figures/*.pdf`, generated by an updated
`06_all_figures_500dpi.py` plus dedicated regeneration for the three
figures outside that script's scope) alongside the existing PNGs, and
the manuscript now includes the PDF versions directly, satisfying the
target journal's preference for vector line-art over raster images at
a fixed DPI.

**Still not done, noted honestly:** a Zenodo (or equivalent
DOI-bearing) archive of this repository has not been created; GitHub
alone does not satisfy a journal's formal research-data-deposit
requirement in the way a versioned, citable archive does, and creating
one remains an action item for the corresponding author at submission
time, not something completed in this repository as delivered.

**A further pass (pass 37) fixed the remaining documentation/wording items**
found by an independent AI peer review of pass 36, all confirmed by the
corresponding author before being changed: (1) two stale "all six
figures" / "6 figures" references in Main.tex and this script's own
comments, corrected to "eight" (the manuscript has had eight figures
since pass 27; these were leftover pre-pass 27 wording); (2) a genuine logical
inconsistency in Main.tex's continuation methodology description --
"allowing the stability of every point on both branches to be
classified" -- when no second branch exists anywhere in this study;
corrected to refer to the single computed branch; (3) `continuation.py`
docstrings and comments that described the pseudo-arclength stage as
tracing "through the fold and beyond", overstating what the reported
runs actually achieve (the manuscript itself has, since pass 33, correctly
described this as an independent diagnostic that reaches but does not
cross the regularity boundary); the code comments now match; (4) the
funding statement was changed to ICHMT's own recommended wording for
studies with no specific grant; (5) the suggested-reviewers list (in
the Overleaf package's supplementary/ folder, not part of this
repository) was revised to remove reviewers who are co-authors of the
one published study this manuscript reports an unresolved discrepancy
against, keeping the list strategically independent.

**Still not done, noted honestly:** the Zenodo/DOI-bearing archive
remains an action item for the corresponding author, not something
completed here (see the pass 36 note above). The manuscript's prose word
count (approximately 10,300 by a texcount-style measure, excluding
equations, tables, and the bibliography) remains above ICHMT's
suggested 7,000-8,000-word range for a Communication; the cover letter
in the supplementary/ folder now explicitly acknowledges this to the
Editor-in-Chief, per the journal's own stated practice for this
situation, and a small further trim was made in this pass, but a
further, more aggressive compression pass -- which would need to cut
into the Worked Examples and stability-methodology sections rather
than administrative or repeated material -- was not attempted, since
doing so risks the accuracy of verified numerical claims under the
scope of a documentation-focused revision round.

**A further pass (pass 39) closed the Zenodo research-data requirement**
flagged by an independent AI peer review of pass 38: this repository's
v1.0.0 GitHub release is now permanently archived at Zenodo with a
persistent DOI, **10.5281/zenodo.22009268** (badge added at the top of
this README), and that DOI is cited in the manuscript's Data and Code
Availability section, in `References.bib` as a `[dataset]` entry, and
in the supplementary declarations file. Also fixed in pass 39, all
confirmed by the corresponding author before being changed: (1) the
regularity-boundary claim in the governing-equations section was
slightly stronger than what had actually been proved (it asserted no
physically admissible solution can reach $F=1/(2K)$, when what is
actually established is that the explicit closed-form representation
is singular there, and that the computed branch is observed to
terminate as that point is approached -- a real, if subtle,
overclaim, now corrected); (2) the thermal-conductivity closure is now
identified explicitly as the Maxwell relation / Hamilton-Crosser
relation with spherical-particle shape factor $n=3$, and the
mass-diffusivity closure is now identified explicitly as an assumed
closure rather than phrased in a way that could read as an
experimentally established correlation; (3) the graphical abstract
was regenerated at the journal's requested 531:1328 (height:width)
aspect ratio exactly (previously 2.685:1, now 2.500:1); (4) the
Introduction's nonlinear-dynamics positioning sentence was extended
with one clause making the connection to heat- and mass-transfer
predictions explicit, addressing a scope-fit concern that the
manuscript's central contribution could read as more
dynamical-systems-flavoured than heat/mass-transfer-flavoured to an
ICHMT editor.

**A further pass (pass 40) made an author-list change and one literature
addition, both requested directly by the corresponding author rather
than flagged by review:** (1) the fourth author was changed from
S. Parvathi to Dr. Pavan Kumar Cintaginjala (Associate Professor,
Department of Basic Sciences and Humanities (Mathematics), Vignan
Institute of Technology and Science (Autonomous), Deshmukhi,
Telangana), updated in this README, the manuscript frontmatter, the
CRediT statement, and all supplementary files; (2) a directly relevant
reference was added and cited in the Introduction --
Gangadhar, Sujana Sree, Wakif and Subbarao (2024), "Stefan blowing
impact and chemical response of Rivlin-Reiner fluid through rotating
convective disk," Pramana -- Journal of Physics 98(4):160, DOI
10.1007/s12043-024-02836-w -- a single-phase Reiner-Rivlin
rotating-disk study with Stefan blowing, thermal radiation, and a
chemical reaction, co-authored by K. Subbarao. This reference was
checked before being added: it uses a different problem configuration
(non-shrinking disk, no nanofluid, no multiplicity or stability
treatment, different physical effects), so it introduces no numerical
or scientific conflict with the present manuscript's claims, and adds
a genuinely relevant, closely related citation to the literature
review. pass 40 also completed the remaining wording refinements from the
most recent independent AI peer review: the regularity-boundary claim
in the governing-equations section was tightened to state precisely
what has been proved (the explicit closed-form momentum
representation is singular at F=1/(2K); whether the original
unreduced system could cross this set via a compatibility condition
is not established), the block-triangular spectrum statement was
qualified to the well-posed/finite-dimensional-discretization case,
the Puspanathan comparison was reworded to avoid implying the two
models are demonstrably identical (the permeability difference
between the two formulations is now stated explicitly), and two
"previously unreported" priority claims about the momentum-eigenvalue
dip were softened to "not observed in the range resolved," both in
the manuscript and in the supplementary abstract files.

**A further pass (pass 43) closed the remaining documentation items**
flagged by an independent AI peer review, including a genuine
numerical error caught while addressing one of them: (1) "nine SymPy
scripts" corrected to "eleven" in the manuscript (a real undercount --
this repository's symbolic_derivation/ folder has always had 11
scripts; the manuscript text was simply never updated after the last
two were added); (2) added `data/K03_S_sweep.csv`, the K=0.3,
S-in-{0.5,1.0,2.0,3.0} continuation-failure sweep the manuscript's
Example 2 describes in prose -- computing this data freshly to build
the CSV turned up a real error in the manuscript: the S=3.0 value had
been reported as lambda_c=-1.35, but an independent, carefully
cross-checked recomputation gives lambda_c=-6.61 (S=0.5, 1.0, and 2.0
were confirmed close to their previously reported values). The
corrected, monotonic sequence -4.02, -4.40, -5.55, -6.61 (lambda_c
moving further from zero as suction S increases, which is the
physically sensible direction) is now what both the manuscript and
this CSV report; the old -1.35 value did not appear anywhere in this
repository's code or other data files, only in the manuscript prose,
so this was a manuscript-only error now corrected at the source; (3)
the References.bib Zenodo entry's `[dataset]` tag was moved from the
`note` field (which typically renders at the end of a reference) to
the front of the `title` field, matching Elsevier's convention of
`[dataset]` appearing immediately before the reference rather than
after it; (4) a separate Figure Captions document was added to the
Overleaf package's supplementary/ folder for submission systems that
request figure captions outside the main manuscript file, alongside
the existing embedded captions (which were correct already and remain
unchanged). Also trimmed a modest amount of further prose (Physical
Interpretation Framework, the consolidated-protocol and discretization
subsection openings) without touching any equation, number, or finding.

**A further pass (pass 43.4) fixed the remaining documentation and
disclosure items** flagged by an independent AI peer review of pass 43.3,
each checked before being applied: (1) `data/K03_S_sweep.csv` --
`data/momentum_eigenvalue_trend.csv`'s three "not_yet_resolved" rows
(lambda = -3.5, -4.0, -4.2) previously carried numeric gamma1_M
placeholder values alongside that status label; those values are now
blanked, since the manuscript itself treats them as unreliable and a
reviewer opening the raw CSV should not see numbers that look like
results; (2) `REFERENCE_VERIFICATION.md` still said "37 entries" in
two places, left over from an earlier reference count before three
more citations were added; corrected to 39, matching the actual
`References.bib` entry count; (3) the manuscript's AI declaration
named only Claude, but an independent reviewer (ChatGPT) also
materially contributed to this project by conducting repeated review
passes that led to real, substantive revisions throughout the pass 32-pass 43
series; the declaration now names both tools, since a disclosure that
covers only implementation assistance while omitting extensive review
assistance materially used in revising the manuscript is an
inaccurate declaration, not a stylistic choice. The reviewer's report
also named several additional AI/research tools (Jenni, SciSpace,
Consensus, PeerGenius.ai, Thesify.ai) as allegedly used in manuscript
preparation; this could not be independently verified from anything
in the available project history, so those tools were deliberately
not added to the declaration -- a disclosure should not include tools
whose use cannot actually be confirmed, any more than it should omit
tools known to have been used. Also fixed: a duplicate LaTeX label on
the thermophysical-properties table (`tab:thermophys` and
`tab:props-vals` both pointed at the same table); added explicit
in-text citations for Figures 1-4 and 8 in Results and Discussion
(Figures 5-7 were already cited there); and added an explicit
in-text citation for the Nomenclature table.


---


---

# Reference verification log of 18 August 2026 (superseded)


This note records an independent, source-by-source verification pass over
the reference list of `Main.tex` (39 entries), carried out on 18 Aug 2026
using live web search against publisher pages, DOI resolvers, and
independent citing works. It is kept here as a transparency record, in the
same spirit as the numerical-finding remarks in the manuscript itself --
not because a journal requires it, but because the manuscript's own
standard for itself is that a claim should be checked before it is relied
upon, and a citation is a claim.

## What was checked directly in this pass (highest-risk entries: newest,
least-established, or previously flagged as unverified in an earlier
review round)

All ten entries below were independently located via publisher pages
(ScienceDirect, Wiley, Springer, Emerald, MDPI, PMC) or DOI resolution,
and confirmed to match the manuscript's citation (authors, journal,
volume/pages, and DOI) exactly, unless noted:

- **Naganthran, Mustafa, Mushtaq & Nazar (2020)**, *Physica A* 556:124773,
  doi:10.1016/j.physa.2020.124773 -- confirmed exactly.
- **Mandal & Pal (2025)**, *Int. J. Numerical Methods for Heat & Fluid
  Flow* 35(8):2982-3025, doi:10.1108/HFF-02-2025-0112 -- journal, DOI,
  and author names confirmed via three independent citing sources;
  volume/page range confirmed via a Discover Nano citing article.
- **Jamrus, Ishak, Waini, Khan, Hussain, Madhukesh & Abed (2025)**, *ZAMM*
  105(2):e202300926, doi:10.1002/zamm.202300926 -- confirmed exactly,
  including the full seven-author list.
- **Paul, Patgiri & Sarma (2025)**, *World Journal of Engineering*
  22(3):612-626, doi:10.1108/WJE-01-2024-0025 -- confirmed exactly on the
  Emerald abstract page; the paper's own abstract describes a Casson
  tri-hybrid nanofluid with chemical reaction and activation energy on a
  rotating (non-shrinking) disk, matching Table 1's characterisation.
- **Sadiya & Sucharitha (2025)**, *Discover Nano*,
  doi:10.1186/s11671-025-04353-0 -- confirmed exactly via PMC full text.
- **Bartwal, Rawat, Yaseen & Pant (2025/2026)**, *J. Braz. Soc. Mech. Sci.
  Eng.* 48(1):4, doi:10.1007/s40430-025-05993-0 -- confirmed exactly.
  (Note: the journal's own cover date is 2026 for this article number;
  the manuscript's bib entry already carries `year = 2025`, matching how
  the paper is dated on ResearchGate/Google-Scholar-style listings. This
  is a cosmetic year-vs-cover-date ambiguity common to "online first"
  articles, not a wrong-paper risk, and is left as is.)
- **Algehyne, Alzubaidi, Alzahrani, Saeed & Bognár (2026)**, *Discover
  Nano* 21:41, doi:10.1186/s11671-026-04460-6 -- confirmed exactly via
  PMC full text and DOI resolution.
- **Gizewu & Ibrahim (2025)**, *J. Nonlinear Mathematical Physics* 32:72,
  doi:10.1007/s44198-025-00335-z -- confirmed exactly.
- **Asma, Othman & Muhammad (2019a)**, *Mathematics* 7(10):921,
  doi:10.3390/math7100921 -- confirmed exactly (3 authors).
- **Rauf, Mushtaq, Omar, Akif, Anwar & Shehzad (2025)**, *ZAMM*
  105(8):e70176, doi:10.1002/zamm.70176 -- journal, volume, DOI, and
  title confirmed exactly via the Wiley abstract page, **but the
  author list in the manuscript's `References.bib` had been
  shortened to four authors (Rauf, Mushtaq, Omar, Akif), dropping
  M. I. Anwar and Sabir Ali Shehzad, who are both listed on the
  publisher's own page.** This has been corrected in this revision;
  the bib entry now lists all six authors. This appears to have been
  an accidental truncation introduced during a previous editing pass
  (an earlier draft of this manuscript did carry the full six-author
  list), not a fabrication, but it is exactly the kind of small,
  independently-checkable error a citation-management tool or a
  sharp-eyed referee would catch, so it is recorded here rather than
  silently fixed.

## What was not re-verified in this specific pass

The remaining ~27 entries were not individually re-searched in this
session. Of these:

- **Puspanathan et al. (2024)**, **Mandal & Pal (2024)**, **Algehyne**
  (cross-checked again above), and the **Vitagliano & Lyons (1956)** and
  **Turkyilmazoglu & Senel (2013)** benchmark sources were independently
  verified in an earlier review round of this manuscript (documented in
  the peer-review correspondence, not reproduced here); nothing in this
  pass contradicts those earlier checks.
- The pre-1970 foundational citations (von Kármán 1921, Cochran 1934,
  Reiner 1945, Rivlin 1948, Rogers & Lance 1960, Mitschka & Ulbrecht
  1965) are textbook-standard results cited in essentially every rotating-
  disk or Reiner-Rivlin paper in this literature; they were not
  individually re-searched here, since the risk of a citation error is
  low relative to the more recent, less-indexed entries prioritised
  above, but a reader who wants that final layer of certainty should
  still check them independently before relying on this manuscript's
  reference list as a citation source in its own right.
- The remaining mid-list entries (Mustafa et al. 2016; Tabassum &
  Mustafa 2018; Lv et al. 2021; Sahoo & Kumar 2020; Mousavi et al. 2019;
  Elnaqeeb et al. 2021; Manjunatha et al. 2022; Merkin 1986; Weidman
  et al. 2006; Harris et al. 2009; Turkyilmazoglu 2012, 2014; Tu et al.
  2021; Boyd 2001; Shevchuk 2009) were likewise not individually
  re-searched in this specific pass.

## Bottom line

Ten of the manuscript's most recently added and highest-risk citations
were independently confirmed against publisher records in this pass, one
genuine (non-fabrication) error was found and corrected, and roughly
three-quarters of the full 39-entry list remain unverified by this
specific session -- stated plainly here rather than rounded up to "the
references have been checked." A full DOI-by-DOI pass over every
remaining entry is a mechanical task well suited to a reference manager
(Zotero/EndNote DOI lookup, or a Crossref batch query) and is recommended
as a final pre-submission step, exactly as several of the peer-review
reports already suggested.


## Version 2.1.0

- The conditioning statement was narrowed. The bordered Newton Jacobian replaces the
  wall equation by F(0) = lambda, so it cannot certify the singular coefficient; the
  text now reports only that no progressive deterioration is observed.
- `run_conditioning_check.py` now also reports the coefficient-level quantity
  min |1 - 2KF| (the smallest singular value of the matrix multiplying (F'', G'')) and
  writes `data/conditioning_check.csv`.
- Correction: the wall curvature F''(0) quoted earlier at lambda = -1.66 (-76.7, "resolved
  to 0.2%") does not reproduce across N = 160, 200, 240 (about -74 to -77, spread 3.9%).
  It is resolved to 0.2% only down to lambda = -1.6 (F''(0) = -40.4).
- Graphical abstract regenerated by `python/make_graphical_abstract.py`; K = 0 is now
  identified as the Newtonian control rather than counted among the K values tested.
- README opening and sign statement made neutral and exact.


## Version 2.2.0

- Bibliography corrected and completed: first author of the Discover Nano paper is A. B. Sadiya;
  the Mitschka & Ulbrecht title, issue and DOI are now exact; the Rogers & Lance key matches the
  1960 year; DOIs added to sixteen entries that lacked them; the Zenodo entry no longer carries a
  version label.
- `REFERENCE_VERIFICATION.md` now covers all 43 entries (the earlier log said 42 and had omitted one).
- Wording made to match what was computed: the time-integration checks confirm the unstable
  regime on both branches; the zero crossings themselves rest on the eigensolve and bisection.
- `run_conditioning_check.py` runs without arguments (default K = -0.3).
- No numerical result changed. The Newtonian validation, the two stability crossings
  (-0.4469 and -0.6036), the time-integration growth rate and the stored negative-K, suction and
  conditioning data were re-run or re-read and agree with the manuscript.
