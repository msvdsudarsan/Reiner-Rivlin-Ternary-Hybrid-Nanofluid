# Exploratory Kremer–Grest NEMD check (LAMMPS)

Generic bead–spring melt (LJ units, chains of fixed length, FENE bonds, WCA pair interaction), SLLOD with `fix deform` at four shear rates (0.02, 0.05, 0.1, 0.2), 300 000 steps each, single seed (11).
Files: `in.equil`, `in.shear`, `run_all.sh`, `kg_chains.data`, `kg_equil.data`, logs, `stress_r*.dat`, `analyze.py`, `kg_results.json`.
Conventions: sigma = -P, N1 = sxx-syy, N2 = syy-szz, eta = sxy/rate, Psi = N/rate^2. First third of each run discarded; 5 block averages give the standard errors.

Scope: exploratory only. This is a generic polymer melt, not PDMS or any fluid in the paper, and Psi2 depends strongly on rate, so it does not validate the Reiner–Rivlin closure. It only shows that the sign (N2 < 0, N2/N1 about -0.12 to -0.17) and order of magnitude of a second normal-stress difference can be obtained from a molecular simulation. Single seed, short runs, no finite-size study.
