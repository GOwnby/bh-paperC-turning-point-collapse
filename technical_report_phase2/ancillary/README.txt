Scripts underlying the Phase 2 progress report, in the order they were used.

eos.py                    Chandrasekhar degenerate-electron-gas EOS + constants.
initial_profile.py        Newtonian hydrostatic equilibrium with the pure EOS;
                           reproduces the Chandrasekhar mass (Sec. 3.1).
neutronization.py         SEMF-based beta-equilibrium Y_e(rho) (Sec. 3.2).
eos_neutronized.py        Fast tabulated neutronized EOS built from the above.
initial_profile_neu.py    Hydrostatic equilibrium with the neutronized EOS;
                           shows the mass-radius turnover (Sec. 3.2).
collapse2.py               Builds inhomogeneous initial data on a Lagrangian grid.  FIRST VERSION of the
                           Lagrangian code (see the note after this list).
evolve2.py                Misner-Sharp evolution RHS (Sec. 2), pure-EOS version.  FIRST VERSION: Sec. 4 of
                           the report shows two errors in it (dU/dt with a^2 instead of a; R'/X frozen at its
                           initial value) and that the initial data are not an equilibrium.
                           (The neutronized-EOS version, evolve_neu.py, and its
                           companion collapse_neu.py, were the same code with
                           eos_neutronized substituted for eos/initial_profile.
                           These two files were NOT retained and are NOT part of
                           this archive; the Lagrangian neutronized-EOS runs of
                           the first version of the report cannot be rerun from it.
                           They are replaced by lagrangian_L/, see below.)
lagrangian_envelope.py     FIRST-VERSION re-creation (cubic-spline derivatives; not an original file) of the lost
                           shell-crossing driver:  python3 lagrangian_envelope.py pure|neu 1e-3 ; logs
                           lagrangian_envelope_pure.log, lagrangian_envelope_neu.log.  It did NOT reproduce the 14.2 / 9.1
                           of the first version, and Sec. 4 of the report shows that its crossings (14.783, pure, shell 0;
                           2.188, neutronized, shell 299) are a DISCRETISATION ARTEFACT: an exact matrix replica reproduces
                           them and traces them to the initial data (acceleration -0.77 g / -0.35 g at node 0, -1.0 g at
                           the outer node) and to the vacuum edge.  No result of the report depends on it.
euler_hydro.py             Eulerian conservative hydrodynamics grid class (Sec. 5).
euler_run.py               Full Eulerian evolution (produces the Sec. 5.3 numbers).
euler_convergence.py       The three-resolution convergence check (Sec. 5.3 table) --
                           this is the most important script in this directory.

Dependencies: numpy, scipy (sympy only for the symbolic checks). No exotic packages required.
Running: the scripts import each other by module name. Run a script of this directory from this
directory (cd ancillary; python3 euler_convergence.py); for the scripts in the subdirectories
linear_perturbation/ and seeded_evolution/ run from inside the subdirectory with PYTHONPATH=..
(for example: cd linear_perturbation; PYTHONPATH=.. python3 stable_crosscheck.py).
Logs: the original archive (Zenodo v1.0.1) contained no logs of these scripts. The <script>.log files in this directory and in the subdirectories are the
output of the scripts as run for the verification (verify.sh --full; see ../../PROVENANCE.md for which file differs from the original and why).

--- Added in the extended resolution study ---
euler_module.py     Refactored, checkpoint-capable Eulerian setup shared by all
                     resolution runs; tracks mass conservation alongside central
                     density and surface radius.
study1.py            Validates the refactor against N=800, 1600 (matches the
                     original convergence table).
study2.py            Pushes to N=3200, 6400 (the extended, still-diverging table
                     in Sec. 5.3).
control_test.py      The critical control: same code, same resolutions, but a
                     STABLE-branch configuration (rho_c=1e11 kg/m^3). Converges
                     properly, which argues against a generic implementation bug and
                     localizes the non-convergence to the marginal configuration.
                     NOTE: its profile is built with rmax = 3e6 m, i.e. the star is
                     CUT at R = 3000 km (M = 0.99 Msun, density still 3.6 % of rho_c
                     there; the complete star has R = 4246 km). The same cut is in
                     linear_perturbation/stable_crosscheck.py (+4.65 s^-2; +1.756 s^-2
                     for the complete star): lagrangian_L/logs/control_star_truncation.log.

--- Added after version 1.0.1: the convergence-verified Lagrangian re-do (Sec. 4 of the report) ---
lagrangian_L/              Staggered conservative Misner-Sharp scheme (src/), tests (a)-(h) of Sec. 4, logs/, results/*.json,
                           results/RESULTS.md, derivations.md, ref/ (byte-identical copies of the report's routines used by it;
                           sha256 in logs/ref_checksums.log), run_all.sh (regenerates everything: about 20 min on two cores,
                           numpy + scipy only; 1177 s in the shipped run).  Its README.md gives the layout and the findings.
                           NEW WORK, not part of Zenodo v1.0.1: generated by an AI assistant, independent of the lost original
                           driver; see ../../PROVENANCE.md, section 7, for what it replaces and for its limits.  It replaces
                           evolve2.py / collapse2.py / lagrangian_envelope.py (first version) as the basis of Sec. 4.
lagrangian_L_extra_checks.py   (+ .log) Statements of Sec. 4 that no log of lagrangian_L contains: first-step kick of an
                           unfloored surface zone for N = 100 ... 1600; 1 - a at the centre; zones with R' <= 0 at the fixed-step
                           'core crossing'.  Run from this directory:  python3 lagrangian_L_extra_checks.py lagrangian_L  (40 s).
