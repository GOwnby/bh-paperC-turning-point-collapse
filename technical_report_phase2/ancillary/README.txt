Scripts underlying the Phase 2 progress report, in the order they were used.

eos.py                    Chandrasekhar degenerate-electron-gas EOS + constants.
initial_profile.py        Newtonian hydrostatic equilibrium with the pure EOS;
                           reproduces the Chandrasekhar mass (Sec. 3.1).
neutronization.py         SEMF-based beta-equilibrium Y_e(rho) (Sec. 3.2).
eos_neutronized.py        Fast tabulated neutronized EOS built from the above.
initial_profile_neu.py    Hydrostatic equilibrium with the neutronized EOS;
                           shows the mass-radius turnover (Sec. 3.2).
collapse2.py               Builds inhomogeneous initial data on a Lagrangian grid.
evolve2.py                Misner-Sharp evolution RHS (Sec. 2), pure-EOS version.
                           (The neutronized-EOS version, evolve_neu.py, and its
                           companion collapse_neu.py, are the same code with
                           eos_neutronized substituted for eos/initial_profile;
                           omitted here for brevity, available on request.)
euler_hydro.py             Eulerian conservative hydrodynamics grid class (Sec. 5).
euler_run.py               Full Eulerian evolution (produces the Sec. 5.3 numbers).
euler_convergence.py       The three-resolution convergence check (Sec. 5.3 table) --
                           this is the most important script in this directory.

Dependencies: numpy, scipy. No exotic packages required.

--- Added in the extended resolution study ---
euler_module.py     Refactored, checkpoint-capable Eulerian setup shared by all
                     resolution runs; tracks mass conservation alongside central
                     density and surface radius.
study1.py            Validates the refactor against N=800, 1600 (matches the
                     original convergence table).
study2.py            Pushes to N=3200, 6400 (the extended, still-diverging table
                     in Sec. 5.3).
control_test.py      The critical control: same code, same resolutions, but a
                     robustly STABLE configuration (rho_c=1e11 kg/m^3). Converges
                     properly, ruling out a generic implementation bug and
                     localizing the non-convergence to the marginal configuration
                     specifically.
