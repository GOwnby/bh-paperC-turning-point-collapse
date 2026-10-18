Scripts for Section 7 (eigenfunction-seeded Eulerian evolution), in order.

get_eigenfunction.py      Re-solves the Galerkin eigenvalue problem (Section 6)
                          and saves the normalized fundamental eigenfunction
                          to eigenfunction_data.pkl.
euler_module_seeded.py    Eulerian hydro setup (same physics as
                          ../euler_module.py) but with the initial density
                          perturbed by the saved eigenfunction at a controlled
                          amplitude epsilon, instead of relying on grid
                          truncation noise. Includes the R_star rescaling fix
                          and the near-center regularization fix described in
                          the report.
run_seeded2.py             Runs the controlled test at N=800 and N=1600,
                          epsilon=1e-4 -- produces the headline divergence
                          result (rarefaction vs collapse).
run_seeded_strong.py       The epsilon=1e-2 (100x stronger) check that drove
                          both resolutions into an unphysical regime --
                          included for completeness/reproducibility, not
                          because it is a trustworthy physical result.

Run order: get_eigenfunction.py first (produces eigenfunction_data.pkl),
then run_seeded2.py. Depends on the base EOS/hydro modules in ../ (eos.py,
eos_neutronized.py, initial_profile_neu.py, euler_hydro.py) and on
galerkin_fem.py, validate_polytrope.py from ../linear_perturbation/.

--- Added: the pure-growing-mode seed (Section 7.4) ---
verify_pure_mode.py        Symbolic derivation/verification of the pure-
                          growing-mode initial condition v(r,0)=sigma*xi(r,0).
euler_module_seeded_v2.py  Eulerian setup with BOTH density and velocity
                          seeded consistently (no decaying-mode admixture),
                          as opposed to euler_module_seeded.py's density-only
                          (zero-velocity) seed.
run_pure_mode.py            Runs the pure-growing-mode test at N=800, 1600 --
                          produces the Section 7.4 result (same qualitative
                          divergence as the simpler seed).
