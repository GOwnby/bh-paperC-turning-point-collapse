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

--- Notes added for the verification (verify.sh) ---
euler_module_seeded.py           is the ORIGINAL module of the archive (byte-identical to Zenodo v1.0.1). Its run() has no stop test other than NaN / negative density;
                                 after the central density has run away the density floor of the conservative-to-primitive step keeps rho > 0, so the run goes on with
                                 dt -> 0 (garbage state, e.g. v_c = -1.9e43 m/s) until the wall-clock budget is used up. The collapse time t = 0.618 s quoted for N = 1600 is
                                 the time reached when that happened (status WALL_LIMIT), i.e. the time at which the state became meaningless, not a measured event.
euler_module_seeded_guarded.py   the same module with ONE added condition: run() stops (status BLOWUP) as soon as max|v| > 1e12 m/s (c = 3e8 m/s). Everything before that
                                 point is unchanged (checked: same trajectory, N = 800, 1600, eps = 1e-4 and 1e-2; see ../../../PROVENANCE.md). verify.sh uses it so that the
                                 end of the run does not depend on the speed of the machine.
run_seeded2.py, run_seeded_strong.py   choose the module with the environment variable SEEDED_MODULE (default: euler_module_seeded, the original);
                                 verify.sh sets SEEDED_MODULE=euler_module_seeded_guarded and WALL_BUDGET=100000.
run_seeded2.log, run_seeded_strong.log                  output of the ORIGINAL module (default wall budgets 150 s and 120 s per resolution; the last lines of the
                                 N = 1600 run show the meaningless state reached after the run-away and depend on the speed of the machine).
run_seeded2_guarded.log, run_seeded_strong_guarded.log  output with the guarded module (what verify.sh checks).
