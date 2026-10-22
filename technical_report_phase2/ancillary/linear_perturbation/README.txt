Scripts for the linear perturbation (Section 6) analysis, in logical order.

verify_exact_poisson.py   Symbolic verification that phi'=-4 pi G rho0 xi solves the perturbed
                          Poisson equation exactly for l=0 (sympy; log: verify_exact_poisson.log).
galerkin_fem.py            The validated Galerkin finite-element assembly (weak
                          form + boundary term + self-gravity closure).
validate_polytrope.py      Polytrope equilibrium builder, WITH the surface-
                          finding fix (auto-expands rmax until the true
                          surface is found; this was bug #2, see the report).
analytic_check2.py         Fully analytic (sympy, zero numerical differentiation)
                          identity check with a nonzero density gradient --
                          the check that definitively isolated bug #1
                          (the missing boundary term) from numerical noise.
validate_galerkin_sg.py    Full Pekeris Gamma1=4/3 validation of the Galerkin
                          solver (produces the clean sign-flip result).
apply_to_real_star.py      Applies the validated solver to the actual marginal
                          white dwarf (rho_c=5e12, neutronized EOS) --
                          produces the headline -18.68 s^-2 result.
sanity_check_real.py       Independent, non-matrix Rayleigh-quotient cross-check
                          of the above (-16.65 s^-2).
stable_crosscheck.py       Same solver on the stable control (rho_c=1e11) --
                          produces the +4.65 s^-2 cross-check.  NOTE: the profile is
                          built with rmax = 3e6 m, i.e. CUT at R = 3000 km (0.99 Msun);
                          for the complete star (R = 4246 km) the same solver gives
                          +1.756 s^-2 (lagrangian_L/logs/control_star_truncation.log).
eigenfunction_shape.py     Examines the fundamental eigenfunction's radial shape.

Dependencies: numpy, scipy, sympy. Run from this directory with PYTHONPATH=.. (the scripts that use the
equation of state import eos.py, eos_neutronized.py and initial_profile_neu.py of the parent directory).
validate_polytrope.py and validate_galerkin_sg.py only need galerkin_fem.py of this directory; the former
once imported a solver module (lawe_solver) that was not retained, and the import and the old self-test that used it were removed.
Run validate_galerkin_sg.py first to
reproduce the core validation before trusting apply_to_real_star.py's result.
