Scripts for the linear perturbation (Section 6) analysis, in logical order.

verify_exact_poisson.py   Symbolic verification that phi'=-4 pi G rho0 xi solves
                          the perturbed Poisson equation exactly for l=0.
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
                          produces the +4.65 s^-2 cross-check.
eigenfunction_shape.py     Examines the fundamental eigenfunction's radial shape.

Dependencies: numpy, scipy, sympy. Run validate_galerkin_sg.py first to
reproduce the core validation before trusting apply_to_real_star.py's result.
