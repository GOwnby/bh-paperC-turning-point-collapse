# Verification report: bh-paperC-turning-point-collapse

Mode: **quick**.  Date: 2026-10-22.  Total run time of this invocation: 0.2 min.

**886 PASS, 0 FAIL, 0 SKIPPED** of 886 checks.

## Environment

- python: 3.10.12
- platform: Linux-6.8.0-138-generic-x86_64-with-glibc2.35
- machine: x86_64
- numpy: 2.2.6
- scipy: None
- sympy: 1.14.0
- mpmath: 1.3.0
- matplotlib: 3.10.9
- gcc: gcc (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0

## Checks

Tolerances: half a unit in the last quoted digit unless stated; the reason is given in `expected.json`. "carried" marks a check that cannot run on the machine of this invocation (missing module) and shows its last result.

| id | where | quantity | paper | computed | tol | status |
|---|---|---|---|---|---|---|
| SYM-sympy_eos_checks | Sec. 2.1, 2.3, 5.1, 5.2, 5.4, 5.7 | EOS identities: Fermi integrals, f', Gamma1 limits 5/3 and 4/3, first law of eps (exact), failure for eps=rho c^2, Salpeter lattice, du/drho=P/rho^2, linearisation 3P(2-3Gamma1), g_eff identity, nu |  | 15 PASS, 0 FAIL lines |  | **PASS** |
| SYM-sympy_weak_form | Sec. 3.2 (Eq. ACD, Eq. boundary-term) | weak form A,C,D from the strong form, symmetry, boundary term |  | 4 PASS, 0 FAIL lines |  | **PASS** |
| SYM-sympy_fold_scaling | Sec. 5.6-5.7 | Painleve-I scaling, slow-manifold -1/8s^2, Boutroux map, overdamped scaling, Airy law |  | 10 PASS, 0 FAIL lines |  | **PASS** |
| SYM-verify_exact_poisson | Sec. 3.1 (Eq. exact-poisson-main) | phi' = -4 pi G rho0 xi solves the perturbed Poisson equation exactly (residual 0) |  | residual printed: 0 |  | **PASS** |
| SYM-analytic_check2 | Sec. 3.2 (analytic test case) | weak form = strong form for rho0 = rho_c(1-(r/R)^2), residual 0 |  | LHS_exact - (weak_no_boundary + boundary_term) = 0 |  | **PASS** |
| SYM-analytic_check_uniform | Sec. 3.2 (analytic test case with boundary term) | weak form = strong form for uniform density, only with the boundary term |  | 2 PASS |  | **PASS** |
| SYM-verify_pure_mode | tech. report Sec. 7.4 | pure growing mode v = sigma xi |  | residual 0; C2=0 <=> v0 = sigma xi0 |  | **PASS** |
| TOV-M | Sec. 2.3 | M at rho_c = 5e13 kg/m^3 [Msun] | 1.4143 | 1.41435 | 5e-05 | **PASS** |
| TOV-R | Sec. 2.3 | R [km] | 815.8 | 815.817 | 0.05 | **PASS** |
| TOV-C | Sec. 2.3 | compactness 2GM/(Rc^2) | 0.00512 | 0.00512 | 6e-06 | **PASS** |
| TOV-G1c | Sec. 2.3, 3.3 | Gamma1 at the centre | 1.33409 | 1.33409 | 5e-06 | **PASS** |
| TOV-G1min | Sec. 3.3 | Gamma1 min over the star (= centre value) | 1.33409 | 1.33409 | 5e-06 | **PASS** |
| FIG5-peak | Sec. 2.3, Fig. 5 | M_max of the TOV sequence [Msun] | 1.4166 | 1.41661 | 5e-05 | **PASS** |
| FIG5-rhopeak | Sec. 2.3 | rho_c at the mass peak [1e13 kg/m^3] | 2.31 | 2.31057 | 0.005 | **PASS** |
| CT-1.4176 | Sec. 2.3 | agreement with Chandrasekhar-Tooper 1.4176 [%] | 0.07 | 0.0697331 | 0.005 | **PASS** |
| NEWT-1.456 | Sec. 2.1 | Newtonian hydrostatic mass at rho_c = 1e20 (asymptote) [Msun] | 1.456 | 1.4563 | 0.0005 | **PASS** |
| NEWT-noexceed | Sec. 2.1 | Newtonian sequence never exceeds its asymptote (Chandrasekhar mass) |  | max over sequence 1.45630 |  | **PASS** |
| NEWT-analytic | Sec. 2.1 | closed-form Chandrasekhar mass from Lane-Emden n=3 (w_3 sqrt(3pi)/2 (hbar c/G)^1.5/(mu_e m_u)^2) [Msun] | 1.456 | 1.4563 | 0.0005 | **PASS** |
| EOS-G1-NR | Sec. 2.1 | Gamma1(x->0) | 1.6667 | 1.66667 | 5e-05 | **PASS** |
| EOS-G1-UR | Sec. 2.1 | Gamma1(x->oo) | 1.3333 | 1.33333 | 5e-05 | **PASS** |
| EOS-Pc2 | Sec. 2.3, 3.6 | P/(rho c^2) at rho_c = 5e13 | 0.002 | 0.0020205 | 0.0005 | **PASS** |
| EOS-epsrho | Sec. 3.6 | eps/(rho c^2) - 1 at rho_c = 5e13 | 0.006 | 0.0058011 | 0.0005 | **PASS** |
| EOS-margin-fid | Sec. 2.3 | Gamma1 - 4/3 at the centre, rho_c = 5e13 | 0.00076 | 0.00075842 | 5e-06 | **PASS** |
| EOS-margin-max | Sec. 3.6 | Gamma1 - 4/3 at rho_max (upper end of (0.8-1.3)e-3) | 0.0013 | 0.001262 | 5e-05 | **PASS** |
| PEK-sign | Sec. 3.2, Fig. 1 | fundamental mode unstable for Gamma1 < 4/3 and stable above (13 polytropes) |  | 13 points |  | **PASS** |
| PEK-zero | Sec. 3.2 | omega^2 at Gamma1 = 4/3 [s^-2] | 2e-11 | 2.18582e-11 | 1e-11 | **PASS** |
| PEK-scale | Sec. 3.2 | |omega^2| of order 1e-8 s^-2 away from the critical point (and >= 3 orders above the value at 4/3) |  | min away from 4/3: 6.14e-09 |  | **PASS** |
| WSI-good | Sec. 3.2 | weak form = strong form on polytropes for Gamma1 = 1.25 ... 5/3, two trial functions (relative to the integrand scale <= 5e-6, i.e. >= 5-6 digits) |  | max (S-B)/scale = 1.0e-06 |  | **PASS** |
| WSI-trunc | Sec. 3.2 (error 2) | the same identity FAILS on a polytrope truncated before its surface (> 1e-4 of the scale) |  | see log |  | **PASS** |
| PEK-overtones | Sec. 3.2 | overtones correctly ordered (node counts 0,1,2,3,4 for increasing eigenvalue) and positive, fundamental negative for Gamma1 < 4/3 and positive above, for the six polytropes Gamma1 = 1.25 ... 5/3 (validate_galerkin_sg; the row Gamma1 = 1.2 is the n = 5 polytrope, which has no finite surface, and is correctly reported NOT FOUND) |  | 6 polytropes |  | **PASS** |
| RR-n1 | Sec. 3.6, Fig. 6 | Rayleigh-Ritz sigma^2 with 1 basis function(s) [s^-2] | -27.85 | -27.8524 | 0.005 | **PASS** |
| RR-n2 | Sec. 3.6, Fig. 6 | Rayleigh-Ritz sigma^2 with 2 basis function(s) [s^-2] | -28.52 | -28.525 | 0.005 | **PASS** |
| RR-n3 | Sec. 3.6, Fig. 6 | Rayleigh-Ritz sigma^2 with 3 basis function(s) [s^-2] | -28.55 | -28.5517 | 0.005 | **PASS** |
| RR-n4 | Sec. 3.6, Fig. 6 | Rayleigh-Ritz sigma^2 with 4 basis function(s) [s^-2] | -28.55 | -28.5528 | 0.005 | **PASS** |
| RR-n5 | Sec. 3.6, Fig. 6 | Rayleigh-Ritz sigma^2 with 5 basis function(s) [s^-2] | -28.55 | -28.5534 | 0.005 | **PASS** |
| RR-xi-r | Sec. 3.5 | trial xi = r: sigma^2 | -27.85 | -27.8524 | 0.005 | **PASS** |
| RR-xi-rea | Sec. 3.5 | trial xi = r e^a: sigma^2 | -27.67 | -27.6712 | 0.005 | **PASS** |
| RR-symm | Sec. 3.6 | matrices A, B symmetric to machine precision (relative asymmetry < 1e-12) |  | max asymmetry 2.0e-16 |  | **PASS** |
| RR-sigma | Sec. 3.6 | sigma = sqrt(-sigma^2) [1/s] | 5.343 | 5.34354 | 0.0015 | **PASS** |
| RR-efold | Sec. 3.6 | e-folding time 1/sigma [s] | 0.187 | 0.187142 | 0.0005 | **PASS** |
| RR-cmp-12k | Sec. 3.7 | sigma^2 at rho_c=5e13, n=5, 12000 pts [s^-2] | -28.553 | -28.553 | 0.0005 | **PASS** |
| RR-grid0 | Sec. 3.7 | sigma^2(5e13), n=5, 3000 points [s^-2] | -28.511 | -28.511 | 0.0005 | **PASS** |
| RR-grid1 | Sec. 3.7 | sigma^2(5e13), n=5, 6000 points [s^-2] | -28.546 | -28.546 | 0.0005 | **PASS** |
| RR-grid2 | Sec. 3.7 | sigma^2(5e13), n=5, 12000 points [s^-2] | -28.553 | -28.553 | 0.0005 | **PASS** |
| RR-grid3 | Sec. 3.7 | sigma^2(5e13), n=5, 24000 points [s^-2] | -28.555 | -28.555 | 0.0005 | **PASS** |
| RR-extrap | Sec. 3.7 | Richardson extrapolation (O(N^-2)) of the 12000/24000-point values [s^-2] | -28.556 | -28.5557 | 0.0006 | **PASS** |
| TP-idea-peak | Sec. 3.7, table | ideal (mu_e=2): rho_c at dM/drho_c = 0 [1e13 kg/m^3] | 2.3106 | 2.3106 | 5e-06 | **PASS** |
| TP-idea-zero | Sec. 3.7, table | ideal (mu_e=2): rho_c at sigma^2 = 0 [1e13 kg/m^3] | 2.31074 | 2.31074 | 5e-06 | **PASS** |
| TP-idea-diff | Sec. 3.7, table | ideal (mu_e=2): difference [%] | 0.006 | 0.0061 | 0.0007 | **PASS** |
| TP-idea-Mmax | Sec. 3.7, 5.5 | ideal (mu_e=2): M_max [Msun] | 1.4166 | 1.41661 | 5e-05 | **PASS** |
| TP-He4-peak | Sec. 3.7, table | He4: rho_c at dM/drho_c = 0 [1e13 kg/m^3] | 2.34621 | 2.34621 | 5e-06 | **PASS** |
| TP-He4-zero | Sec. 3.7, table | He4: rho_c at sigma^2 = 0 [1e13 kg/m^3] | 2.34635 | 2.34635 | 5e-06 | **PASS** |
| TP-He4-diff | Sec. 3.7, table | He4: difference [%] | 0.006 | 0.0061 | 0.0007 | **PASS** |
| TP-He4-Mmax | Sec. 3.7, 5.5 | He4: M_max [Msun] | 1.4013 | 1.40134 | 5e-05 | **PASS** |
| TP-C12-peak | Sec. 3.7, table | C12: rho_c at dM/drho_c = 0 [1e13 kg/m^3] | 2.38557 | 2.38557 | 5e-06 | **PASS** |
| TP-C12-zero | Sec. 3.7, table | C12: rho_c at sigma^2 = 0 [1e13 kg/m^3] | 2.38572 | 2.38572 | 5e-06 | **PASS** |
| TP-C12-diff | Sec. 3.7, table | C12: difference [%] | 0.006 | 0.0062 | 0.0007 | **PASS** |
| TP-C12-Mmax | Sec. 3.7, 5.5 | C12: M_max [Msun] | 1.3849 | 1.3849 | 5e-05 | **PASS** |
| TP-O16-peak | Sec. 3.7, table | O16: rho_c at dM/drho_c = 0 [1e13 kg/m^3] | 2.40187 | 2.40187 | 5e-06 | **PASS** |
| TP-O16-zero | Sec. 3.7, table | O16: rho_c at sigma^2 = 0 [1e13 kg/m^3] | 2.40202 | 2.40202 | 5e-06 | **PASS** |
| TP-O16-diff | Sec. 3.7, table | O16: difference [%] | 0.006 | 0.0061 | 0.0007 | **PASS** |
| TP-O16-Mmax | Sec. 3.7, 5.5 | O16: M_max [Msun] | 1.3782 | 1.37823 | 5e-05 | **PASS** |
| TP-Fe56-peak | Sec. 3.7, table | Fe56: rho_c at dM/drho_c = 0 [1e13 kg/m^3] | 2.93211 | 2.93211 | 5e-06 | **PASS** |
| TP-Fe56-zero | Sec. 3.7, table | Fe56: rho_c at sigma^2 = 0 [1e13 kg/m^3] | 2.9323 | 2.9323 | 5e-06 | **PASS** |
| TP-Fe56-diff | Sec. 3.7, table | Fe56: difference [%] | 0.007 | 0.0065 | 0.0007 | **PASS** |
| TP-Fe56-Mmax | Sec. 3.7, 5.5 | Fe56: M_max [Msun] | 1.1507 | 1.15071 | 5e-05 | **PASS** |
| TP-abstract | Abstract, Sec. 3.7, Sec. 8 | turning-point agreement for five EOS (largest difference) [%] | 0.007 | 0.0065 | 0.0006 | **PASS** |
| TP-grid3000 | Sec. 3.7 | sigma^2(rho_max), n=5, 3000 points [s^-2] | 0.021 | 0.02129 | 0.0005 | **PASS** |
| TP-grid12000 | Sec. 3.7 | sigma^2(rho_max), n=5, 12000 points [s^-2] | 0.0012 | 0.00115 | 0.0001 | **PASS** |
| TP-grid48000 | Sec. 3.7 | sigma^2(rho_max), n=5, 48000 points [s^-2] | 0.0001 | 0.00011 | 6e-05 | **PASS** |
| TP-basis | Sec. 3.7 | varying the basis size from 3 to 7 functions changes sigma^2(rho_max) negligibly (all within 1.2e-3 s^-2 of each other) |  | [0.00218, 0.00131, 0.00115, 0.00115, 0.00114] |  | **PASS** |
| TP-3000-offset | Sec. 3.7 | zero-crossing offset from rho_max at 3000 points [%] | 0.11 | 0.114 | 0.005 | **PASS** |
| FIRST-zero | Sec. 3.7 | first implementation (eps=rho c^2, 8piG term off): sigma^2 = 0 at [1e13 kg/m^3] | 1.86 | 1.8601 | 0.005 | **PASS** |
| FIRST-peak | Sec. 3.7 | first implementation: mass peak at [1e13 kg/m^3] | 3.56 | 3.5553 | 0.005 | **PASS** |
| FIRST-sign | Sec. 3.7 | opposite sign of the 8piG term: zero crossing at [1e13] | 1.57 | 1.5683 | 0.005 | **PASS** |
| FIRST-noeps | Sec. 3.7 | 8piG term on, eps = rho c^2: zero crossing at [1e13] | 2.2952 | 2.2952 | 5e-05 | **PASS** |
| FIRST-both | Sec. 3.7 | both corrections: zero crossing at [1e13] | 2.3107 | 2.3107 | 5e-05 | **PASS** |
| FIRST-sigma-old | Sec. 3.7 | first implementation: sigma^2(5e13) [s^-2] | -41.09 | -41.092 | 0.005 | **PASS** |
| FIRST-efold-old | Sec. 3.7 | first implementation: e-folding time [s] | 0.156 | 0.155999 | 0.0005 | **PASS** |
| FIRST-Mmax-old | Sec. 2.3 | earlier version with eps = rho c^2: M_max [Msun] | 1.4252 | 1.4254 | 0.0002 | **PASS** |
| AP-s2 | Sec. 3.7 (independent solver) | shooting A: sigma^2 at rho_c = 5e13 [s^-2] | -28.556 | -28.5559 | 0.001 | **PASS** |
| AP-s2-B | Sec. 3.7 | shooting B (linearised Euler + Einstein constraints): sigma^2 at 5e13 | -28.556 | -28.5559 | 0.001 | **PASS** |
| AP-AB | Sec. 3.7 | two independent shooting formulations agree to better than 1e-4 relative at all five densities |  | max rel diff 0.0e+00 |  | **PASS** |
| AP-nodes | Sec. 3.7 | fundamental mode without nodes at every density tested |  | ['0', '0', '0', '0', '0'] |  | **PASS** |
| AP-zero | Sec. 3.7 | zero of sigma^2 (shooting A) [1e13 kg/m^3] | 2.3106 | 2.3106 | 2e-06 | **PASS** |
| AP-peak | Sec. 3.7 | mass peak (independent TOV, mscan.py) [1e13 kg/m^3] | 2.3106 | 2.3106 | 3e-06 | **PASS** |
| AP-agree | Sec. 3.7 | zero of sigma^2 and mass peak agree to 1e-6 relative |  | rel diff 1.0e-07 |  | **PASS** |
| AP-RR1 | Sec. 3.7 | independent solver with the RR basis, 1 functions [s^-2] | -27.855 | -27.8549 | 0.0005 | **PASS** |
| AP-RR2 | Sec. 3.7 | independent solver with the RR basis, 2 functions [s^-2] | -28.528 | -28.5276 | 0.0005 | **PASS** |
| AP-RR3 | Sec. 3.7 | independent solver with the RR basis, 3 functions [s^-2] | -28.554 | -28.5542 | 0.0005 | **PASS** |
| AP-RR5 | Sec. 3.7 | independent solver with the RR basis, 5 functions [s^-2] | -28.5559 | -28.5559 | 5e-05 | **PASS** |
| AP-conv | Sec. 3.7 (independent solver) | spread of the shooting root under variations of r0, rtol, x_end, relative (upper bound) | 0.0001 | 1.57158e-07 | 5e-05 | **PASS** |
| AP-fem | Sec. 3.7 | finite-element eigenvalue, Richardson extrapolation of N = 800, 1600 [s^-2] | -28.556 | -28.5558 | 0.001 | **PASS** |
| AP-fem-agree | Sec. 3.7 | finite-element result agrees with shooting to < 1e-4 relative (Richardson-extrapolated) |  | rel diff 6.2e-06 |  | **PASS** |
| NW-ratio1 | Sec. 3.7 (Newtonian limit, independent solver) | sigma^2_GR/sigma^2_Newt - 1 at compactness 2.7e-5 | -0.0001 | -9.797e-05 | 2e-05 | **PASS** |
| NW-ratio2 | Sec. 3.7 | ratio GR/Newt at compactness 1.1e-4 | 0.9996 | 0.999556 | 0.00015 | **PASS** |
| FIG5-ratio1 | Sec. 3.7, Fig. 5 left | sigma^2_GR/sigma^2_Newt at compactness 2.7e-5 | 0.9999 | 0.999894 | 5e-05 | **PASS** |
| FIG5-ratio2 | Sec. 3.7, Fig. 5 left | ratio at compactness 1.1e-4 | 0.9996 | 0.999551 | 5e-05 | **PASS** |
| FIG5-trend | Sec. 3.7 | ratio decreases steadily with compactness (departure grows) |  | 11 points, 0.9999 ... 0.5398 |  | **PASS** |
| SALP-C12 | Sec. 3.7, 5.5 | C12 Salpeter-GR M_max [Msun] | 1.3849 | 1.3849 | 5e-05 | **PASS** |
| SALP-C12-rel | Sec. 3.7 | C12 vs Rotondo et al. 1.38603 [%] | 0.08 | 0.0815278 | 0.005 | **PASS** |
| SALP-He4-rel | Sec. 3.7 | He4 vs Rotondo et al. 1.40906 [%] | 0.55 | 0.547883 | 0.005 | **PASS** |
| SALP-rhoC12 | Sec. 3.7 | C12 critical density, Salpeter-GR [1e13 kg/m^3] (vs 2.12e13 of Rotondo et al.) | 2.39 | 2.386 | 0.005 | **PASS** |
| FLAT-p25 | Sec. 3.7 | relative mass deficit at +25 % in density (ideal gas) | 0.00014 | 0.00013901 | 5e-06 | **PASS** |
| FLAT-m25 | Sec. 3.7 | relative mass deficit at -25 % | 0.00024 | 0.000241239 | 5e-06 | **PASS** |
| FLAT-p50 | Sec. 3.7 | relative mass deficit at +50 % | 0.00045 | 0.000452735 | 5e-06 | **PASS** |
| FR-ratio-lo | Sec. 3.8 | omega^2_virial/(E/I) at Gamma1 = 1.25 | 5.8 | 5.7969 | 0.005 | **PASS** |
| FR-ratio-hi | Sec. 3.8 | omega^2_virial/(E/I) at Gamma1 = 5/3 | 9.19 | 9.1884 | 0.005 | **PASS** |
| FR-zero | Sec. 3.8 | both E and the virial pulsation constant vanish at exactly Gamma1 = 4/3 (|value| < 1e-12) |  | Gamma1=1.3333  omega2_virial=+0.00000e+00  E/I=+0.00000e+00 |  | **PASS** |
| FR-sign | Sec. 3.8 | both negative for Gamma1 < 4/3 and positive above, in agreement at every value tested |  | 10 values |  | **PASS** |
| FR-fig | Sec. 3.8, Fig. 7 | left panel: each curve normalised to its own maximum magnitude (the normalisation of the original figure is documented; see README) |  | script divides each curve by its own max |value|, as stated in its docstring and in the figure caption |  | **PASS** |
| HOM-kN | Sec. 5.2 | Newtonian-only linearisation k [s^-2] | -120.8 | -120.796 | 0.05 | **PASS** |
| HOM-kF | Sec. 5.2 | uncalibrated frozen-background k [s^-2] | 21 | 21.025 | 0.05 | **PASS** |
| HOM-gap | Sec. 5.2 | gap of the uncalibrated k to 27.85 [%] | 25 | 24.5063 | 0.5 | **PASS** |
| HOM-fcal | Sec. 5.2 | calibration factor f_cal | 1.048 | 1.0481 | 0.0005 | **PASS** |
| HOM-fcal-ddot | Sec. 5.2 | equilibrium check lambda_ddot(1) [s^-2] | 0 | 0 | 1e-06 | **PASS** |
| HOM-dkdf | Sec. 5.2 | dk/df_cal [s^-2] (from k(f=0) and k(f=1) if linear) | 142 | 141.821 | 0.5 | **PASS** |
| HOM-xiratio | Sec. 5.1 | xi_GR(r)/r centre/surface | 1.133 | 1.13343 | 0.0005 | **PASS** |
| HOM-lamcrit | Sec. 5.3 | lambda_crit = max 2GM_r/(r c^2) | 0.00893 | 0.00892998 | 5e-06 | **PASS** |
| HOM-lamcrit-r | Sec. 5.3 | location of lambda_crit [r/R] | 0.42 | 0.419542 | 0.005 | **PASS** |
| HOM-surfC | Sec. 5.3 | surface compactness | 0.00512 | 0.00511992 | 5e-06 | **PASS** |
| HOM-Tdyn | Sec. 5.3 | T_dyn = sqrt(R^3/GM) [s] | 0.0538 | 0.0537843 | 5e-05 | **PASS** |
| HOM-I | Sec. 5.6 | model moment of inertia at the fold [1e41 kg m^2] | 4.2 | 4.2403 | 0.05 | **PASS** |
| HOM-equip | Sec. 5.6 | equipartition sqrt(kT/I) at 1e8 K [1e-29 s^-1] | 6 | 5.70615 | 0.5 | **PASS** |
| HOM-equip-seed | Sec. 5.6 | delta lambda over a time beta (1e-7 Msun/yr): upper bound 1e-26 | 1e-26 | 6.14496e-27 | 1e-26 | **PASS** |
| EXE-k | Sec. 5.2 | exact-energy model: -sigma^2 [s^-2] | 27.857 | 27.8574 | 0.0005 | **PASS** |
| EXE-k-RR | Sec. 5.2 | Rayleigh-Ritz n=1: -sigma^2 | 27.852 | 27.8524 | 0.0005 | **PASS** |
| EXE-k-rel | Sec. 5.2 | agreement exact-energy / RR [%] | 0.02 | 0.0178273 | 0.005 | **PASS** |
| EXE-t | Sec. 5.2 | exact-energy model: time to lambda = 0.5 [s] | 2.484 | 2.48368 | 0.0005 | **PASS** |
| EXE-dt | Sec. 5.2 | exact-energy plunge 0.5 -> 0.15 [s] | 0.02713 | 0.0271332 | 5e-06 | **PASS** |
| FS-t-cal | Sec. 5.2 | calibrated model: time to lambda = 0.5 from the same seed [s] | 2.498 | 2.49789 | 0.0005 | **PASS** |
| FS-dt-cal | Sec. 5.2, 5.3 | calibrated plunge 0.5 -> 0.15 [s] | 0.02776 | 0.0277623 | 5e-06 | **PASS** |
| FS-dt-rel | Sec. 5.2 | exact-energy vs calibrated: time to 0.5 [%] | 0.6 | 0.568746 | 0.05 | **PASS** |
| FS-dt-rel2 | Sec. 5.2 | exact-energy vs calibrated: plunge [%] | 2.3 | 2.26607 | 0.05 | **PASS** |
| FS-fcal-n5 | Sec. 5.2 | calibration against 28.55 instead of 27.85: change of f_cal [%] | 0.47 | 0.47159 | 0.005 | **PASS** |
| FS-t-n5 | Sec. 5.2 | same: change of the time to lambda = 0.5 [%] | 1.2 | 1.19509 | 0.05 | **PASS** |
| FS-dt-n5 | Sec. 5.2 | same: change of the plunge duration [%] | 0.5 | 0.488236 | 0.05 | **PASS** |
| FS-sigma | Sec. 5.3 | calibrated linear rate sigma = sqrt(27.85) [1/s] | 5.278 | 5.27754 | 0.0006 | **PASS** |
| FS-efold | Sec. 5.3 | 1/sigma [s] | 0.189 | 0.189482 | 0.0006 | **PASS** |
| FS-t05-0.001 | Sec. 5.3 | time to lambda = 0.5 for delta_0 = 0.001 [s] | 1.189 | 1.18873 | 0.0005 | **PASS** |
| FS-t05-1e-06 | Sec. 5.3 | time to lambda = 0.5 for delta_0 = 1e-06 [s] | 2.498 | 2.49789 | 0.0005 | **PASS** |
| FS-t05-1e-09 | Sec. 5.3 | time to lambda = 0.5 for delta_0 = 1e-09 [s] | 3.807 | 3.80679 | 0.0005 | **PASS** |
| FS-t05-1e-12 | Sec. 5.3 | time to lambda = 0.5 for delta_0 = 1e-12 [s] | 5.117 | 5.11675 | 0.0005 | **PASS** |
| FS-inc0 | Sec. 5.3 | increase per factor 1e3 in 1/delta_0 [s] (0.001 -> 1e-06) | 1.309 | 1.30916 | 0.0006 | **PASS** |
| FS-inc1 | Sec. 5.3 | increase per factor 1e3 in 1/delta_0 [s] (1e-06 -> 1e-09) | 1.309 | 1.3089 | 0.0006 | **PASS** |
| FS-inc2 | Sec. 5.3 | increase per factor 1e3 in 1/delta_0 [s] (1e-09 -> 1e-12) | 1.31 | 1.30996 | 0.0006 | **PASS** |
| FS-ln | Sec. 5.3 | ln(1e3)/sigma [s] | 1.309 | 1.3089 | 0.0006 | **PASS** |
| FS-plunge-indep | Sec. 5.3 | plunge independent of delta_0: max relative spread over 1e-3 ... 1e-12 (upper bound; "six significant figures" = spread below 1e-6) | 1e-06 | 1.52405e-07 | 5e-07 | **PASS** |
| FS-plunge-end | Sec. 5.3 | plunge 0.5 -> end of integration [s] | 0.0295 | 0.0295011 | 5e-05 | **PASS** |
| FS-lamcrit | Sec. 5.3 | lambda_crit | 0.00893 | 0.00892998 | 5e-06 | **PASS** |
| FS-lamend | Sec. 5.3, Fig. 8 | lambda at the end of integration (stop above 1.02 lambda_crit) | 0.0091 | 0.00910858 | 5e-05 | **PASS** |
| FS-rhoamp | Sec. 5.3, Fig. 8 | central density amplification at the end [1e6] | 1.3 | 1.32327 | 0.05 | **PASS** |
| FS-plunge-Tdyn | Sec. 5.3 | plunge / T_dyn (about one half) | 0.5 | 0.516179 | 0.03 | **PASS** |
| FS-M | Sec. 2.3 | M of the fiducial star [Msun] | 1.4143 | 1.41435 | 5e-05 | **PASS** |
| FS-C | Sec. 2.3 | compactness of the fiducial star | 0.00512 | 0.00511992 | 5e-06 | **PASS** |
| FIG8-c | Sec. 5.3, Fig. 8 | offset c = 0.63 of t(0.5) = [ln(1/delta_0) - c]/sigma (mean of the four seeds) | 0.63 | 0.631741 | 0.005 | **PASS** |
| FIG8-amp | Fig. 8 | central-density amplification ~1.3e6 at lambda = 9.1e-3 | 1.3 | 1.32327 | 0.05 | **PASS** |
| IC-sigma | Sec. 5.3 | integrate_collapse: sigma [1/s] | 5.278 | 5.2775 | 0.0006 | **PASS** |
| IC-1mC | Sec. 5.3 | background 1 - 2GM/(Rc^2) | 0.99107 | 0.99107 | 6e-06 | **PASS** |
| TM-s2 | Sec. 5.4 | two-mode linearisation: -sigma^2 [s^-2] | 28.34 | 28.3366 | 0.005 | **PASS** |
| TM-RR2 | Sec. 5.4 | Rayleigh-Ritz n = 2: -sigma^2 | 28.52 | 28.525 | 0.005 | **PASS** |
| TM-diff | Sec. 5.4 | difference of the two [%] | 0.7 | 0.660416 | 0.05 | **PASS** |
| TM-1s | Sec. 5.4 | two-mode 1/sigma [s] | 0.1879 | 0.187856 | 6e-05 | **PASS** |
| TM-1s-RR | Sec. 5.4 | RR n=2 1/sigma [s] | 0.1872 | 0.187235 | 6e-05 | **PASS** |
| TM-mulam | Sec. 5.4 | eigenvector mu/lambda | 0.087 | 0.0869128 | 0.0005 | **PASS** |
| TM-second | Sec. 5.4 | second eigenvalue: frequency^2 [s^-2] | 3617 | 3617.07 | 0.5 | **PASS** |
| TM-t1e-3 | Sec. 5.4 | time to lambda=0.5, delta_0 = 1e-3 [s] | 1.171 | 1.17114 | 0.0005 | **PASS** |
| TM-t1e-6 | Sec. 5.4 | time to lambda=0.5, delta_0 = 1e-6 [s] | 2.469 | 2.46908 | 0.0005 | **PASS** |
| TM-shift1 | Sec. 5.4 | shift w.r.t. single mode, delta_0 = 1e-3 [%] | 1.5 | 1.48031 | 0.05 | **PASS** |
| TM-shift2 | Sec. 5.4 | shift w.r.t. single mode, delta_0 = 1e-6 [%] | 1.2 | 1.15357 | 0.05 | **PASS** |
| TM-mulam05 | Sec. 5.4 | mu/lambda at lambda = 0.5 | -0.076 | -0.075705 | 0.0005 | **PASS** |
| TM-mu-evol | Sec. 5.4 | mu/lambda from +0.087 (start) to -0.076 at lambda = 0.5 | 0.087 | 0.0869128 | 0.0005 | **PASS** |
| TM-plunge | Sec. 5.4 | two-mode time from lambda=0.5 to 0.173 [s] | 0.0247 | 0.0247305 | 5e-05 | **PASS** |
| TM-plunge1 | Sec. 5.4 | single-mode time from lambda = 0.5 to 0.173 [s] | 0.0271 | 0.0271239 | 5e-05 | **PASS** |
| TM-plunge-rel | Sec. 5.4 | two-mode plunge shorter by [%] | 9 | 8.8238 | 0.5 | **PASS** |
| TM-lamend | Sec. 5.4 | two-mode integration ends at lambda | 0.173 | 0.173016 | 0.0005 | **PASS** |
| TM-rhoc | Sec. 5.4 | rho_c/rho_0 at lambda = 0.29 | 79 | 78.61 | 0.6 | **PASS** |
| TM-rhohom | Sec. 5.4 | homologous rho_c/rho_0 = lambda^-3 at lambda = 0.29 | 41 | 40.7992 | 0.5 | **PASS** |
| SV-idea-Mmax | Table 2 | ideal-C(Z=0,mue=2): M_max [Msun] | 1.4166 | 1.41661 | 5e-05 | **PASS** |
| SV-idea-1.10-rho | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: rho_c [1e10 g/cm^3] | 2.54 | 2.54157 | 0.005 | **PASS** |
| SV-idea-1.10-s2 | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: -sigma^2 (n=5) [s^-2] | 1.94 | 1.93572 | 0.005 | **PASS** |
| SV-idea-1.10-fcal | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: f_cal | 1.0523 | 1.05225 | 5e-05 | **PASS** |
| SV-idea-1.10-esig | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: 1/sigma (n=5) [s] | 0.719 | 0.718752 | 0.0005 | **PASS** |
| SV-idea-1.10-dtpl | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: plunge 0.5 -> 0.15 [s] | 0.049 | 0.0489908 | 5e-05 | **PASS** |
| SV-idea-1.10-dtT | Table 2 | ideal-C(Z=0,mue=2), offset 1.10: plunge/T_dyn | 0.67 | 0.672232 | 0.005 | **PASS** |
| SV-idea-1.50-rho | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: rho_c [1e10 g/cm^3] | 3.47 | 3.46578 | 0.005 | **PASS** |
| SV-idea-1.50-s2 | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: -sigma^2 (n=5) [s^-2] | 10.8 | 10.8006 | 0.005 | **PASS** |
| SV-idea-1.50-fcal | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: f_cal | 1.0503 | 1.0503 | 5e-05 | **PASS** |
| SV-idea-1.50-esig | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: 1/sigma (n=5) [s] | 0.304 | 0.304282 | 0.0005 | **PASS** |
| SV-idea-1.50-dtpl | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: plunge 0.5 -> 0.15 [s] | 0.0371 | 0.0370517 | 5e-05 | **PASS** |
| SV-idea-1.50-dtT | Table 2 | ideal-C(Z=0,mue=2), offset 1.50: plunge/T_dyn | 0.58 | 0.584013 | 0.005 | **PASS** |
| SV-idea-2.00-rho | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: rho_c [1e10 g/cm^3] | 4.62 | 4.62104 | 0.005 | **PASS** |
| SV-idea-2.00-s2 | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: -sigma^2 (n=5) [s^-2] | 23.87 | 23.872 | 0.005 | **PASS** |
| SV-idea-2.00-fcal | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: f_cal | 1.0486 | 1.04859 | 5e-05 | **PASS** |
| SV-idea-2.00-esig | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: 1/sigma (n=5) [s] | 0.205 | 0.204671 | 0.0005 | **PASS** |
| SV-idea-2.00-dtpl | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: plunge 0.5 -> 0.15 [s] | 0.0295 | 0.0294674 | 5e-05 | **PASS** |
| SV-idea-2.00-dtT | Table 2 | ideal-C(Z=0,mue=2), offset 2.00: plunge/T_dyn | 0.53 | 0.5287 | 0.005 | **PASS** |
| SV-idea-2.16-rho | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: rho_c [1e10 g/cm^3] | 5 | 5 | 0.005 | **PASS** |
| SV-idea-2.16-s2 | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: -sigma^2 (n=5) [s^-2] | 28.55 | 28.5534 | 0.005 | **PASS** |
| SV-idea-2.16-fcal | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: f_cal | 1.0481 | 1.04814 | 5e-05 | **PASS** |
| SV-idea-2.16-esig | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: 1/sigma (n=5) [s] | 0.187 | 0.187142 | 0.0005 | **PASS** |
| SV-idea-2.16-dtpl | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: plunge 0.5 -> 0.15 [s] | 0.0277 | 0.0276714 | 5e-05 | **PASS** |
| SV-idea-2.16-dtT | Table 2 | ideal-C(Z=0,mue=2), offset 2.16: plunge/T_dyn | 0.51 | 0.514488 | 0.005 | **PASS** |
| SV-idea-3.00-rho | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: rho_c [1e10 g/cm^3] | 6.93 | 6.93156 | 0.005 | **PASS** |
| SV-idea-3.00-s2 | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: -sigma^2 (n=5) [s^-2] | 54.87 | 54.8703 | 0.005 | **PASS** |
| SV-idea-3.00-fcal | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: f_cal | 1.0463 | 1.0463 | 5e-05 | **PASS** |
| SV-idea-3.00-esig | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: 1/sigma (n=5) [s] | 0.135 | 0.134999 | 0.0005 | **PASS** |
| SV-idea-3.00-dtpl | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: plunge 0.5 -> 0.15 [s] | 0.0215 | 0.0214895 | 5e-05 | **PASS** |
| SV-idea-3.00-dtT | Table 2 | ideal-C(Z=0,mue=2), offset 3.00: plunge/T_dyn | 0.46 | 0.463447 | 0.005 | **PASS** |
| SV-He4-Mmax | Table 2 | He4: M_max [Msun] | 1.4013 | 1.40134 | 5e-05 | **PASS** |
| SV-He4-1.10-rho | Table 2 | He4, offset 1.10: rho_c [1e10 g/cm^3] | 2.58 | 2.58076 | 0.005 | **PASS** |
| SV-He4-1.10-s2 | Table 2 | He4, offset 1.10: -sigma^2 (n=5) [s^-2] | 1.96 | 1.96199 | 0.005 | **PASS** |
| SV-He4-1.10-fcal | Table 2 | He4, offset 1.10: f_cal | 1.0523 | 1.05226 | 5e-05 | **PASS** |
| SV-He4-1.10-esig | Table 2 | He4, offset 1.10: 1/sigma (n=5) [s] | 0.714 | 0.713923 | 0.0005 | **PASS** |
| SV-He4-1.10-dtpl | Table 2 | He4, offset 1.10: plunge 0.5 -> 0.15 [s] | 0.0488 | 0.0487721 | 5e-05 | **PASS** |
| SV-He4-1.10-dtT | Table 2 | He4, offset 1.10: plunge/T_dyn | 0.67 | 0.674614 | 0.005 | **PASS** |
| SV-He4-1.50-rho | Table 2 | He4, offset 1.50: rho_c [1e10 g/cm^3] | 3.52 | 3.51922 | 0.005 | **PASS** |
| SV-He4-1.50-s2 | Table 2 | He4, offset 1.50: -sigma^2 (n=5) [s^-2] | 10.95 | 10.9459 | 0.005 | **PASS** |
| SV-He4-1.50-fcal | Table 2 | He4, offset 1.50: f_cal | 1.0503 | 1.05031 | 5e-05 | **PASS** |
| SV-He4-1.50-esig | Table 2 | He4, offset 1.50: 1/sigma (n=5) [s] | 0.302 | 0.302255 | 0.0005 | **PASS** |
| SV-He4-1.50-dtpl | Table 2 | He4, offset 1.50: plunge 0.5 -> 0.15 [s] | 0.037 | 0.0369709 | 5e-05 | **PASS** |
| SV-He4-1.50-dtT | Table 2 | He4, offset 1.50: plunge/T_dyn | 0.59 | 0.587401 | 0.005 | **PASS** |
| SV-He4-2.00-rho | Table 2 | He4, offset 2.00: rho_c [1e10 g/cm^3] | 4.69 | 4.69229 | 0.005 | **PASS** |
| SV-He4-2.00-s2 | Table 2 | He4, offset 2.00: -sigma^2 (n=5) [s^-2] | 24.19 | 24.1922 | 0.005 | **PASS** |
| SV-He4-2.00-fcal | Table 2 | He4, offset 2.00: f_cal | 1.0486 | 1.0486 | 5e-05 | **PASS** |
| SV-He4-2.00-esig | Table 2 | He4, offset 2.00: 1/sigma (n=5) [s] | 0.203 | 0.203312 | 0.0005 | **PASS** |
| SV-He4-2.00-dtpl | Table 2 | He4, offset 2.00: plunge 0.5 -> 0.15 [s] | 0.0292 | 0.0292454 | 5e-05 | **PASS** |
| SV-He4-2.00-dtT | Table 2 | He4, offset 2.00: plunge/T_dyn | 0.53 | 0.528896 | 0.005 | **PASS** |
| SV-He4-2.13-rho | Table 2 | He4, offset 2.13: rho_c [1e10 g/cm^3] | 5 | 5 | 0.005 | **PASS** |
| SV-He4-2.13-s2 | Table 2 | He4, offset 2.13: -sigma^2 (n=5) [s^-2] | 27.97 | 27.9714 | 0.005 | **PASS** |
| SV-He4-2.13-fcal | Table 2 | He4, offset 2.13: f_cal | 1.0482 | 1.04824 | 5e-05 | **PASS** |
| SV-He4-2.13-esig | Table 2 | He4, offset 2.13: 1/sigma (n=5) [s] | 0.189 | 0.189079 | 0.0005 | **PASS** |
| SV-He4-2.13-dtpl | Table 2 | He4, offset 2.13: plunge 0.5 -> 0.15 [s] | 0.0276 | 0.0275906 | 5e-05 | **PASS** |
| SV-He4-2.13-dtT | Table 2 | He4, offset 2.13: plunge/T_dyn | 0.51 | 0.513502 | 0.005 | **PASS** |
| SV-He4-3.00-rho | Table 2 | He4, offset 3.00: rho_c [1e10 g/cm^3] | 7.04 | 7.03843 | 0.005 | **PASS** |
| SV-He4-3.00-s2 | Table 2 | He4, offset 3.00: -sigma^2 (n=5) [s^-2] | 55.6 | 55.604 | 0.005 | **PASS** |
| SV-He4-3.00-fcal | Table 2 | He4, offset 3.00: f_cal | 1.0463 | 1.04632 | 5e-05 | **PASS** |
| SV-He4-3.00-esig | Table 2 | He4, offset 3.00: 1/sigma (n=5) [s] | 0.134 | 0.134106 | 0.0005 | **PASS** |
| SV-He4-3.00-dtpl | Table 2 | He4, offset 3.00: plunge 0.5 -> 0.15 [s] | 0.0213 | 0.0213118 | 5e-05 | **PASS** |
| SV-He4-3.00-dtT | Table 2 | He4, offset 3.00: plunge/T_dyn | 0.46 | 0.463256 | 0.005 | **PASS** |
| SV-C12-Mmax | Table 2 | C12: M_max [Msun] | 1.3849 | 1.3849 | 5e-05 | **PASS** |
| SV-C12-1.10-rho | Table 2 | C12, offset 1.10: rho_c [1e10 g/cm^3] | 2.62 | 2.62408 | 0.005 | **PASS** |
| SV-C12-1.10-s2 | Table 2 | C12, offset 1.10: -sigma^2 (n=5) [s^-2] | 1.99 | 1.99096 | 0.005 | **PASS** |
| SV-C12-1.10-fcal | Table 2 | C12, offset 1.10: f_cal | 1.0523 | 1.05226 | 5e-05 | **PASS** |
| SV-C12-1.10-esig | Table 2 | C12, offset 1.10: 1/sigma (n=5) [s] | 0.709 | 0.70871 | 0.0005 | **PASS** |
| SV-C12-1.10-dtpl | Table 2 | C12, offset 1.10: plunge 0.5 -> 0.15 [s] | 0.0485 | 0.0484574 | 5e-05 | **PASS** |
| SV-C12-1.10-dtT | Table 2 | C12, offset 1.10: plunge/T_dyn | 0.68 | 0.676126 | 0.005 | **PASS** |
| SV-C12-1.50-rho | Table 2 | C12, offset 1.50: rho_c [1e10 g/cm^3] | 3.58 | 3.57829 | 0.005 | **PASS** |
| SV-C12-1.50-s2 | Table 2 | C12, offset 1.50: -sigma^2 (n=5) [s^-2] | 11.11 | 11.1062 | 0.005 | **PASS** |
| SV-C12-1.50-fcal | Table 2 | C12, offset 1.50: f_cal | 1.0503 | 1.05031 | 5e-05 | **PASS** |
| SV-C12-1.50-esig | Table 2 | C12, offset 1.50: 1/sigma (n=5) [s] | 0.3 | 0.300066 | 0.0005 | **PASS** |
| SV-C12-1.50-dtpl | Table 2 | C12, offset 1.50: plunge 0.5 -> 0.15 [s] | 0.0368 | 0.0367959 | 5e-05 | **PASS** |
| SV-C12-1.50-dtT | Table 2 | C12, offset 1.50: plunge/T_dyn | 0.59 | 0.589711 | 0.005 | **PASS** |
| SV-C12-2.00-rho | Table 2 | C12, offset 2.00: rho_c [1e10 g/cm^3] | 4.77 | 4.77105 | 0.005 | **PASS** |
| SV-C12-2.00-s2 | Table 2 | C12, offset 2.00: -sigma^2 (n=5) [s^-2] | 24.55 | 24.5453 | 0.005 | **PASS** |
| SV-C12-2.00-fcal | Table 2 | C12, offset 2.00: f_cal | 1.0486 | 1.04861 | 5e-05 | **PASS** |
| SV-C12-2.00-esig | Table 2 | C12, offset 2.00: 1/sigma (n=5) [s] | 0.202 | 0.201844 | 0.0005 | **PASS** |
| SV-C12-2.00-dtpl | Table 2 | C12, offset 2.00: plunge 0.5 -> 0.15 [s] | 0.029 | 0.0290257 | 5e-05 | **PASS** |
| SV-C12-2.00-dtT | Table 2 | C12, offset 2.00: plunge/T_dyn | 0.53 | 0.529473 | 0.005 | **PASS** |
| SV-C12-2.10-rho | Table 2 | C12, offset 2.10: rho_c [1e10 g/cm^3] | 5 | 5 | 0.005 | **PASS** |
| SV-C12-2.10-s2 | Table 2 | C12, offset 2.10: -sigma^2 (n=5) [s^-2] | 27.34 | 27.3395 | 0.005 | **PASS** |
| SV-C12-2.10-fcal | Table 2 | C12, offset 2.10: f_cal | 1.0483 | 1.04835 | 5e-05 | **PASS** |
| SV-C12-2.10-esig | Table 2 | C12, offset 2.10: 1/sigma (n=5) [s] | 0.191 | 0.191251 | 0.0005 | **PASS** |
| SV-C12-2.10-dtpl | Table 2 | C12, offset 2.10: plunge 0.5 -> 0.15 [s] | 0.0278 | 0.0277706 | 5e-05 | **PASS** |
| SV-C12-2.10-dtT | Table 2 | C12, offset 2.10: plunge/T_dyn | 0.52 | 0.51742 | 0.005 | **PASS** |
| SV-C12-3.00-rho | Table 2 | C12, offset 3.00: rho_c [1e10 g/cm^3] | 7.16 | 7.15658 | 0.005 | **PASS** |
| SV-C12-3.00-s2 | Table 2 | C12, offset 3.00: -sigma^2 (n=5) [s^-2] | 56.41 | 56.413 | 0.005 | **PASS** |
| SV-C12-3.00-fcal | Table 2 | C12, offset 3.00: f_cal | 1.0463 | 1.04633 | 5e-05 | **PASS** |
| SV-C12-3.00-esig | Table 2 | C12, offset 3.00: 1/sigma (n=5) [s] | 0.133 | 0.133141 | 0.0005 | **PASS** |
| SV-C12-3.00-dtpl | Table 2 | C12, offset 3.00: plunge 0.5 -> 0.15 [s] | 0.0212 | 0.0211564 | 5e-05 | **PASS** |
| SV-C12-3.00-dtT | Table 2 | C12, offset 3.00: plunge/T_dyn | 0.46 | 0.463845 | 0.005 | **PASS** |
| SV-O16-Mmax | Table 2 | O16: M_max [Msun] | 1.3782 | 1.37823 | 5e-05 | **PASS** |
| SV-O16-2.00-rho | Table 2 | O16, offset 2.00: rho_c [1e10 g/cm^3] | 4.8 | 4.80367 | 0.005 | **PASS** |
| SV-O16-2.00-s2 | Table 2 | O16, offset 2.00: -sigma^2 (n=5) [s^-2] | 24.69 | 24.6913 | 0.005 | **PASS** |
| SV-O16-2.00-fcal | Table 2 | O16, offset 2.00: f_cal | 1.0486 | 1.04862 | 5e-05 | **PASS** |
| SV-O16-2.00-esig | Table 2 | O16, offset 2.00: 1/sigma (n=5) [s] | 0.201 | 0.201246 | 0.0005 | **PASS** |
| SV-O16-2.00-dtpl | Table 2 | O16, offset 2.00: plunge 0.5 -> 0.15 [s] | 0.0289 | 0.0289413 | 5e-05 | **PASS** |
| SV-O16-2.00-dtT | Table 2 | O16, offset 2.00: plunge/T_dyn | 0.53 | 0.529803 | 0.005 | **PASS** |
| SV-O16-2.08-rho | Table 2 | O16, offset 2.08: rho_c [1e10 g/cm^3] | 5 | 5 | 0.005 | **PASS** |
| SV-O16-2.08-s2 | Table 2 | O16, offset 2.08: -sigma^2 (n=5) [s^-2] | 27.08 | 27.0812 | 0.005 | **PASS** |
| SV-O16-2.08-fcal | Table 2 | O16, offset 2.08: f_cal | 1.0484 | 1.04839 | 5e-05 | **PASS** |
| SV-O16-2.08-esig | Table 2 | O16, offset 2.08: 1/sigma (n=5) [s] | 0.192 | 0.192161 | 0.0005 | **PASS** |
| SV-O16-2.08-dtpl | Table 2 | O16, offset 2.08: plunge 0.5 -> 0.15 [s] | 0.0279 | 0.027901 | 5e-05 | **PASS** |
| SV-O16-2.08-dtT | Table 2 | O16, offset 2.08: plunge/T_dyn | 0.52 | 0.520083 | 0.005 | **PASS** |
| SV-O16-3.00-rho | Table 2 | O16, offset 3.00: rho_c [1e10 g/cm^3] | 7.21 | 7.2055 | 0.005 | **PASS** |
| SV-O16-3.00-s2 | Table 2 | O16, offset 3.00: -sigma^2 (n=5) [s^-2] | 56.75 | 56.7475 | 0.005 | **PASS** |
| SV-O16-3.00-fcal | Table 2 | O16, offset 3.00: f_cal | 1.0463 | 1.04634 | 5e-05 | **PASS** |
| SV-O16-3.00-esig | Table 2 | O16, offset 3.00: 1/sigma (n=5) [s] | 0.133 | 0.132748 | 0.0005 | **PASS** |
| SV-O16-3.00-dtpl | Table 2 | O16, offset 3.00: plunge 0.5 -> 0.15 [s] | 0.0211 | 0.0210915 | 5e-05 | **PASS** |
| SV-O16-3.00-dtT | Table 2 | O16, offset 3.00: plunge/T_dyn | 0.46 | 0.464049 | 0.005 | **PASS** |
| SV-Fe56-Mmax | Table 2 | Fe56: M_max [Msun] | 1.1507 | 1.15071 | 5e-05 | **PASS** |
| SV-Fe56-1.71-rho | Table 2 | Fe56, offset 1.71: rho_c [1e10 g/cm^3] | 5 | 5 | 0.005 | **PASS** |
| SV-Fe56-1.71-s2 | Table 2 | Fe56, offset 1.71: -sigma^2 (n=5) [s^-2] | 19.07 | 19.0731 | 0.005 | **PASS** |
| SV-Fe56-1.71-fcal | Table 2 | Fe56, offset 1.71: f_cal | 1.0497 | 1.0497 | 5e-05 | **PASS** |
| SV-Fe56-1.71-esig | Table 2 | Fe56, offset 1.71: 1/sigma (n=5) [s] | 0.229 | 0.228975 | 0.0005 | **PASS** |
| SV-Fe56-1.71-dtpl | Table 2 | Fe56, offset 1.71: plunge 0.5 -> 0.15 [s] | 0.0305 | 0.0304581 | 5e-05 | **PASS** |
| SV-Fe56-1.71-dtT | Table 2 | Fe56, offset 1.71: plunge/T_dyn | 0.57 | 0.571929 | 0.005 | **PASS** |
| SV-Fe56-2.00-rho | Table 2 | Fe56, offset 2.00: rho_c [1e10 g/cm^3] | 5.86 | 5.8642 | 0.005 | **PASS** |
| SV-Fe56-2.00-s2 | Table 2 | Fe56, offset 2.00: -sigma^2 (n=5) [s^-2] | 28.58 | 28.5764 | 0.005 | **PASS** |
| SV-Fe56-2.00-fcal | Table 2 | Fe56, offset 2.00: f_cal | 1.0488 | 1.04882 | 5e-05 | **PASS** |
| SV-Fe56-2.00-esig | Table 2 | Fe56, offset 2.00: 1/sigma (n=5) [s] | 0.187 | 0.187067 | 0.0005 | **PASS** |
| SV-Fe56-2.00-dtpl | Table 2 | Fe56, offset 2.00: plunge 0.5 -> 0.15 [s] | 0.0269 | 0.0268878 | 5e-05 | **PASS** |
| SV-Fe56-2.00-dtT | Table 2 | Fe56, offset 2.00: plunge/T_dyn | 0.54 | 0.542556 | 0.005 | **PASS** |
| SV-Fe56-3.00-rho | Table 2 | Fe56, offset 3.00: rho_c [1e10 g/cm^3] | 8.8 | 8.7963 | 0.005 | **PASS** |
| SV-Fe56-3.00-s2 | Table 2 | Fe56, offset 3.00: -sigma^2 (n=5) [s^-2] | 65.67 | 65.6668 | 0.005 | **PASS** |
| SV-Fe56-3.00-fcal | Table 2 | Fe56, offset 3.00: f_cal | 1.0467 | 1.04665 | 5e-05 | **PASS** |
| SV-Fe56-3.00-esig | Table 2 | Fe56, offset 3.00: 1/sigma (n=5) [s] | 0.123 | 0.123403 | 0.0005 | **PASS** |
| SV-Fe56-3.00-dtpl | Table 2 | Fe56, offset 3.00: plunge 0.5 -> 0.15 [s] | 0.0197 | 0.0196775 | 5e-05 | **PASS** |
| SV-Fe56-3.00-dtT | Table 2 | Fe56, offset 3.00: plunge/T_dyn | 0.48 | 0.47736 | 0.005 | **PASS** |
| SV-i-C12 | Sec. 5.5 (i) | fiducial density: -sigma^2 ideal gas | 28.55 | 28.5534 | 0.005 | **PASS** |
| SV-ii | Sec. 5.5 (ii) | f_cal lies between 1.046 and 1.052 across all equations of state and offsets (0.6 % spread) |  | range 1.0463-1.0523 |  | **PASS** |
| SV-iii-a | Sec. 5.5 (iii) | -sigma^2 grows by ~5 between offsets 1.5 and 3 (ideal gas) | 5 | 5.08031 | 0.5 | **PASS** |
| SV-iii-b | Sec. 5.5 (iii) | ideal gas, He4, C12, O16 agree at offset 2.0 to [%] | 3.5 | 3.43222 | 0.1 | **PASS** |
| SV-iii-c | Sec. 5.5 | Fe56 at offset 2.0: -sigma^2 above the ideal gas by [%] | 20 | 19.7069 | 0.5 | **PASS** |
| SV-iii-d | Sec. 5.5 | Fe56 at the fiducial density: -sigma^2 below the ideal gas by [%] | 33 | 33.2018 | 0.5 | **PASS** |
| SV-iv-lo | Sec. 5.5 (iv) | minimum of plunge/T_dyn over the survey | 0.46 | 0.463256 | 0.005 | **PASS** |
| SV-iv-hi | Sec. 5.5 (iv) | maximum of plunge/T_dyn over the survey | 0.69 | 0.69339 | 0.005 | **PASS** |
| SV-lat | Sec. 5.5 (i) | lattice decrease of -sigma^2 at fixed rho_c: He4, C12, O16 relative to the ideal gas [%] (min of 3) | 2 | 2.03811 | 0.1 | **PASS** |
| SV-lat2 | Sec. 5.5 (i) | ... maximum of the 3 [%] | 5 | 5.15565 | 0.5 | **PASS** |
| SV-spot | Sec. 5.5 | two rows of the survey recomputed now (survey_spot.py) agree with survey4.json to 1e-4 relative |  | max rel diff [1.345434648269972e-05, 3.791384769580297e-06] |  | **PASS** |
| FOLD-rhomax | Sec. 5.6 | rho_max [1e13 kg/m^3] | 2.3106 | 2.3106 | 5e-05 | **PASS** |
| FOLD-Mmax | Sec. 5.6 | M_max [Msun] | 1.41661 | 1.41661 | 5e-06 | **PASS** |
| FOLD-kappa | Sec. 5.6 | kappa [s^-2] (y = 0.002 ... 0.02) | 18.733 | 18.7334 | 0.0005 | **PASS** |
| FOLD-kappa-const | Sec. 5.6 | kappa constant over y = +-0.002 ... 0.02 to (upper bound) [%] | 0.01 | 0.00483577 | 0.005 | **PASS** |
| FOLD-a | Sec. 5.6 | a [1e-3 Msun] | 4.027 | 4.02681 | 0.0005 | **PASS** |
| FOLD-fcal | Sec. 5.6 | f_cal at the fold (k = 0) | 1.059 | 1.05866 | 0.0005 | **PASS** |
| FOLD-quad | Sec. 5.6 | model quadratic coefficient [s^-2] | 26.09 | 26.0865 | 0.005 | **PASS** |
| FOLD-3k2 | Sec. 5.6 | exact-GR 3 kappa/2 [s^-2] | 28.1 | 28.0995 | 0.005 | **PASS** |
| FOLD-quad-rel | Sec. 5.6 | model coefficient below exact-GR [%] | 7.2 | 7.16393 | 0.05 | **PASS** |
| FOLD-rescale | Sec. 5.6 | force rescaling to the exact-GR coefficient | 1.077 | 1.07717 | 0.0005 | **PASS** |
| FT-1e-09-beta | Table 3 | Mdot = 1e-09: beta [s] | 270.5 | 270.512 | 0.05 | **PASS** |
| FT-1e-09-alpha | Table 3 | Mdot = 1e-09: alpha | 1.5e-06 | 1.45898e-06 | 5e-08 | **PASS** |
| FT-1e-09-t | Table 3 | Mdot = 1e-09: t_coll (rescaled model, to lambda = 0.5) [s] | 922 | 921.902 | 0.5 | **PASS** |
| FT-1e-09-raw | Table 3 | Mdot = 1e-09: raw model [s] | 950 | 949.763 | 0.6 | **PASS** |
| FT-1e-09-dM | Table 3 | Mdot = 1e-09: Delta M [Msun] | 2.9e-14 | 2.92133e-14 | 1.16e-15 | **PASS** |
| FT-1e-09-pl | Table 3 | Mdot = 1e-09: plunge [s] | 0.052 | 0.0518233 | 0.0005 | **PASS** |
| FT-1e-08-beta | Table 3 | Mdot = 1e-08: beta [s] | 170.7 | 170.681 | 0.05 | **PASS** |
| FT-1e-08-alpha | Table 3 | Mdot = 1e-08: alpha | 3.7e-06 | 3.6648e-06 | 5e-08 | **PASS** |
| FT-1e-08-t | Table 3 | Mdot = 1e-08: t_coll (rescaled model, to lambda = 0.5) [s] | 581 | 581.314 | 0.5 | **PASS** |
| FT-1e-08-raw | Table 3 | Mdot = 1e-08: raw model [s] | 599 | 598.893 | 0.6 | **PASS** |
| FT-1e-08-dM | Table 3 | Mdot = 1e-08: Delta M [Msun] | 1.8e-13 | 1.84207e-13 | 7.2e-15 | **PASS** |
| FT-1e-08-pl | Table 3 | Mdot = 1e-08: plunge [s] | 0.052 | 0.0518233 | 0.0005 | **PASS** |
| FT-1e-07-beta | Table 3 | Mdot = 1e-07: beta [s] | 107.7 | 107.693 | 0.05 | **PASS** |
| FT-1e-07-alpha | Table 3 | Mdot = 1e-07: alpha | 9.2e-06 | 9.20556e-06 | 5e-08 | **PASS** |
| FT-1e-07-t | Table 3 | Mdot = 1e-07: t_coll (rescaled model, to lambda = 0.5) [s] | 366 | 366.417 | 0.5 | **PASS** |
| FT-1e-07-raw | Table 3 | Mdot = 1e-07: raw model [s] | 378 | 377.516 | 0.6 | **PASS** |
| FT-1e-07-dM | Table 3 | Mdot = 1e-07: Delta M [Msun] | 1.2e-12 | 1.1611e-12 | 4.8e-14 | **PASS** |
| FT-1e-07-pl | Table 3 | Mdot = 1e-07: plunge [s] | 0.052 | 0.0518233 | 0.0005 | **PASS** |
| FT-1e-06-beta | Table 3 | Mdot = 1e-06: beta [s] | 67.9 | 67.9495 | 0.05 | **PASS** |
| FT-1e-06-alpha | Table 3 | Mdot = 1e-06: alpha | 2.3e-05 | 2.31233e-05 | 5e-07 | **PASS** |
| FT-1e-06-t | Table 3 | Mdot = 1e-06: t_coll (rescaled model, to lambda = 0.5) [s] | 231 | 230.824 | 0.5 | **PASS** |
| FT-1e-06-raw | Table 3 | Mdot = 1e-06: raw model [s] | 238 | 237.847 | 0.6 | **PASS** |
| FT-1e-06-dM | Table 3 | Mdot = 1e-06: Delta M [Msun] | 7.3e-12 | 7.31437e-12 | 2.92e-13 | **PASS** |
| FT-1e-06-pl | Table 3 | Mdot = 1e-06: plunge [s] | 0.052 | 0.0518233 | 0.0005 | **PASS** |
| FT-1e-05-beta | Table 3 | Mdot = 1e-05: beta [s] | 42.9 | 42.8732 | 0.05 | **PASS** |
| FT-1e-05-alpha | Table 3 | Mdot = 1e-05: alpha | 5.8e-05 | 5.80832e-05 | 5e-07 | **PASS** |
| FT-1e-05-t | Table 3 | Mdot = 1e-05: t_coll (rescaled model, to lambda = 0.5) [s] | 145 | 145.268 | 0.5 | **PASS** |
| FT-1e-05-raw | Table 3 | Mdot = 1e-05: raw model [s] | 150 | 149.74 | 0.6 | **PASS** |
| FT-1e-05-dM | Table 3 | Mdot = 1e-05: Delta M [Msun] | 4.6e-11 | 4.60325e-11 | 1.84e-12 | **PASS** |
| FT-1e-05-pl | Table 3 | Mdot = 1e-05: plunge [s] | 0.052 | 0.0518233 | 0.0005 | **PASS** |
| FT-plunge-raw | Table 3 caption | raw plunge [s] | 0.054 | 0.0537856 | 0.0005 | **PASS** |
| FT-raw-3pc | Table 3 caption | raw model longer than rescaled by [%] | 3 | 3.02907 | 0.5 | **PASS** |
| FT-Mdot-pow | Sec. 5.6 | four decades in Mdot move the delay by a factor of [6] | 6 | 6.34623 | 0.5 | **PASS** |
| FT-overshoot | Sec. 5.6 | mass accreted up to lambda = 0.5 (Delta M of Table 3), maximum over Mdot <= [1e-11 Msun] | 4.6 | 4.60325 | 0.1 | **PASS** |
| FT-y1e2-resc | Sec. 5.6 | rescaled model reaches y = 1e-2 at s = t/beta (Mdot = 1e-7) | 3.337 | 3.33679 | 0.001 | **PASS** |
| FT-y1e2-u | Sec. 5.6 | u = y/alpha at y = 1e-2 | 1086 | 1086.3 | 1.5 | **PASS** |
| FT-start-20 | Sec. 5.6 | raw model start independence: t(0.5)/beta for s_0 = -20 | 3.5056 | 3.50564 | 0.0002 | **PASS** |
| FT-start-30 | Sec. 5.6 | raw model start independence: t(0.5)/beta for s_0 = -30 | 3.5055 | 3.50546 | 0.0002 | **PASS** |
| FT-start-45 | Sec. 5.6 | raw model start independence: t(0.5)/beta for s_0 = -45 | 3.5055 | 3.50549 | 0.0002 | **PASS** |
| FT-restart | Sec. 5.6 | a rest start shifts t by [%] | 0.6 | 0.553555 | 0.05 | **PASS** |
| FT-rate | Sec. 5.6 | accretion-rate dependence of the raw t/beta (spread) [%] | 0.6 | 0.525989 | 0.1 | **PASS** |
| FT-ratediff | Sec. 5.6 | seed shifts delay 3.30/3.46: pole estimator of the PI solution (see pi_checks) | 3.41167 | 3.41167 | 1e-05 | **PASS** |
| PI-const | Sec. 5.6 | s_* = 3.41167 constant to 2e-6 for s_0 = -40 ... -320 |  | values [3.411670444543177, 3.4116716220567964, 3.4116718135826676, 3.4116717943242243] |  | **PASS** |
| PI-lit-x | Sec. 5.6 | literature check: first real pole of Y''=6Y^2+x with the Boutroux data, x_* | 2.38417 | 2.38417 | 2e-06 | **PASS** |
| PI-lit-s | Sec. 5.6 | literature check: s_* = 6^(1/5) x_* | 3.41167 | 3.41167 | 5e-06 | **PASS** |
| PI-fold-s | Sec. 5.6 | universal solution reaches u = 1086 (y=1e-2, Mdot=1e-7) at s | 3.337 | 3.3374 | 0.0006 | **PASS** |
| PI-NF-t | Sec. 5.6 | normal form reaches lambda=0.5 at t [s] (Mdot = 1e-7) | 367 | 367.1 | 0.5 | **PASS** |
| PI-off+ | Sec. 5.6 | pole with offset A = +0.1 at s = -100 | 3.3 | 3.302 | 0.005 | **PASS** |
| PI-off- | Sec. 5.6 | pole with offset A = -0.1 | 3.46 | 3.464 | 0.005 | **PASS** |
| PI-off1 | Sec. 5.6 | pole with offset A = +1 | 2.56 | 2.565 | 0.006 | **PASS** |
| PI-decay | Sec. 5.6 | oscillation about the slow manifold decays as (-s)^(-1/8): fitted exponent | -0.125 | -0.129 | 0.01 | **PASS** |
| PI-dd-lo | Sec. 5.6 | delta rho/rho of the 3 % shift (A=0.1, alpha=1.5e-6): lower end [1e-7] | 1 | 1.459 | 0.5 | **PASS** |
| PI-dd-hi | Sec. 5.6 | ... upper end (alpha = 5.8e-5) [1e-6] | 6 | 5.80838 | 0.2 | **PASS** |
| PI-next-mode | Sec. 5.6 | next radial mode omega at the fold (N=200) [rad/s] | 37 | 37.03 | 0.5 | **PASS** |
| PI-next-ratio | Sec. 5.6 | ratio of that mode to 1/beta (Mdot = 1e-7) [1e3] | 4 | 3.98776 | 0.3 | **PASS** |
| DAMP-0-s | Table 4 | nu~ = 0: s_* = t_collapse/beta | 3.412 | 3.41167 | 0.0005 | **PASS** |
| DAMP-0.1-s | Table 4 | nu~ = 0.1: s_* = t_collapse/beta | 3.519 | 3.5193 | 0.0005 | **PASS** |
| DAMP-0.1-airy | Table 4 | nu~ = 0.1: 2.338 nu~^(2/3) | 0.504 | 0.50373 | 0.0005 | **PASS** |
| DAMP-0.3-s | Table 4 | nu~ = 0.3: s_* = t_collapse/beta | 3.732 | 3.73162 | 0.0005 | **PASS** |
| DAMP-0.3-airy | Table 4 | nu~ = 0.3: 2.338 nu~^(2/3) | 1.048 | 1.0478 | 0.0005 | **PASS** |
| DAMP-1-s | Table 4 | nu~ = 1: s_* = t_collapse/beta | 4.447 | 4.44695 | 0.0005 | **PASS** |
| DAMP-1-airy | Table 4 | nu~ = 1: 2.338 nu~^(2/3) | 2.338 | 2.33811 | 0.0005 | **PASS** |
| DAMP-1-ratio | Table 4 / text | nu~ = 1: ratio to the Airy law | 1.9 | 1.90194 | 0.005 | **PASS** |
| DAMP-3-s | Table 4 | nu~ = 3: s_* = t_collapse/beta | 6.307 | 6.30668 | 0.0005 | **PASS** |
| DAMP-3-airy | Table 4 | nu~ = 3: 2.338 nu~^(2/3) | 4.863 | 4.86346 | 0.0005 | **PASS** |
| DAMP-3-ratio | Table 4 / text | nu~ = 3: ratio to the Airy law | 1.3 | 1.29675 | 0.005 | **PASS** |
| DAMP-10-s | Table 4 | nu~ = 10: s_* = t_collapse/beta | 11.63 | 11.6279 | 0.005 | **PASS** |
| DAMP-10-airy | Table 4 | nu~ = 10: 2.338 nu~^(2/3) | 10.85 | 10.8525 | 0.005 | **PASS** |
| DAMP-10-ratio | Table 4 / text | nu~ = 10: ratio to the Airy law | 1.07 | 1.07144 | 0.005 | **PASS** |
| DAMP-airy0 | Sec. 5.7 | first zero of Ai: tau_* = -a_1 | 2.33811 | 2.33811 | 5e-06 | **PASS** |
| DAMP-ex1-t | Sec. 5.7 | Mdot = 1e-7, nu = 10: t_collapse [1e4 s] | 2.6 | 2.65 | 0.06 | **PASS** |
| DAMP-ex1-nt | Sec. 5.7 | ... nu~ | 1100 | 1077 | 30 | **PASS** |
| DAMP-ex1-dM | Sec. 5.7 | ... Delta M [1e-11 Msun] | 8 | 8.38 | 0.5 | **PASS** |
| DAMP-ex2-dM | Sec. 5.7 | Mdot = 1e-5, nu = 100: Delta M [1e-9 Msun] | 8 | 8.38 | 0.5 | **PASS** |
| VISC-coef | Sec. 5.7 | nu = 9 Gamma1bar W_P/I Delta_r tau: coefficient [1e3 s^-2] | 6.5 | 6.53996 | 0.05 | **PASS** |
| VISC-G1 | Sec. 5.7 | Gamma1bar at the fold | 1.337 | 1.33749 | 0.0006 | **PASS** |
| VISC-WP | Sec. 5.7 | W_P = int P dV at the fold [1e44 J] | 2.3 | 2.30377 | 0.005 | **PASS** |
| VISC-I | Sec. 5.7 | I at the fold [1e41 kg m^2] | 4.24 | 4.2403 | 0.005 | **PASS** |
| VISC-1e-09 | Sec. 5.7 | Delta_r tau needed for nu~ = 1 at Mdot = 1e-09 [s] | 5.7e-07 | 5.652e-07 | 5e-09 | **PASS** |
| VISC-1e-07 | Sec. 5.7 | Delta_r tau needed for nu~ = 1 at Mdot = 1e-07 [s] | 1.4e-06 | 1.42e-06 | 5e-08 | **PASS** |
| VISC-1e-05 | Sec. 5.7 | Delta_r tau needed for nu~ = 1 at Mdot = 1e-05 [s] | 3.6e-06 | 3.566e-06 | 5e-08 | **PASS** |
| VISC-ms | Sec. 5.7 | Delta_r = 1e-3 needs tau >= [ms] at Mdot = 1e-7 | 1.4 | 1.42 | 0.05 | **PASS** |
| THERM-EF | Sec. 5.7 | electron Fermi kinetic energy at rho_max [MeV] | 11.2 | 11.1555 | 0.05 | **PASS** |
| THERM-kT | Sec. 5.7 | (k_B T/E_F)^2 at 1e8 K | 6e-07 | 5.96852e-07 | 5e-08 | **PASS** |
| THERM-ion4 | Sec. 5.7 | ion thermal pressure / electron pressure, A = 4 | 0.0015 | 0.00148168 | 5e-05 | **PASS** |
| THERM-ion12 | Sec. 5.7 | ion thermal pressure / electron pressure, A = 12 | 0.00049 | 0.000493894 | 5e-06 | **PASS** |
| THERM-ion16 | Sec. 5.7 | ion thermal pressure / electron pressure, A = 16 | 0.00037 | 0.000370421 | 5e-06 | **PASS** |
| THERM-ion56 | Sec. 5.7 | ion thermal pressure / electron pressure, A = 56 | 0.00011 | 0.000105835 | 5e-06 | **PASS** |
| FC-idea-rho | Table 5 | ideal(Z=0): rho_max [1e10 g/cm^3] | 2.311 | 2.3106 | 0.0006 | **PASS** |
| FC-idea-M | Table 5 | ideal(Z=0): M_max [Msun] | 1.4166 | 1.41661 | 5e-05 | **PASS** |
| FC-idea-kappa | Table 5 | ideal(Z=0): kappa [s^-2] | 18.73 | 18.7334 | 0.005 | **PASS** |
| FC-idea-a | Table 5 | ideal(Z=0): a [1e-3 Msun] | 4.027 | 4.02689 | 0.0005 | **PASS** |
| FC-idea-beta | Table 5 | ideal(Z=0): beta (Mdot = 1e-7) [s] | 107.7 | 107.691 | 0.05 | **PASS** |
| FC-idea-t | Table 5 | ideal(Z=0): t_coll = 3.41167 beta [s] | 367 | 367.407 | 0.6 | **PASS** |
| FC-He4-rho | Table 5 | He4: rho_max [1e10 g/cm^3] | 2.346 | 2.34621 | 0.0006 | **PASS** |
| FC-He4-M | Table 5 | He4: M_max [Msun] | 1.4013 | 1.40134 | 5e-05 | **PASS** |
| FC-He4-kappa | Table 5 | He4: kappa [s^-2] | 18.99 | 18.9861 | 0.005 | **PASS** |
| FC-He4-a | Table 5 | He4: a [1e-3 Msun] | 3.976 | 3.97619 | 0.0005 | **PASS** |
| FC-He4-beta | Table 5 | He4: beta (Mdot = 1e-7) [s] | 106.8 | 106.844 | 0.05 | **PASS** |
| FC-He4-t | Table 5 | He4: t_coll = 3.41167 beta [s] | 364 | 364.518 | 0.6 | **PASS** |
| FC-C12-rho | Table 5 | C12: rho_max [1e10 g/cm^3] | 2.386 | 2.38557 | 0.0006 | **PASS** |
| FC-C12-M | Table 5 | C12: M_max [Msun] | 1.3849 | 1.3849 | 5e-05 | **PASS** |
| FC-C12-kappa | Table 5 | C12: kappa [s^-2] | 19.26 | 19.2649 | 0.005 | **PASS** |
| FC-C12-a | Table 5 | C12: a [1e-3 Msun] | 3.922 | 3.92172 | 0.0005 | **PASS** |
| FC-C12-beta | Table 5 | C12: beta (Mdot = 1e-7) [s] | 105.9 | 105.931 | 0.05 | **PASS** |
| FC-C12-t | Table 5 | C12: t_coll = 3.41167 beta [s] | 361 | 361.401 | 0.6 | **PASS** |
| FC-O16-rho | Table 5 | O16: rho_max [1e10 g/cm^3] | 2.402 | 2.40187 | 0.0006 | **PASS** |
| FC-O16-M | Table 5 | O16: M_max [Msun] | 1.3782 | 1.37823 | 5e-05 | **PASS** |
| FC-O16-kappa | Table 5 | O16: kappa [s^-2] | 19.38 | 19.3801 | 0.005 | **PASS** |
| FC-O16-a | Table 5 | O16: a [1e-3 Msun] | 3.9 | 3.89962 | 0.0005 | **PASS** |
| FC-O16-beta | Table 5 | O16: beta (Mdot = 1e-7) [s] | 105.6 | 105.559 | 0.05 | **PASS** |
| FC-O16-t | Table 5 | O16: t_coll = 3.41167 beta [s] | 360 | 360.132 | 0.6 | **PASS** |
| FC-Fe56-rho | Table 5 | Fe56: rho_max [1e10 g/cm^3] | 2.932 | 2.93211 | 0.0006 | **PASS** |
| FC-Fe56-M | Table 5 | Fe56: M_max [Msun] | 1.1507 | 1.15071 | 5e-05 | **PASS** |
| FC-Fe56-kappa | Table 5 | Fe56: kappa [s^-2] | 22.44 | 22.4387 | 0.005 | **PASS** |
| FC-Fe56-a | Table 5 | Fe56: a [1e-3 Msun] | 3.094 | 3.09378 | 0.0005 | **PASS** |
| FC-Fe56-beta | Table 5 | Fe56: beta (Mdot = 1e-7) [s] | 95 | 95.0455 | 0.05 | **PASS** |
| FC-Fe56-t | Table 5 | Fe56: t_coll = 3.41167 beta [s] | 324 | 324.264 | 0.6 | **PASS** |
| FC-dkappa | Sec. 5.9 | kappa rises by 1 ... 20 %: maximum [%] | 20 | 19.7791 | 0.5 | **PASS** |
| FC-da | Sec. 5.9 | a falls by 1 ... 23 %: maximum [%] | 23 | 23.172 | 0.5 | **PASS** |
| FC-dt | Sec. 5.9 | delay changes by 1-2 % (He4, C12, O16): maximum [%] | 2 | 1.98 | 0.1 | **PASS** |
| FC-dt-Fe | Sec. 5.9 | delay change for Fe56 [%] | 12 | 11.7425 | 0.5 | **PASS** |
| FC-dalpha | Sec. 5.9 | seed alpha changes by at most about [%] | 7 | 7.18052 | 0.5 | **PASS** |
| FC-rhomax-g | Sec. 5.6 (iii) | rho_max of the ideal gas [1e10 g/cm^3] | 2.3 | 2.3106 | 0.05 | **PASS** |
| LE-g100 | Sec. 5.8 | unstable growth rate, N = 100 [1/s] | 5.013 | 5.0133 | 0.0005 | **PASS** |
| LE-s100 | Sec. 5.8 | lowest discrete sigma^2 at the turning point, N = 100 [s^-2] | 1.585 | 1.58504 | 0.0005 | **PASS** |
| LE-w100 | Sec. 5.8 | stable-branch fundamental frequency, N = 100 [rad/s] | 2.922 | 2.92243 | 0.0005 | **PASS** |
| LE-g200 | Sec. 5.8 | unstable growth rate, N = 200 [1/s] | 5.237 | 5.2374 | 0.0005 | **PASS** |
| LE-s200 | Sec. 5.8 | lowest discrete sigma^2 at the turning point, N = 200 [s^-2] | 0.521 | 0.52142 | 0.0005 | **PASS** |
| LE-w200 | Sec. 5.8 | stable-branch fundamental frequency, N = 200 [rad/s] | 2.842 | 2.84199 | 0.0005 | **PASS** |
| LE-g400 | Sec. 5.8 | unstable growth rate, N = 400 [1/s] | 5.309 | 5.3094 | 0.0005 | **PASS** |
| LE-s400 | Sec. 5.8 | lowest discrete sigma^2 at the turning point, N = 400 [s^-2] | 0.17 | 0.16965 | 0.0005 | **PASS** |
| LE-w400 | Sec. 5.8 | stable-branch fundamental frequency, N = 400 [rad/s] | 2.815 | 2.81486 | 0.0005 | **PASS** |
| LE-g800 | Sec. 5.8 | unstable growth rate, N = 800 [1/s] | 5.333 | 5.3326 | 0.0005 | **PASS** |
| LE-s800 | Sec. 5.8 | lowest discrete sigma^2 at the turning point, N = 800 [s^-2] | 0.055 | 0.05502 | 0.0005 | **PASS** |
| LE-re | Sec. 5.8 | all Jacobian eigenvalues of the stable-branch configuration have |Re| <= 4e-10 s^-1 (N = 100 ... 800) |  | max |Re| = 4.24e-10 |  | **PASS** |
| LE-ratio | Sec. 5.8 | successive differences of the growth rate shrink by a factor (per doubling) | 3.1 | 3.1125 | 0.05 | **PASS** |
| LE-pred800 | Sec. 5.8 | N = 800 growth rate predicted from the first three levels (ratio r) | 5.3326 | 5.33253 | 0.0001 | **PASS** |
| LE-found800 | Sec. 5.8 | N = 800 growth rate found | 5.3326 | 5.3326 | 6e-05 | **PASS** |
| LE-limit | Sec. 5.8 | limit with the observed ratio (from N = 100, 200, 400) [1/s] | 5.3435 | 5.34348 | 0.0006 | **PASS** |
| LE-limit-2nd | Sec. 5.8 | limit assuming exactly second order [1/s] | 5.34 | 5.34033 | 0.001 | **PASS** |
| LE-RR | Sec. 5.8 | Rayleigh-Ritz sqrt(28.55) [1/s] | 5.343 | 5.34322 | 0.0015 | **PASS** |
| LE-tp-extrap | Sec. 5.8 | turning-point sigma^2 extrapolated to N -> oo [s^-2] (Richardson with the observed ratio) | -0.0004 | -0.000390462 | 5e-05 | **PASS** |
| LE-tp-order | Sec. 5.8 | ratio of successive turning-point sigma^2 values (about 3, i.e. 1st order) | 3.1 | 3.08342 | 0.2 | **PASS** |
| LE-RR-stable | Sec. 5.8 | stable branch (rho_c = 1e10 g/cm^3): Rayleigh-Ritz fundamental frequency [rad/s] | 2.802 | 2.8018 | 0.0005 | **PASS** |
| LE-fold-N400 | Sec. 5.8 | discrete fold at sigma^2_N = 0.17 (N = 400) <=> y = sigma^2/kappa | 0.009 | 0.00905621 | 0.0005 | **PASS** |
| LE-y01 | Sec. 5.8 | hand-off at y = 0.1: sigma^2 = -kappa y [s^-2] | -1.9 | -1.8733 | 0.05 | **PASS** |
| LAG-200-t05 | Table 6 | Lagrangian N = 200: t(lambda_h = 0.5) [s] | 1.9209 | 1.9209 | 5e-05 | **PASS** |
| LAG-200-t015 | Table 6 | Lagrangian N = 200: t(0.15) [s] | 1.9745 | 1.9745 | 5e-05 | **PASS** |
| LAG-200-plunge | Table 6 | Lagrangian N = 200: plunge [s] | 0.0536 | 0.0536535 | 0.0001 | **PASS** |
| LAG-200-tAH | Table 6 | Lagrangian N = 200: t_AH [s] | 1.9777 | 1.97771 | 5e-05 | **PASS** |
| LAG-200-mAH | Table 6 | Lagrangian N = 200: m_AH/M | 0.178 | 0.177769 | 0.0005 | **PASS** |
| LAG-400-t05 | Table 6 | Lagrangian N = 400: t(lambda_h = 0.5) [s] | 1.7622 | 1.7622 | 5e-05 | **PASS** |
| LAG-400-t015 | Table 6 | Lagrangian N = 400: t(0.15) [s] | 1.8155 | 1.8155 | 5e-05 | **PASS** |
| LAG-400-plunge | Table 6 | Lagrangian N = 400: plunge [s] | 0.0533 | 0.0533558 | 0.0001 | **PASS** |
| LAG-400-tAH | Table 6 | Lagrangian N = 400: t_AH [s] | 1.8187 | 1.81872 | 5e-05 | **PASS** |
| LAG-400-mAH | Table 6 | Lagrangian N = 400: m_AH/M | 0.175 | 0.17532 | 0.0005 | **PASS** |
| LAG-800-t05 | Table 6 | Lagrangian N = 800: t(lambda_h = 0.5) [s] | 1.72 | 1.72 | 5e-05 | **PASS** |
| LAG-800-t015 | Table 6 | Lagrangian N = 800: t(0.15) [s] | 1.7731 | 1.7731 | 5e-05 | **PASS** |
| LAG-800-plunge | Table 6 | Lagrangian N = 800: plunge [s] | 0.0532 | 0.053176 | 0.0001 | **PASS** |
| LAG-800-tAH | Table 6 | Lagrangian N = 800: t_AH [s] | 1.7763 | 1.7763 | 5e-05 | **PASS** |
| LAG-800-mAH | Table 6 | Lagrangian N = 800: m_AH/M | 0.177 | 0.177439 | 0.0005 | **PASS** |
| LAG-1600-t05 | Table 6 | Lagrangian N = 1600: t(lambda_h = 0.5) [s] | 1.7071 | 1.7071 | 5e-05 | **PASS** |
| LAG-1600-t015 | Table 6 | Lagrangian N = 1600: t(0.15) [s] | 1.7602 | 1.7602 | 5e-05 | **PASS** |
| LAG-1600-plunge | Table 6 | Lagrangian N = 1600: plunge [s] | 0.0531 | 0.0531229 | 0.0001 | **PASS** |
| LAG-1600-tAH | Table 6 | Lagrangian N = 1600: t_AH [s] | 1.7633 | 1.76334 | 5e-05 | **PASS** |
| LAG-1600-mAH | Table 6 | Lagrangian N = 1600: m_AH/M | 0.176 | 0.17607 | 0.0005 | **PASS** |
| LAG-ext-t05 | Table 6 | extrapolated t(0.5) [s] (difference ratios 3 and 4) | 1.702 | 1.70173 | 0.002 | **PASS** |
| LAG-ext-t015 | Table 6 | extrapolated t(0.15) [s] | 1.755 | 1.75482 | 0.003 | **PASS** |
| LAG-ratios | Table 6 caption | observed ratios of differences of t(0.5): 3.8, 3.3 | 3.8 | 3.76066 | 0.1 | **PASS** |
| LAG-ratio2 | Table 6 caption | ... second ratio | 3.3 | 3.27132 | 0.1 | **PASS** |
| LAG-order | Sec. 5.8 | convergence order of t(0.5): log2 of the ratio (1.7-1.9) | 1.8 | 1.70987 | 0.15 | **PASS** |
| LAG-diff1 | Sec. 5.8 | differences of t(0.5): N=200-400 [s] | 0.159 | 0.1587 | 0.0005 | **PASS** |
| LAG-diff2 | Sec. 5.8 | differences of t(0.5): N=400-800 [s] | 0.042 | 0.0422 | 0.0005 | **PASS** |
| LAG-diff3 | Sec. 5.8 | differences of t(0.5): N=800-1600 [s] | 0.013 | 0.0129 | 0.0005 | **PASS** |
| LAG-hom-resc-t05 | Table 6 | homologous, rescaled: t(0.5) [s] | 1.7083 | 1.7083 | 5e-05 | **PASS** |
| LAG-hom-resc-t015 | Table 6 | homologous, rescaled: t(0.15) [s] | 1.7601 | 1.7601 | 5e-05 | **PASS** |
| LAG-hom-resc-pl | Table 6 | homologous, rescaled: plunge [s] | 0.0518 | 0.0518 | 5e-05 | **PASS** |
| LAG-hom-raw-t05 | Table 6 | homologous, raw: t(0.5) [s] | 1.7595 | 1.7595 | 5e-05 | **PASS** |
| LAG-hom-raw-t015 | Table 6 | homologous, raw: t(0.15) [s] | 1.8133 | 1.8133 | 5e-05 | **PASS** |
| LAG-hom-raw-pl | Table 6 | homologous, raw: plunge [s] | 0.0538 | 0.0538 | 5e-05 | **PASS** |
| LAG-cmp1 | Sec. 5.8 | Lagrangian extrapolated t(0.5) vs rescaled homologous [%] | 0.4 | 0.386373 | 0.05 | **PASS** |
| LAG-cmp2 | Sec. 5.8 | ... vs raw [%] | 3.4 | 3.39508 | 0.05 | **PASS** |
| LAG-cmp-vel | Sec. 5.8 | rescaled homologous vs Lagrangian extrapolation, with the hand-off velocity read as dr/dt [%] | 0.5 | 0.498335 | 0.05 | **PASS** |
| LE-quick | Sec. 5.8 | regression: N = 100, 200 Jacobian eigenvalues recomputed now agree with the shipped lag_eigs.log (N = 100 ... 800) to 1e-5 relative |  | max rel diff 0.0e+00 |  | **PASS** |
| LAG-pl-resc | Sec. 5.8 | plunge (N = 1600) vs rescaled homologous [%] | -2.5 | -2.44821 | 0.15 | **PASS** |
| LAG-pl-raw | Sec. 5.8 | plunge (N = 1600) vs raw homologous [%] | 1.3 | 1.31827 | 0.15 | **PASS** |
| LAG-handoff-vel | Sec. 5.8 | hand-off velocity read as dr/dt shifts t(0.5) by [%] | -0.11 | -0.111407 | 0.01 | **PASS** |
| LAG-85 | Sec. 5.8 | fraction of t(lambda_h = 0.5) spent at lambda > 0.8 (homologous, Lagrangian N = 400, 800) | 0.85 | 0.85 | 0.03 | **PASS** |
| LAG-85-L | Sec. 5.8 | ... Lagrangian N = 800 | 0.85 | 0.848 | 0.03 | **PASS** |
| LAG-shell-min | Sec. 5.8 | plunge for shells enclosing up to 90 % of the mass: minimum [s] | 0.0522 | 0.0522 | 0.001 | **PASS** |
| LAG-shell-max | Sec. 5.8 | ... maximum [s] | 0.0541 | 0.0541 | 0.001 | **PASS** |
| LAG-resid | Sec. 5.8 | static residual at the outermost node as a fraction of local gravity | 0.5 | 0.52 | 0.1 | **PASS** |
| LAG-outer-mass | Sec. 5.8 | mass of the outermost zone / M (< 1e-4) | 0.0001 | 2.5e-05 | 5e-05 | **PASS** |
| TRAP-t05-AH | Sec. 5.8 | time from lambda_h = 0.5 to the first trapped shell, N = 1600 [s] | 0.0563 | 0.0562871 | 5e-05 | **PASS** |
| TRAP-dt200 | Sec. 5.8 | ... N = 200 [s] | 0.0568 | 0.0568257 | 5e-05 | **PASS** |
| TRAP-dt400 | Sec. 5.8 | ... N = 400 [s] | 0.0565 | 0.0565234 | 5e-05 | **PASS** |
| TRAP-dt800 | Sec. 5.8 | ... N = 800 [s] | 0.0563 | 0.0563401 | 5e-05 | **PASS** |
| TRAP-dt1600 | Sec. 5.8 | ... N = 1600 [s] | 0.0563 | 0.0562871 | 5e-05 | **PASS** |
| TRAP-after015 | Sec. 5.8 | t_AH - t(0.15), N = 1600 [s] | 0.003 | 0.00316412 | 0.0005 | **PASS** |
| TRAP-mAH | Sec. 5.8 | m_AH/M (N = 1600; paper 0.176 +- 0.003) | 0.176 | 0.17607 | 0.003 | **PASS** |
| TRAP-mAH-msun | Sec. 5.8 | m_AH [Msun] | 0.249 | 0.249422 | 0.004 | **PASS** |
| TRAP-remaining | Sec. 5.8 (iii) | mass not followed [Msun] | 1.17 | 1.16719 | 0.01 | **PASS** |
| TRAP-rhoc200 | Sec. 5.8 | central rest-mass density amplification, N = 200 [1e7] | 2.5 | 2.52949 | 0.06 | **PASS** |
| TRAP-rhoc400 | Sec. 5.8 | central rest-mass density amplification, N = 400 [1e7] | 2.9 | 2.85099 | 0.06 | **PASS** |
| TRAP-rhoc800 | Sec. 5.8 | central rest-mass density amplification, N = 800 [1e7] | 3.1 | 3.06749 | 0.06 | **PASS** |
| TRAP-rhoc1600 | Sec. 5.8 | central rest-mass density amplification, N = 1600 [1e7] | 3.2 | 3.20808 | 0.06 | **PASS** |
| TRAP-rhoc-g | Sec. 5.8 (i) | central density at the first trapped shell (N = 1600) [1e17 g/cm^3] | 7 | 7.41258 | 1.5 | **PASS** |
| TRAP-shell0.05 | Sec. 5.8 | radius of the shell enclosing 0.05 of the rest mass relative to its initial radius (N = 1600) | 0.0036 | 0.00357962 | 0.0006 | **PASS** |
| TRAP-shell0.2 | Sec. 5.8 | radius of the shell enclosing 0.2 of the rest mass relative to its initial radius (N = 1600) | 0.0054 | 0.00543904 | 0.0006 | **PASS** |
| TRAP-shell0.5 | Sec. 5.8 | radius of the shell enclosing 0.5 of the rest mass relative to its initial radius (N = 1600) | 0.0127 | 0.012749 | 0.0006 | **PASS** |
| TRAP-shell0.8 | Sec. 5.8 | radius of the shell enclosing 0.8 of the rest mass relative to its initial radius (N = 1600) | 0.03 | 0.0303667 | 0.0006 | **PASS** |
| TRAP-shell0.95 | Sec. 5.8 | radius of the shell enclosing 0.95 of the rest mass relative to its initial radius (N = 1600) | 0.06 | 0.0597006 | 0.0006 | **PASS** |
| TRAP-shell0.99 | Sec. 5.8 | radius of the shell enclosing 0.99 of the rest mass relative to its initial radius (N = 1600) | 0.091 | 0.0914813 | 0.0006 | **PASS** |
| TRAP-surf-lo | Sec. 5.8 | surface radius ratio, N = 200 | 0.167 | 0.167266 | 0.0005 | **PASS** |
| TRAP-surf-hi | Sec. 5.8 | surface radius ratio, N = 1600 | 0.171 | 0.171028 | 0.0035 | **PASS** |
| TRAP-viol200 | Sec. 5.8 | Hamiltonian-constraint violation (lambda_h > 0.2), N = 200 | 8e-08 | 8.0374e-08 | 4e-09 | **PASS** |
| TRAP-viol400 | Sec. 5.8 | Hamiltonian-constraint violation (lambda_h > 0.2), N = 400 | 2.7e-08 | 2.69653e-08 | 1.35e-09 | **PASS** |
| TRAP-viol800 | Sec. 5.8 | Hamiltonian-constraint violation (lambda_h > 0.2), N = 800 | 8.8e-09 | 8.81583e-09 | 4.4e-10 | **PASS** |
| TRAP-viol1600 | Sec. 5.8 | Hamiltonian-constraint violation (lambda_h > 0.2), N = 1600 | 3.9e-09 | 3.93842e-09 | 1.95e-10 | **PASS** |
| TRAP-violAH200 | Sec. 5.8 | ... at the first trapped shell, N = 200 | 2.6e-05 | 2.62228e-05 | 1.3e-06 | **PASS** |
| TRAP-violAH400 | Sec. 5.8 | ... at the first trapped shell, N = 400 | 5.7e-06 | 5.70914e-06 | 2.85e-07 | **PASS** |
| TRAP-violAH800 | Sec. 5.8 | ... at the first trapped shell, N = 800 | 1.2e-06 | 1.19255e-06 | 6e-08 | **PASS** |
| TRAP-violAH1600 | Sec. 5.8 | ... at the first trapped shell, N = 1600 | 2.4e-07 | 2.42707e-07 | 1.2e-08 | **PASS** |
| TRAP-gam2 | Sec. 5.8 | minimum Gamma^2 in every run (>= 0.41) | 0.41 | 0.415 | 0.005 | **PASS** |
| TRAP-cfl-t | Sec. 5.8 (scheme) | N = 200: t(lambda_h = 0.5), t(0.15), t_AH at Courant numbers 0.2, 0.4, 1.2 agree to (upper bound) [s] | 4e-08 | 3.64554e-08 | 5e-09 | **PASS** |
| TRAP-cfl-m | Sec. 5.8 (scheme) | N = 200: interpolated m_AH changes with the Courant number by [%] | 2 | 2.07426 | 0.5 | **PASS** |
| LAG-minDr | Sec. 5.8 (scheme) | minimum zone width stays positive in every run (N = 200 ... 1600), and Gamma^2 >= 0.41 |  | min Dr = 4.93 m, min Gamma^2 = 0.415 |  | **PASS** |
| LE-norm | Sec. 5.8 (scheme) | the linear tests fix the normalisation of sigma (extrapolated Lagrangian limit vs Rayleigh-Ritz) to (upper bound) [%] | 0.1 | 0.00103515 | 0.05 | **PASS** |
| LAG-fig-pdf | Figures | all six figures of the paper are regenerated by scripts (fig1_pekeris, figs [2], fig8_collapse_timescale, figs_fold, figs_friedman) |  | six PNG present |  | **PASS** |
| FIG6-last | Fig. 6 | Rayleigh-Ritz figure: five-function value [s^-2] | -28.55 | -28.5534 | 0.005 | **PASS** |
| FIG6-first | Fig. 6 | Rayleigh-Ritz figure: one-function value [s^-2] | -27.85 | -27.8524 | 0.005 | **PASS** |
| AP-eos | Sec. 3.7 (independent solver) | independent EOS (series for small x): first law, Gamma1, dP/dx checked to rounding (audit_pulsation/check_eos.py) |  | 6 points |  | **PASS** |
| TR-mch | tech. report Sec. 3.1 | Newtonian sequence, pure electron-gas EOS: maximum mass in the scan [Msun] | 1.42 | 1.4204 | 0.005 | **PASS** |
| TR-mch-an | tech. report Sec. 3.1 | analytic Chandrasekhar mass 5.836/mu_e^2 [Msun] | 1.459 | 1.459 | 0.0005 | **PASS** |
| TR-mch-3pc | tech. report Sec. 3.1 | scan maximum below the analytic value by (upper bound 3) [%] | 3 | 2.64565 | 0.5 | **PASS** |
| TR-Z0 | tech. report Sec. 3.2 | SEMF zero-density beta-stability charge for A = 12 | 5.77 | 5.773 | 0.005 | **PASS** |
| TR-Ye-lo | tech. report Sec. 3.2 | Y_e at low density | 0.488 | 0.4876 | 0.0005 | **PASS** |
| TR-Ye-hi | tech. report Sec. 3.2 | Y_e at 1.7e14 kg/m^3 | 0.35 | 0.3545 | 0.005 | **PASS** |
| TR-neu-max | tech. report Sec. 3.2 | neutronized Newtonian sequence: maximum mass [Msun] | 1.214 | 1.2136 | 0.0005 | **PASS** |
| TR-neu-rho | tech. report Sec. 3.2 | ... at rho_c [1e12 kg/m^3] | 2.3 | 2.305 | 0.05 | **PASS** |
| TR-neu-branch | tech. report Sec. 3.2 | rho_c = 5e12 lies on the unstable branch dM/drho_c < 0 of the neutronized sequence |  | M(5e12) = 1.1977, max = 1.2136 |  | **PASS** |
| TR-lin | tech. report Sec. 6 | fundamental Galerkin eigenvalue of the rho_c = 5e12 star [s^-2] | -18.68 | -18.6773 | 0.005 | **PASS** |
| TR-lin-ef | tech. report Sec. 6 | e-folding time [s] | 0.23 | 0.23139 | 0.005 | **PASS** |
| TR-lin-sig | tech. report Sec. 6 | growth rate sqrt(18.68) [1/s] | 4.32 | 4.32173 | 0.005 | **PASS** |
| TR-lin-sig4 | tech. report Sec. 7.4 | sigma used for the pure-mode seed [1/s] | 4.3217 | 4.32173 | 5e-05 | **PASS** |
| TR-lin-over | tech. report Sec. 6 | every overtone through node number 14 is stable (positive, node counts 1 ... 14) |  | 15 modes |  | **PASS** |
| TR-RQ | tech. report Sec. 6 | independent Rayleigh quotient with xi = alpha r [s^-2] | -16.65 | -16.654 | 0.0051 | **PASS** |
| TR-RQ-bound | tech. report Sec. 6 | the trial-function value is less negative than the Galerkin eigenvalue (Rayleigh-Ritz upper bound) |  | -16.654 > -18.67734 |  | **PASS** |
| TR-stable | tech. report Sec. 6 | stable control rho_c = 1e11: fundamental omega^2 [s^-2] | 4.65 | 4.65322 | 0.005 | **PASS** |
| TR-stable-T | tech. report Sec. 6 | ... oscillation period [s] | 2.91 | 2.9127 | 0.005 | **PASS** |
| TR-stable-over | tech. report Sec. 6 | all overtones of the stable control are positive |  | 8 modes |  | **PASS** |
| TR-evolve2-rhs | tech. report Sec. 2 | right-hand side of the pure-EOS Lagrangian code is finite (no NaN), dR/dt = 0 and dU/dt of Newtonian size on the initial data |  | see tr_evolve2.log |  | **PASS** |
| TR-eul-400 | tech. report Sec. 5.3 (table) | Eulerian resolution scan, N = 400: rho_c(t)/rho_c(0) at t = 400 ms | 0.289 | 0.28898 | 0.0005 | **PASS** |
| TR-eul-800 | tech. report Sec. 5.3 (table) | ... N = 800 (euler_convergence.py) | 0.581 | 0.58063 | 0.0005 | **PASS** |
| TR-eul-800b | tech. report Sec. 5.3 (table) | ... N = 800 (refactored euler_module.py, study1.py) | 0.581 | 0.58111 | 0.0005 | **PASS** |
| TR-eul-1600 | tech. report Sec. 5.3 (table) | ... N = 1600 (refactored euler_module.py, study1.py) | 1.079 | 1.07939 | 0.0006 | **PASS** |
| TR-eul-mass | tech. report Sec. 5.3 (table) | mass conservation exact to 6 digits (study1.py, N = 800, 1600) |  | [1.0, 1.0] |  | **PASS** |
| TR-eul-3200 | tech. report Sec. 5.3 (table) | ... N = 3200 (study2.py) | 1.781 | 1.78056 | 0.0005 | **PASS** |
| TR-eul-6400 | tech. report Sec. 5.3 (table) | ... N = 6400 (study2.py) | 2.565 | 2.5646 | 0.0005 | **PASS** |
| TR-eul-mass2 | tech. report Sec. 5.3 (table) | mass conservation exact to 6 digits (study2.py, N = 3200, 6400) |  | [1.0, 1.0] |  | **PASS** |
| TR-ctl-800 | tech. report Sec. 5.3 (control) | stable control rho_c = 1e11, N = 800: rho_c(t)/rho_c(0) at t = 400 ms | 0.977 | 0.977353 | 0.0005 | **PASS** |
| TR-ctl-1600 | tech. report Sec. 5.3 (control) | ... N = 1600 | 0.99 | 0.989743 | 0.0005 | **PASS** |
| TR-ctl-3200 | tech. report Sec. 5.3 (control) | ... N = 3200 | 0.996 | 0.996136 | 0.0005 | **PASS** |
| TR-ctl-conv | tech. report Sec. 5.3 (control) | the control converges towards 1 with the deviation shrinking by a consistent factor per doubling |  | deviations 0.0226 0.0103 0.0039, factors 2.21 2.65 |  | **PASS** |
| TR-eul-dir | tech. report Sec. 5.3 | the marginal star changes direction with resolution: net rarefaction at N = 400, 800, net compression at N = 1600, 3200, 6400 (monotonic in N, 16-fold range) |  | 0.289 0.581 1.079 1.781 2.565 |  | **PASS** |
| TR-s2-800 | tech. report Sec. 7.2 | seeded (eps = 1e-4, v = 0), N = 800: runs to t = 3 s; central density falls to ~1e-4 of its initial value by t ~ 2.6 s, then shows renewed structure |  | DONE, rho_c/rho_c0 = 1.27e-04 at 2.6 s, 2.29e-04 at 3.00 s |  | **PASS** |
| TR-s2-1600-t | tech. report Sec. 7.2 | seeded, N = 1600: the run ends in collapse at t [s] | 0.618 | 0.618 | 0.0006 | **PASS** |
| TR-s2-1600 | tech. report Sec. 7.2 | seeded, N = 1600: reverses within ~0.4 s and runs away; central density rises by a factor of about 4e6 (>= 1e6) up to the stop |  | BLOWUP, rho_c/rho_c0 = 4.06e+06 at 0.618 s, 1.071 at 0.3 s |  | **PASS** |
| TR-s2-split | tech. report Sec. 7.2 | identical seed, N = 800 and N = 1600 end in qualitatively different fates (rarefaction vs collapse) |  | N = 800: 2.29e-04, N = 1600: 4.06e+06 |  | **PASS** |
| TR-s2-strong | tech. report Sec. 7.2 | 100x stronger seed (eps = 1e-2): both resolutions are driven into an unphysical regime (|v| > 1e12 m/s) within ~0.2 s |  | t = 0.232 s (N = 800), 0.221 s (N = 1600) |  | **PASS** |
| TR-pm-800 | tech. report Sec. 7.4 | pure growing-mode seed, N = 800: same rarefaction as the v = 0 seed (rho_c/rho_c0 ~ 1e-4 at 2.6 s, within 10 % of the v = 0 value) |  | 1.270e-04 vs 1.270e-04 |  | **PASS** |
| TR-pm-1600 | tech. report Sec. 7.4 | pure growing-mode seed, N = 1600: runs away into collapse (t ~ 0.62 s) and passes a central compression of 2.45 |  | BLOWUP at 0.616 s, max sampled compression 3.75e+06 |  | **PASS** |
| TR-pm-sigma | tech. report Sec. 7.4 | sigma of the pure-mode seed [1/s] | 4.3217 | 4.32173 | 5e-05 | **PASS** |
| FT-json-regen | Table 3 | the 14 rows of the shipped fold_dyn2.json (source of Table 3) are reproduced by the shipped fold_dyn2.py (relative differences < 1e-9) |  | 14 rows, largest relative difference 0.0e+00 |  | **PASS** |
| PI-old-estimator | Sec. 5.6 | the early pole estimator (initial slope with the wrong sign), which is still the PI block of fold_dyn2.py, gives pole positions between 3.3 and 3.5 that depend on the start point (shipped PI field and regenerated copy), not 3.41167 |  | shipped: [3.3157, 3.3454, 3.4659, 3.4534]; regenerated: [3.3206, 3.3503, 3.4708, 3.4583] |  | **PASS** |
| TR-guard-equiv-800 | tech. report Sec. 7.2 (README of seeded_evolution) | guarded and original seeded module: identical logged states up to the stop (eps = 1e-2, N = 800) |  | see log |  | **PASS** |
| TR-guard-equiv-1600 | tech. report Sec. 7.2 (README of seeded_evolution) | guarded and original seeded module: identical logged states up to the stop (eps = 1e-4, N = 1600, the run of Sec. 7.2) |  | see log |  | **PASS** |
| LG-L-pure | tech. report Sec. 4 | L = Newtonian radius of the report profile, pure EOS [km] | 1601.9 | 1601.9 | 0.05 | **PASS** |
| LG-L-neu | tech. report Sec. 4 | L, neutronised EOS [km] | 1615.8 | 1615.8 | 0.05 | **PASS** |
| LG-Lc-pure | tech. report Sec. 4 | L/c, pure EOS [ms] | 5.343 | 5.34336 | 0.0005 | **PASS** |
| LG-Lc-neu | tech. report Sec. 4 | L/c, neutronised EOS [ms] | 5.39 | 5.38973 | 0.0005 | **PASS** |
| LG-wall | tech. report Sec. 4 | wall time of the shipped run_all.sh run on two cores [min] (machine dependent: always the SHIPPED log) | 20 | 19.6167 | 0.5 | **PASS** |
| LG-ms-resid | tech. report Sec. 4.1 | residuals of the identities I1-I3 (Einstein tensor, exact 2-jets) [round-2 note: the log prints 1.51e-13 for I2] |  | 1.51e-13 <= 1.5e-13 (+ half unit 5e-15) |  | **PASS** |
| LG-ms-neg | tech. report Sec. 4.1 | a deliberately wrong coefficient is detected at (smallest residual over the points) | 0.0004 | 0.000369 | 5e-05 | **PASS** |
| LG-ms-I56 | tech. report Sec. 4.1 | residuals of I5 (lapse equation) and I6 (continuity equation) from the divergence of the perfect-fluid stress tensor |  | 3.55e-15 <= 4e-15 (+ half unit 5e-16) |  | **PASS** |
| LG-ms-neg5 | tech. report Sec. 4.1 | a coefficient 1 % wrong in the lapse equation is detected at (median residual) | 0.003 | 0.0033 | 0.0005 | **PASS** |
| LG-ms-I4 | tech. report Sec. 4.1 | identity I4 (momentum equation follows from the others) residual | 6e-17 | 5.55e-17 | 5e-18 | **PASS** |
| LG-ms-40 | tech. report Sec. 4.1 | identities checked at 40 random points, log ends with PASS |  | see logs/verify_ms_equations.log |  | **PASS** |
| LG-ms-gdot | tech. report Sec. 4.1 | dGamma/dt / Gamma at generic points | 0.14 | 0.14 | 0.005 | **PASS** |
| LG-2mR-pure | tech. report Sec. 4.1 | 2m/R of the pure star (quoted 2.6e-3 as the size of the whole relativistic correction) | 0.0026 | 0.002604 | 5e-05 | **PASS** |
| LG-gamma-iter | tech. report Sec. 4.2 | Gamma^2 = 1 + U^2 - 2m/R by three fixed-point iterations, error 2e-9 (change at the third iteration, both EOS: 1.5e-9 ... 2.6e-9) |  | [2.5023677650344357e-09, 1.6568202543721537e-09] |  | **PASS** |
| LG-cfl | tech. report Sec. 4.2 | CFL number 0.4 and viscosity q = 2 rho_N (dU)^2 are the defaults / the labels of the runs |  | src/harness.py, RESULTS.md |  | **PASS** |
| LG-tov-M-pure | tech. report Sec. 4.3 (a) | TOV mass, rho_c = 5e12 [Msun], pure | 1.4101 | 1.4101 | 5e-05 | **PASS** |
| LG-tov-R-pure | tech. report Sec. 4.3 (a) | TOV radius [km], pure | 1599.3 | 1599.3 | 0.05 | **PASS** |
| LG-tov-c-pure | tech. report Sec. 4.3 (a) | 2M/R, pure | 0.0026 | 0.002604 | 5e-06 | **PASS** |
| LG-tov-M-neu | tech. report Sec. 4.3 (a) | TOV mass, rho_c = 5e12 [Msun], neu | 1.1901 | 1.1901 | 5e-05 | **PASS** |
| LG-tov-R-neu | tech. report Sec. 4.3 (a) | TOV radius [km], neu | 1613.7 | 1613.7 | 0.05 | **PASS** |
| LG-tov-c-neu | tech. report Sec. 4.3 (a) | 2M/R, neu | 0.00218 | 0.002179 | 5e-06 | **PASS** |
| LG-static-50 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, pure EOS: max |R_j(10) - R_j(0)|/L, N = 50 | 0.00021 | 0.000214748 | 5e-06 | **PASS** |
| LG-static-100 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, pure EOS: max |R_j(10) - R_j(0)|/L, N = 100 | 5.4e-05 | 5.38192e-05 | 5e-07 | **PASS** |
| LG-static-200 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, pure EOS: max |R_j(10) - R_j(0)|/L, N = 200 | 1.3e-05 | 1.34378e-05 | 5e-07 | **PASS** |
| LG-static-400 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, pure EOS: max |R_j(10) - R_j(0)|/L, N = 400 | 3.3e-06 | 3.33358e-06 | 5e-08 | **PASS** |
| LG-static-800 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, pure EOS: max |R_j(10) - R_j(0)|/L, N = 800 | 8.1e-07 | 8.06324e-07 | 5e-09 | **PASS** |
| LG-static-neu-50 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, neutronised EOS: max |R_j(10) - R_j(0)|/L, N = 50 | 0.00022 | 0.000219767 | 5e-06 | **PASS** |
| LG-static-neu-100 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, neutronised EOS: max |R_j(10) - R_j(0)|/L, N = 100 | 5.5e-05 | 5.51512e-05 | 5e-07 | **PASS** |
| LG-static-neu-200 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, neutronised EOS: max |R_j(10) - R_j(0)|/L, N = 200 | 1.4e-05 | 1.37863e-05 | 5e-07 | **PASS** |
| LG-static-neu-400 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, neutronised EOS: max |R_j(10) - R_j(0)|/L, N = 400 | 3.4e-06 | 3.42658e-06 | 5e-08 | **PASS** |
| LG-static-neu-800 | tech. report Sec. 4.3 (a) [after edit R2-5] | pointwise TOV data, neutronised EOS: max |R_j(10) - R_j(0)|/L, N = 800 | 8.4e-07 | 8.36036e-07 | 5e-09 | **PASS** |
| LG-static-ord-disp | tech. report Sec. 4.3 (a) | observed order of the displacement (finest pair; pure/neutronised, T = 10, 30) |  | computed 2.0351 ... 2.1032, quoted 2.0 ... 2.1 |  | **PASS** |
| LG-static-ord-rho | tech. report Sec. 4.3 (a) | observed order of the central density (finest pair) |  | computed 1.7126 ... 1.8027, quoted 1.7 ... 1.8 |  | **PASS** |
| LG-static-db-R | tech. report Sec. 4.3 (a) | discretely balanced equilibrium, radius deviation in units of L (both EOS, all N, t = 10 and 30) |  | 9.691e-15 <= 1e-14 (+ half unit 5e-15) |  | **PASS** |
| LG-static-db-rho30 | tech. report Sec. 4.3 (a) | discretely balanced equilibrium, central-density deviation at t <= 30 (both EOS, all N) |  | 1.1656e-11 <= 1.2e-11 (+ half unit 5e-13) |  | **PASS** |
| LG-static-db-rho150 | tech. report Sec. 4.3 (a) | discretely balanced equilibrium, central-density deviation at t <= 150 (test (g), both EOS, all N) |  | 5.6861e-11 <= 5.7e-11 (+ half unit 5e-13) |  | **PASS** |
| LG-lin-rel-5e12 | tech. report Sec. 4.3 (b) | scheme vs independent Sturm-Liouville solver, fundamental mode, rho_c = 5e12 (both EOS): relative difference |  | computed 2.1677e-05 ... 3.4183e-05, quoted 2e-5 ... 3.4e-5 |  | **PASS** |
| LG-lin-rel-1e11 | tech. report Sec. 4.3 (b) | same, rho_c = 1e11 (both EOS): 6e-6 |  | [5.715650177353814e-06, 5.756442193928895e-06] |  | **PASS** |
| LG-lin-order | tech. report Sec. 4.3 (b) | observed order of the fundamental eigenvalue after Richardson extrapolation: 1.99 (all four cases between 1.98 and 2.01; the neutronised 5e12 case is 2.00) |  | 2.0018 1.9987 1.9938 1.9798 1.9946 1.9889 1.9944 1.9899 |  | **PASS** |
| LG-lin-newt | tech. report Sec. 4.3 (b) | Newtonian limit vs the report Galerkin solver: relative difference of the fundamental mode (4 cases) |  | computed 2.3363e-05 ... 7.6406e-05, quoted 2e-5 ... 8e-5 |  | **PASS** |
| LG-lin-time | tech. report Sec. 4.3 (b) | time-domain frequency of the stable 1e11 star vs eigenvalue, relative difference (both EOS, all N, dt) |  | 1.225e-06 <= 1.3e-6 (+ half unit 5e-08) |  | **PASS** |
| LG-sigma-200 | tech. report Sec. 4.3 (b) | fitted growth rate of the seeded unstable mode, N = 200, eps = 1e-4 [c/L] | 0.02419 | 0.0241894 | 5e-06 | **PASS** |
| LG-sigma-400 | tech. report Sec. 4.3 (b) | ... N = 400 | 0.02417 | 0.0241712 | 5e-06 | **PASS** |
| LG-sigma-ref | tech. report Sec. 4.3 (b) | growth rate of the independent solver, sigma = sqrt(-omega^2) [c/L] | 0.024189 | 0.0241894 | 5e-07 | **PASS** |
| LG-sigma-dev | tech. report Sec. 4.3 (b) | relative deviation of the fitted growth rate from the independent solver (N = 200, 400; eps = 1e-4) |  | 0.00075018 <= 1e-3 (+ half unit 0.0005) |  | **PASS** |
| LG-dt-half-rho | tech. report Sec. 4.3 (c) | halving dt changes rho_c (relative) at t = 20, 40, 50, N = 200, 400, both EOS |  | 2.8667e-05 <= 3e-5 (+ half unit 5e-06) |  | **PASS** |
| LG-dt-half-cross | tech. report Sec. 4.3 (c) | halving dt changes the edge-crossing time [%] (pure, N = 200, 400) |  | 0.259 <= 0.26 (+ half unit 0.005) |  | **PASS** |
| LG-ord-times | tech. report Sec. 4.3 (c) | observed orders of the times to reach 2, 5, 10, 100, 300 x rho_c0 (RESULTS.md (h)): 2.0 ... 2.2 |  | orders 2.14 2.03 2.02 2.04 2.00 2.14 2.19 |  | **PASS** |
| LG-ord-rho40-pure | tech. report Sec. 4.3 (c) [after edit R2-3] | observed order of rho_c(t = 40), pure, N = 400, 800, 1600 (RESULTS.md (c); Richardson limits assume p = 2) | 1.65 | 1.65 | 0.005 | **PASS** |
| LG-ord-rho40-neu | tech. report Sec. 4.3 (c) [after edit R2-3] | ... neutronised | 1.62 | 1.62 | 0.005 | **PASS** |
| LG-mass-drift | tech. report Sec. 4.3 (d) | gravitational mass drift over 50 L/c (report data, N = 100 ... 800, both EOS) [max 2.10e-11, printed with 3 digits] |  | 2.0995e-11 <= 2e-11 (+ half unit 5e-12) |  | **PASS** |
| LG-energy-drift | tech. report Sec. 4.3 (d) | Newtonian-limit energy drift / |W| |  | 4.916e-08 <= 5e-8 (+ half unit 5e-09) |  | **PASS** |
| LG-energy-order | tech. report Sec. 4.3 (d) | observed orders of the discrete energy-equation defect |  | computed 1.9924 ... 2.0928, quoted 1.99 ... 2.09 |  | **PASS** |
| LG-gamma-res | tech. report Sec. 4.3 (d) | Gamma^2 = 1 + U^2 - 2m/R holds to (residual; quoted as 1e-16, printed 1.11e-16) |  | 1.1102e-16 <= 1.5e-16 (+ half unit 5e-18) |  | **PASS** |
| LG-ctl-rmax | tech. report Sec. 4.3 (control) | the report's stable_crosscheck.py (and control_test.py, euler_module.py) build the star with rmax = 3e6 m (source of ref/ and of the repository) |  | ref/stable_crosscheck.py, linear_perturbation/stable_crosscheck.py, control_test.py |  | **PASS** |
| LG-ctl-dens-neu | tech. report Sec. 4.3 (control) | density at the cut as % of rho_c (neutronised) | 3.6 | 3.56 | 0.05 | **PASS** |
| LG-ctl-dens-pure | tech. report Sec. 4.3 (control) | ... (pure) | 4.1 | 4.13 | 0.05 | **PASS** |
| LG-ctl-R-neu | tech. report Sec. 4.3 (control) | radius of the full star [km] (neutronised), quoted to 10 km | 4250 | 4246.3 | 5 | **PASS** |
| LG-ctl-R-pure | tech. report Sec. 4.3 (control) | ... (pure) | 4330 | 4327.9 | 5 | **PASS** |
| LG-ctl-R-trunc | tech. report Sec. 4.3 (control) [after edit R2-1] | radius at which the report control profile is cut (rmax = 3e6 m) [km] | 3000 | 3000 | 0.5 | **PASS** |
| LG-ctl-M | tech. report Sec. 4.3 (control) | mass of the truncated (neutronised) object [Msun] | 0.99 | 0.989 | 0.005 | **PASS** |
| LG-ctl-w-trunc | tech. report Sec. 4.3 (control) | report's Galerkin solver on the truncated (neutronised) profile: omega^2 [s^-2] | 4.638 | 4.638 | 0.0005 | **PASS** |
| LG-ctl-w-full-neu | tech. report Sec. 4.3 (control) | same solver, full star, neutronised [s^-2] | 1.756 | 1.7557 | 0.0005 | **PASS** |
| LG-ctl-w-full-pure | tech. report Sec. 4.3 (control) | same solver, full star, pure [s^-2] | 1.953 | 1.9533 | 0.0005 | **PASS** |
| LG-ctl-w-gr-neu | tech. report Sec. 4.3 (control) | scheme (GR), full star, neutronised [s^-2] | 1.743 | 1.74305 | 0.0005 | **PASS** |
| LG-ctl-w-gr-pure | tech. report Sec. 4.3 (control) | scheme (GR), full star, pure [s^-2] | 1.94 | 1.9397 | 0.0005 | **PASS** |
| LG-ctl-period-full | tech. report Sec. 4.3 (control) [after edit R2-1] | oscillation period of the full (neutronised) 1e11 star [s] | 4.74 | 4.742 | 0.005 | **PASS** |
| LG-eq-nocross | tech. report Sec. 4.4 | exact (discretely balanced) TOV configurations keep R' = 1 to rounding: t <= 150 for N = 200, 400 (test g), t <= 60 for N = 100 ... 800 (test c); no crossing |  | 12 runs, largest |min R' - 1| = 4.1e-10 |  | **PASS** |
| LG-w2-pure | tech. report Sec. 4.4 | omega^2 of the pure 5e12 equilibrium = +9.273 s^-2 (scheme 9.2726, independent solver 9.2723: the quoted value agrees with one of the two to half a unit of the last digit and the two agree to 5e-4) |  | scheme 9.27264, Sturm-Liouville 9.27232 |  | **PASS** |
| LG-w2-neu | tech. report Sec. 4.4 | omega^2 of the neutronised 5e12 equilibrium = -20.144 s^-2 (independent solver -20.1438, scheme -20.1434) |  | scheme -20.14337, Sturm-Liouville -20.14380 |  | **PASS** |
| LG-efold-s | tech. report Sec. 4.4 / abstract | e-folding time 1/sqrt(-omega^2) of the neutronised equilibrium [s] | 0.223 | 0.222807 | 0.0005 | **PASS** |
| LG-efold-Lc | tech. report Sec. 4.4 | ... in units of L/c (neutronised) | 41.3 | 41.3392 | 0.05 | **PASS** |
| LG-w2-newt-scheme | tech. report Sec. 4.4 | Newtonian-limit omega^2, neutronised: scheme [s^-2] | -18.721 | -18.7214 | 0.0005 | **PASS** |
| LG-w2-newt-galerkin | tech. report Sec. 4.4 | ... report Galerkin solver (Richardson) [s^-2] | -18.721 | -18.7209 | 0.0005 | **PASS** |
| LG-w2-newt-diff | tech. report Sec. 4.4 | difference between -18.721 and the report's -18.68 [%] | 0.2 | 0.219486 | 0.05 | **PASS** |
| LG-gr-corr | tech. report Sec. 4.4 | general-relativistic correction of omega^2 (neutronised) [%] | 7.6 | 7.60048 | 0.05 | **PASS** |
| LG-2mR-neu | tech. report Sec. 4.4 | 2M/R of the neutronised star | 0.0022 | 0.002179 | 5e-05 | **PASS** |
| LG-noise-38 | tech. report Sec. 4.4 | amplification of rounding noise exp(sigma * 150) over t = 150, sigma = 0.0242 c/L | 38 | 37.7128 | 0.5 | **PASS** |
| LG-drift-p10 | tech. report Sec. 4.4 (initial data) | drift of the report data, pure, t = 10 [L] | 0.0014 | 0.00141388 | 5e-05 | **PASS** |
| LG-drift-n10 | tech. report Sec. 4.4 (initial data) | ... neutronised, t = 10 | 0.0011 | 0.00114295 | 5e-05 | **PASS** |
| LG-drift-p30 | tech. report Sec. 4.4 (initial data) | ... pure, t = 30 | 0.013 | 0.0129438 | 0.0005 | **PASS** |
| LG-drift-n30 | tech. report Sec. 4.4 (initial data) [after edit R2-5] | ... neutronised, t = 30 | 0.011 | 0.010507 | 0.0005 | **PASS** |
| LG-drift-N | tech. report Sec. 4.4 (initial data) | the drift is independent of N (observed order 0): |order| < 0.05 for the four rows (N = 400 -> 800) |  | -0.001 -0.001 -0.001 -0.002 |  | **PASS** |
| LG-rho-p10 | tech. report Sec. 4.4 (initial data) | central density rise at t = 10, pure [%] | 1.76 | 1.75566 | 0.005 | **PASS** |
| LG-rho-n10 | tech. report Sec. 4.4 (initial data) | ... neutronised | 1.59 | 1.58755 | 0.005 | **PASS** |
| LG-rho-p30 | tech. report Sec. 4.4 (initial data) | at t = 30, pure [%] | 8.8 | 8.79081 | 0.05 | **PASS** |
| LG-rho-n30 | tech. report Sec. 4.4 (initial data) | ... neutronised | 8.7 | 8.71575 | 0.05 | **PASS** |
| LG-acc-pure | tech. report Sec. 4.4 (initial data) | median initial acceleration / g of the report data, pure, N = 800 | -0.0046 | -0.00462046 | 5e-05 | **PASS** |
| LG-acc-neu | tech. report Sec. 4.4 (initial data) | ... neutronised | -0.004 | -0.00396686 | 5e-05 | **PASS** |
| LG-acc-pct-lo | tech. report Sec. 4.4 (initial data) | smaller of the two accelerations as % of g ("0.4-0.5%") | 0.4 | 0.396686 | 0.05 | **PASS** |
| LG-acc-pct-hi | tech. report Sec. 4.4 (initial data) | larger ...  | 0.5 | 0.462046 | 0.05 | **PASS** |
| LG-acc-N | tech. report Sec. 4.4 (initial data) | the acceleration of the report data is independent of N (N = 400 and 800 agree to 1 %) |  | [(-0.004628354321393151, -0.004620458257380758), (-0.003973987561675287, -0.003966860780831498)] |  | **PASS** |
| LG-acc-tov | tech. report Sec. 4.4 (initial data) | pointwise TOV data: median acceleration < 1e-4 g and maximum falling as h^2 (ratio 3.5 ... 4.5 per doubling) |  | max median 4.2e-05, ratios ['3.70', '3.84', '3.68', '3.83'] |  | **PASS** |
| LG-eq-dist-pure | tech. report Sec. 4.4 (initial data) [after edit R2-4] | nearest GR equilibrium with the same zone masses (pure): the report central density is 55 % below it; that equilibrium is 122-123 % above the TOV star (N = 200: 122.3, N = 400: 122.7) |  | report/equilibrium - 1 = [-0.5497, -0.5509]; equilibrium/TOV - 1 = [1.2233, 1.2274] |  | **PASS** |
| LG-eq-dist-neu | tech. report Sec. 4.4 (initial data) | neutronised: the report central density is 29 % above the nearest equilibrium, which is 23 % below the TOV star (N = 200, 400) |  | [(0.2911, -0.2274), (0.2938, -0.2276)] |  | **PASS** |
| LG-tx-pure-200 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, pure, N = 200 [L/c], text | 59.5 | 59.5462 | 0.05 | **PASS** |
| LG-tx-pure-400 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, pure, N = 400 [L/c], text | 59.5 | 59.4807 | 0.05 | **PASS** |
| LG-tx-pure-800 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, pure, N = 800 [L/c], text | 59.3 | 59.2802 | 0.05 | **PASS** |
| LG-tx-pure-1600 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, pure, N = 1600 [L/c], text | 59 | 58.9655 | 0.05 | **PASS** |
| LG-tx-neu-200 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, neu, N = 200 [L/c], text | 66.9 | 66.9323 | 0.05 | **PASS** |
| LG-tx-neu-400 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, neu, N = 400 [L/c], text | 66.4 | 66.3932 | 0.05 | **PASS** |
| LG-tx-neu-800 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, neu, N = 800 [L/c], text | 66.4 | 66.3709 | 0.05 | **PASS** |
| LG-tx-neu-1600 | tech. report Sec. 4.4 (first R' <= 0) | time of the first loss of ordering, neu, N = 1600 [L/c], text | 66.1 | 66.139 | 0.05 | **PASS** |
| LG-T2-t-pure-200 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, pure, N = 200 [L/c] | 59.546 | 59.5462 | 0.0005 | **PASS** |
| LG-T2-t-pure-400 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, pure, N = 400 [L/c] | 59.481 | 59.4807 | 0.0005 | **PASS** |
| LG-T2-t-pure-800 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, pure, N = 800 [L/c] | 59.28 | 59.2802 | 0.0005 | **PASS** |
| LG-T2-t-pure-1600 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, pure, N = 1600 [L/c] | 58.965 | 58.9655 | 0.0005 | **PASS** |
| LG-T2-t-neu-200 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, neu, N = 200 [L/c] | 66.932 | 66.9323 | 0.0005 | **PASS** |
| LG-T2-t-neu-400 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, neu, N = 400 [L/c] | 66.393 | 66.3932 | 0.0005 | **PASS** |
| LG-T2-t-neu-800 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, neu, N = 800 [L/c] | 66.371 | 66.3709 | 0.0005 | **PASS** |
| LG-T2-t-neu-1600 | tech. report Table 5 (tab:lag-edge) | inviscid first loss of ordering, neu, N = 1600 [L/c] | 66.139 | 66.139 | 0.0005 | **PASS** |
| LG-T2-zone-pure-200 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, pure, N = 200 | 198 | 198 | 0 | **PASS** |
| LG-T2-zone-pure-400 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, pure, N = 400 | 398 | 398 | 0 | **PASS** |
| LG-T2-zone-pure-800 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, pure, N = 800 | 794 | 794 | 0 | **PASS** |
| LG-T2-zone-pure-1600 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, pure, N = 1600 | 1585 | 1585 | 0 | **PASS** |
| LG-T2-zone-neu-200 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, neu, N = 200 | 199 | 199 | 0 | **PASS** |
| LG-T2-zone-neu-400 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, neu, N = 400 | 397 | 397 | 0 | **PASS** |
| LG-T2-zone-neu-800 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, neu, N = 800 | 792 | 792 | 0 | **PASS** |
| LG-T2-zone-neu-1600 | tech. report Table 5 (tab:lag-edge) | zone of the crossing, neu, N = 1600 | 1586 | 1586 | 0 | **PASS** |
| LG-T2-frac-pure-200 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, pure, N = 200 | 8e-07 | 7.52121e-07 | 5e-08 | **PASS** |
| LG-T2-frac-pure-400 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, pure, N = 400 | 1e-07 | 1.3224e-07 | 5e-08 | **PASS** |
| LG-T2-frac-pure-800 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, pure, N = 800 | 1e-06 | 1.31749e-06 | 5e-07 | **PASS** |
| LG-T2-frac-pure-1600 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, pure, N = 1600 | 3e-06 | 3.07217e-06 | 5e-07 | **PASS** |
| LG-T2-frac-neu-200 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, neu, N = 200 | 0 | 0 | 5e-08 | **PASS** |
| LG-T2-frac-neu-400 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, neu, N = 400 | 7e-07 | 7.43478e-07 | 5e-08 | **PASS** |
| LG-T2-frac-neu-800 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, neu, N = 800 | 3e-06 | 3.03333e-06 | 5e-07 | **PASS** |
| LG-T2-frac-neu-1600 | tech. report Table 5 (tab:lag-edge) | mass fraction outside the crossing zone, neu, N = 1600 | 3e-06 | 2.51738e-06 | 5e-07 | **PASS** |
| LG-T2-ratio-200 | tech. report Table 5 (tab:lag-edge) | ratio t_neu/t_pure, N = 200 | 1.124 | 1.12404 | 0.0005 | **PASS** |
| LG-T2-ratio-400 | tech. report Table 5 (tab:lag-edge) | ratio t_neu/t_pure, N = 400 | 1.116 | 1.11621 | 0.0005 | **PASS** |
| LG-T2-ratio-800 | tech. report Table 5 (tab:lag-edge) | ratio t_neu/t_pure, N = 800 | 1.12 | 1.11961 | 0.0005 | **PASS** |
| LG-T2-ratio-1600 | tech. report Table 5 (tab:lag-edge) | ratio t_neu/t_pure, N = 1600 | 1.122 | 1.12166 | 0.0005 | **PASS** |
| LG-edge-frac | tech. report Sec. 4.4 (first R' <= 0) | mass outside the crossing zone, all runs [less than 3.1e-6] |  | 3.0722e-06 <= 3.1e-6 (+ half unit 5e-08) |  | **PASS** |
| LG-rho40-pure-400 | tech. report Sec. 4.4 (first R' <= 0) | pure rho_c/rho_c0 - 1 at t = 40, N = 400 | 0.1439 | 0.143902 | 5e-06 | **PASS** |
| LG-rho40-pure-800 | tech. report Sec. 4.4 (first R' <= 0) | pure rho_c/rho_c0 - 1 at t = 40, N = 800 | 0.14376 | 0.143759 | 5e-06 | **PASS** |
| LG-rho40-pure-1600 | tech. report Sec. 4.4 (first R' <= 0) | pure rho_c/rho_c0 - 1 at t = 40, N = 1600 | 0.14371 | 0.143713 | 5e-06 | **PASS** |
| LG-rho40-pure-R | tech. report Sec. 4.4 (first R' <= 0) | ... Richardson limit (p = 2) | 0.1437 | 0.143698 | 5e-06 | **PASS** |
| LG-cross-dmin | tech. report Sec. 4.4 (first R' <= 0) | smallest change of the crossing time per doubling of N [L/c] | 0.02 | 0.0222494 | 0.005 | **PASS** |
| LG-cross-dmax | tech. report Sec. 4.4 (first R' <= 0) | largest change per doubling [L/c] | 0.54 | 0.539146 | 0.005 | **PASS** |
| LG-mid-pure | tech. report Sec. 4.4 / Established | t_cross pure: mid-point of the range over N [L/c] | 59.3 | 59.2559 | 0.05 | **PASS** |
| LG-hw-pure | tech. report Sec. 4.4 / Established | ... half width | 0.3 | 0.290372 | 0.05 | **PASS** |
| LG-mid-neu | tech. report Sec. 4.4 / Established | t_cross neutronised: mid-point [L/c] | 66.5 | 66.5357 | 0.05 | **PASS** |
| LG-hw-neu | tech. report Sec. 4.4 / Established | ... half width | 0.4 | 0.396644 | 0.05 | **PASS** |
| LG-mid-RESULTS | tech. report Sec. 4.4 | RESULTS.md (g) prints the same mid-points and half widths (59.26 +- 0.29; 66.54 +- 0.40) |  | RESULTS.md |  | **PASS** |
| LG-ratio-eqN | tech. report Sec. 4.4 | ratio t_neu/t_pure at equal N (1.116 ... 1.124) |  | computed 1.1162 ... 1.124, quoted 1.116 ... 1.124 |  | **PASS** |
| LG-ratio-full | tech. report Sec. 4.4 | ratio over the full ranges (min neutronised / max pure ... max neutronised / min pure) |  | computed 1.1107 ... 1.1351, quoted 1.11 ... 1.14 |  | **PASS** |
| LG-ratio-pm | tech. report Sec. 4.4 / abstract | ratio 1.12 +- 0.02 covers every equal-N and full-range ratio |  | ratios 1.1107 ... 1.1351 |  | **PASS** |
| LG-AV-rho | tech. report Sec. 4.4 (vacuum edge) | viscous vs inviscid rho_c/rho_c0 - 1 at t = 20, 40, 50 (N = 200, 400, both EOS) |  | 1.1413e-05 <= 2e-5 (+ half unit 5e-06) |  | **PASS** |
| LG-AV-radii | tech. report Sec. 4.4 (vacuum edge) [after edit R2-7] | viscous vs inviscid mass radii (all tabulated mass radii at t = 20, 40, 50; N = 200, 400; both EOS) |  | 1.506e-05 <= 2e-5 (+ half unit 5e-06) |  | **PASS** |
| LG-floor-null | tech. report Sec. 4.4 (vacuum edge) | pressure floors up to 1e-10 p_c change nothing: crossing time (1e-9 relative) and rho_c(30) identical to the baseline (both EOS) |  | [True, True, True, True, True, True] |  | **PASS** |
| LG-floor-blow | tech. report Sec. 4.4 (vacuum edge) | floors 1e-8 and 1e-6 p_c destroy the run (rho_c(30) is NaN: overflow, edge flung out) with first-step kicks of 9-11 and ~1e3 g dt |  | ([8.845275316501786, 11.389955857035089], [982.2662254177402, 1236.3097776858708]) |  | **PASS** |
| LG-kick-N400 | tech. report Sec. 4.4 (vacuum edge) | intrinsic first-step kick of an unfloored surface zone / (g dt), N = 400 (test d) | 0.57 | 0.565355 | 0.005 | **PASS** |
| LG-kick-allN | tech. report Sec. 4.4 (vacuum edge) | same kick for every N: 0.56-0.57 g dt for N = 100 ... 1600, both EOS (produced by lagrangian_L_extra_checks.py; test_d_edge.py itself runs N = 400 only) |  | [0.5611, 0.5639, 0.5654, 0.566, 0.5663, 0.5608, 0.5637, 0.5652, 0.5659, 0.5662] |  | **PASS** |
| LG-lastzone | tech. report Sec. 4.4 (vacuum edge) | last-zone mass x 0.25 ... 4: pure crossing between 59.2 and 59.3 (baseline 59.48); no crossing of the neutronised run before t = 60 |  | pure ['59.314', '59.261', '59.211', '59.220'], baseline 59.481, neutronised ['tmax', 'tmax', 'tmax', 'tmax'] |  | **PASS** |
| LG-pext-pure | tech. report Sec. 4.4 (vacuum edge) | exterior pressure 1e-8 p_c: crossing time, pure [L/c] | 1.99 | 1.9878 | 0.005 | **PASS** |
| LG-pext-neu | tech. report Sec. 4.4 (vacuum edge) | ... neutronised | 0.61 | 0.613372 | 0.005 | **PASS** |
| LG-T1-neu-400-1.5 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 1.5x [L/c] | 63.756 | 63.7557 | 0.0005 | **PASS** |
| LG-T1-neu-400-2.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 2.0x [L/c] | 76.24 | 76.2403 | 0.0005 | **PASS** |
| LG-T1-neu-400-3.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 3.0x [L/c] | 86.071 | 86.0712 | 0.0005 | **PASS** |
| LG-T1-neu-400-5.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 5.0x [L/c] | 92.905 | 92.9045 | 0.0005 | **PASS** |
| LG-T1-neu-400-10.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 10.0x [L/c] | 97.786 | 97.7863 | 0.0005 | **PASS** |
| LG-T1-neu-400-100.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 100.0x [L/c] | 102.656 | 102.656 | 0.0005 | **PASS** |
| LG-T1-neu-400-300.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 400: first time rho_c/rho_c0 reaches 300.0x [L/c] | 103.218 | 103.218 | 0.0005 | **PASS** |
| LG-T1-neu-800-1.5 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 1.5x [L/c] | 63.747 | 63.7472 | 0.0005 | **PASS** |
| LG-T1-neu-800-2.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 2.0x [L/c] | 76.221 | 76.2213 | 0.0005 | **PASS** |
| LG-T1-neu-800-3.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 3.0x [L/c] | 86.047 | 86.0469 | 0.0005 | **PASS** |
| LG-T1-neu-800-5.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 5.0x [L/c] | 92.876 | 92.8762 | 0.0005 | **PASS** |
| LG-T1-neu-800-10.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 10.0x [L/c] | 97.756 | 97.7563 | 0.0005 | **PASS** |
| LG-T1-neu-800-100.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 100.0x [L/c] | 102.627 | 102.627 | 0.0005 | **PASS** |
| LG-T1-neu-800-300.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 800: first time rho_c/rho_c0 reaches 300.0x [L/c] | 103.186 | 103.186 | 0.0005 | **PASS** |
| LG-T1-neu-1600-1.5 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 1.5x [L/c] | 63.746 | 63.7457 | 0.0005 | **PASS** |
| LG-T1-neu-1600-2.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 2.0x [L/c] | 76.217 | 76.217 | 0.0005 | **PASS** |
| LG-T1-neu-1600-3.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 3.0x [L/c] | 86.041 | 86.0411 | 0.0005 | **PASS** |
| LG-T1-neu-1600-5.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 5.0x [L/c] | 92.869 | 92.8693 | 0.0005 | **PASS** |
| LG-T1-neu-1600-10.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 10.0x [L/c] | 97.749 | 97.7489 | 0.0005 | **PASS** |
| LG-T1-neu-1600-100.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 100.0x [L/c] | 102.62 | 102.62 | 0.0005 | **PASS** |
| LG-T1-neu-1600-300.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, adaptive step, N = 1600: first time rho_c/rho_c0 reaches 300.0x [L/c] | 103.178 | 103.178 | 0.0005 | **PASS** |
| LG-T1-neu-R-1.5 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 1.5x [L/c] | 63.745 | 63.7452 | 0.0005 | **PASS** |
| LG-T1-neu-R-2.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 2.0x [L/c] | 76.216 | 76.2156 | 0.0005 | **PASS** |
| LG-T1-neu-R-3.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 3.0x [L/c] | 86.039 | 86.0391 | 0.0005 | **PASS** |
| LG-T1-neu-R-5.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 5.0x [L/c] | 92.867 | 92.867 | 0.0005 | **PASS** |
| LG-T1-neu-R-10.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 10.0x [L/c] | 97.746 | 97.7465 | 0.0005 | **PASS** |
| LG-T1-neu-R-100.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 100.0x [L/c] | 102.618 | 102.618 | 0.0005 | **PASS** |
| LG-T1-neu-R-300.0 | tech. report Table 4 (tab:lag-collapse) | neutronised, Richardson limit (p = 2) from N = 400, 800, 1600, 300.0x [L/c] | 103.176 | 103.176 | 0.0005 | **PASS** |
| LG-T1-pure-200-1.5 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 200: first time rho_c/rho_c0 reaches 1.5x [L/c] | 72.485 | 72.4854 | 0.0005 | **PASS** |
| LG-T1-pure-200-2.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 200: first time rho_c/rho_c0 reaches 2.0x [L/c] | 95.896 | 95.8962 | 0.0005 | **PASS** |
| LG-T1-pure-200-3.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 200: first time rho_c/rho_c0 reaches 3.0x [L/c] | 124.668 | 124.668 | 0.0005 | **PASS** |
| LG-T1-pure-200-5.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 200: first time rho_c/rho_c0 reaches 5.0x [L/c] | 168.759 | 168.759 | 0.0005 | **PASS** |
| LG-T1-pure-400-1.5 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 400: first time rho_c/rho_c0 reaches 1.5x [L/c] | 72.43 | 72.4296 | 0.0005 | **PASS** |
| LG-T1-pure-400-2.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 400: first time rho_c/rho_c0 reaches 2.0x [L/c] | 95.754 | 95.754 | 0.0005 | **PASS** |
| LG-T1-pure-400-3.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 400: first time rho_c/rho_c0 reaches 3.0x [L/c] | 124.396 | 124.396 | 0.0005 | **PASS** |
| LG-T1-pure-400-5.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 400: first time rho_c/rho_c0 reaches 5.0x [L/c] | 167.494 | 167.494 | 0.0005 | **PASS** |
| LG-T1-pure-800-1.5 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 800: first time rho_c/rho_c0 reaches 1.5x [L/c] | 72.419 | 72.4189 | 0.0005 | **PASS** |
| LG-T1-pure-800-2.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 800: first time rho_c/rho_c0 reaches 2.0x [L/c] | 95.722 | 95.7217 | 0.0005 | **PASS** |
| LG-T1-pure-800-3.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 800: first time rho_c/rho_c0 reaches 3.0x [L/c] | 124.328 | 124.328 | 0.0005 | **PASS** |
| LG-T1-pure-800-5.0 | tech. report Table 4 (tab:lag-collapse) | pure, fixed step, N = 800: first time rho_c/rho_c0 reaches 5.0x [L/c] | 167.217 | 167.217 | 0.0005 | **PASS** |
| LG-T1-pure-none | tech. report Table 4 (tab:lag-collapse) | the pure star never reaches 10x, 100x, 300x (N = 200, 400, 800) |  | hits None |  | **PASS** |
| LG-T1-pure-adaptive | tech. report Table 4 caption | pure: fixed step identical to 1e-3 with the adaptive step (N = 400) |  | [(72.4296, 72.4297), (95.754, 95.7541), (124.3961, 124.3961), (167.4937, 167.4937)] |  | **PASS** |
| LG-Rich-neu-2.0 | tech. report Sec. 4.4 (viscous continuation) | neutronised: Richardson limit of the time to 2.0x (RESULTS.md (h)) [L/c] | 76.216 | 76.216 | 0.0005 | **PASS** |
| LG-Rich-neu-5.0 | tech. report Sec. 4.4 (viscous continuation) | neutronised: Richardson limit of the time to 5.0x (RESULTS.md (h)) [L/c] | 92.867 | 92.867 | 0.0005 | **PASS** |
| LG-Rich-neu-10.0 | tech. report Sec. 4.4 (viscous continuation) | neutronised: Richardson limit of the time to 10.0x (RESULTS.md (h)) [L/c] | 97.746 | 97.746 | 0.0005 | **PASS** |
| LG-Rich-neu-100.0 | tech. report Sec. 4.4 (viscous continuation) | neutronised: Richardson limit of the time to 100.0x (RESULTS.md (h)) [L/c] | 102.618 | 102.618 | 0.0005 | **PASS** |
| LG-Rich-neu-300.0 | tech. report Sec. 4.4 (viscous continuation) | neutronised: Richardson limit of the time to 300.0x (RESULTS.md (h)) [L/c] | 103.176 | 103.176 | 0.0005 | **PASS** |
| LG-Rich-pure-2.0 | tech. report Sec. 4.4 (viscous continuation) | pure: Richardson limit of the time to 2.0x (RESULTS.md (h)) [L/c] | 95.711 | 95.711 | 0.0005 | **PASS** |
| LG-Rich-pure-5.0 | tech. report Sec. 4.4 (viscous continuation) | pure: Richardson limit of the time to 5.0x (RESULTS.md (h)) [L/c] | 167.12 | 167.124 | 0.005 | **PASS** |
| LG-tables-1e16 | tech. report Sec. 4.4 (viscous continuation) | adaptive runs stop at 300 rho_c0 = 1.5e15 kg/m^3; the EOS tables end at 1e16 kg/m^3 (RESULTS.md) |  | RESULTS.md (h) header |  | **PASS** |
| LG-300x-s | tech. report Sec. 4.4 (viscous continuation) | time of the 300x event in seconds (103.18 L/c x L/c of the neutronised star) [s] | 0.556 | 0.556088 | 0.0005 | **PASS** |
| LG-fixed-dt1 | tech. report Sec. 4.4 (viscous continuation) | fixed step neu|AV|N400: spurious core crossing time [L/c] | 98.1 | 98.1488 | 0.05 | **PASS** |
| LG-fixed-dt0.5 | tech. report Sec. 4.4 (viscous continuation) | fixed step neu|AV|dt0.5|N400: spurious core crossing time [L/c] | 102.3 | 102.275 | 0.05 | **PASS** |
| LG-fixed-dt0.25 | tech. report Sec. 4.4 (viscous continuation) | fixed step neu|AV|dt0.25|N400: spurious core crossing time [L/c] | 103.4 | 103.445 | 0.05 | **PASS** |
| LG-fixed-zone | tech. report Sec. 4.4 (viscous continuation) [after edit R2-6] | zone of the minimum gap at the core crossing of the three fixed-step runs: 3-5; alternate zones out to 13 have R' <= 0 (lagrangian_L_extra_checks.py) |  | min-gap zones [4, 5, 3] (json); zones with R'<=0 up to 13 |  | **PASS** |
| LG-pure-max | tech. report Sec. 4.4 / Table 4 caption / abstract | pure star: maximum of rho_c/rho_c0 (N = 800) | 5.44 | 5.44174 | 0.005 | **PASS** |
| LG-pure-tmax | tech. report Sec. 4.4 / Table 4 caption | ... reached at t [L/c] | 189 | 189 | 0.5 | **PASS** |
| LG-pure-300 | tech. report Sec. 4.4 | ... falls back to +63 % at t = 300 | 63 | 62.6024 | 0.5 | **PASS** |
| LG-pure-max-abs | tech. report abstract | maximum "5.4x" of the pure star | 5.4 | 5.44174 | 0.05 | **PASS** |
| LG-pure-halfT | tech. report Sec. 4.4 | half the linear period of the pure star, pi/omega in units of L/c (omega^2 = 9.273 s^-2) [L/c] | 193 | 193.078 | 0.5 | **PASS** |
| LG-pure-halfT-agree | tech. report Sec. 4.4 | the time of the maximum (189) is within 3 % of half the linear period (193.1): "agrees" |  | t_max 189.0 vs 193.1 |  | **PASS** |
| LG-pure-conv | tech. report Sec. 4.4 | the maximum of rho_c is converged to 0.3 % (N = 400 -> 800, relative change [%]) |  | 0.19952 <= 0.3 (+ half unit 0.05) |  | **PASS** |
| LG-pure-nocross | tech. report Sec. 4.4 | the pure star has no crossing to t = 300 (viscous runs, N = 200, 400, 800) |  | status tmax, tend >= 299.9 |  | **PASS** |
| LG-same-comp-1.5 | tech. report Sec. 4.4 / summary | ratio t_neu/t_pure of the times to reach 1.5x (neutronised adaptive N = 1600 / pure N = 800) | 0.88 | 0.880235 | 0.005 | **PASS** |
| LG-same-comp-2.0 | tech. report Sec. 4.4 / summary | ratio t_neu/t_pure of the times to reach 2.0x (neutronised adaptive N = 1600 / pure N = 800) | 0.8 | 0.796235 | 0.005 | **PASS** |
| LG-same-comp-3.0 | tech. report Sec. 4.4 / summary | ratio t_neu/t_pure of the times to reach 3.0x (neutronised adaptive N = 1600 / pure N = 800) | 0.69 | 0.69205 | 0.005 | **PASS** |
| LG-same-comp-5.0 | tech. report Sec. 4.4 / summary | ratio t_neu/t_pure of the times to reach 5.0x (neutronised adaptive N = 1600 / pure N = 800) | 0.56 | 0.555383 | 0.005 | **PASS** |
| LG-same-comp-R | tech. report Sec. 4.4 / summary | RESULTS.md (h) prints the ratios 0.880, 0.796, 0.692, 0.555 |  | RESULTS.md |  | **PASS** |
| LG-dlnrho-50 | tech. report Sec. 4.4 | d ln rho_c / dt over t = 50-60 (neutronised, viscous, N = 1600) [c/L] | 0.012 | 0.0121808 | 0.0005 | **PASS** |
| LG-dlnrho-80 | tech. report Sec. 4.4 | ... over t = 80-90 | 0.053 | 0.0534777 | 0.0005 | **PASS** |
| LG-sigma-lin | tech. report Sec. 4.4 | linear growth rate sigma [c/L] | 0.0242 | 0.0241894 | 5e-05 | **PASS** |
| LG-rep-rhs | tech. report Sec. 4.5 | replica of lagrangian_envelope.py: right-hand side agrees to (largest relative difference, both EOS) | 3e-10 | 3.01021e-10 | 5e-11 | **PASS** |
| LG-rep-rk | tech. report Sec. 4.5 | ... 300 Runge-Kutta steps agree to | 1e-13 | 1.15214e-13 | 5e-14 | **PASS** |
| LG-rep-pure-t | tech. report Sec. 4.5 / Established | reproduced crossing of the pure re-creation: t [L/c] | 14.783 | 14.783 | 0.0005 | **PASS** |
| LG-rep-pure-shell | tech. report Sec. 4.5 | ... shell | 0 | 0 | 0 | **PASS** |
| LG-rep-pure-rR | tech. report Sec. 4.5 | ... r/R of that shell | 0.01 | 0.01 | 0.005 | **PASS** |
| LG-rep-neu-t | tech. report Sec. 4.5 / Established | reproduced crossing of the neutronised re-creation: t [L/c] | 2.188 | 2.188 | 0.0005 | **PASS** |
| LG-rep-neu-shell | tech. report Sec. 4.5 | ... shell (of 300 shells, labels 0 ... 299) | 299 | 299 | 0 | **PASS** |
| LG-rep-acc0-pure | tech. report Sec. 4.5 | acceleration of node 0 in the initial data of the re-creation, pure [g] | -0.77 | -0.772128 | 0.005 | **PASS** |
| LG-rep-acc0-neu | tech. report Sec. 4.5 | ... neutronised [g] | -0.35 | -0.349004 | 0.005 | **PASS** |
| LG-rep-accN | tech. report Sec. 4.5 | acceleration of the outer node is -1.0 g in both cases (pressure force switched off by the weight) |  | [-1.0, -1.0] |  | **PASS** |
| LG-rep-dr | tech. report Sec. 4.5 | spacing of the outermost pair | 0.0033 | 0.00331 | 5e-05 | **PASS** |
| LG-rep-ff-pure | tech. report Sec. 4.5 | free-fall time sqrt(2 dr/g) of the outermost pair, pure [L/c] | 2.25 | 2.24868 | 0.005 | **PASS** |
| LG-rep-ff-neu | tech. report Sec. 4.5 | ... neutronised | 2.46 | 2.45949 | 0.005 | **PASS** |
| LG-rep-N150-neu | tech. report Sec. 4.5 | N = 150 nodes, neu: crossing time [L/c] | 2.73 | 2.734 | 0.005 | **PASS** |
| LG-rep-N300-neu | tech. report Sec. 4.5 | N = 300 nodes, neu: crossing time [L/c] | 2.19 | 2.188 | 0.005 | **PASS** |
| LG-rep-N600-neu | tech. report Sec. 4.5 | N = 600 nodes, neu: crossing time [L/c] | 2.23 | 2.228 | 0.005 | **PASS** |
| LG-rep-N150-pure | tech. report Sec. 4.5 | N = 150 nodes, pure: crossing time [L/c] | 2.48 | 2.482 | 0.005 | **PASS** |
| LG-rep-N300-pure | tech. report Sec. 4.5 | N = 300 nodes, pure: crossing time [L/c] | 14.78 | 14.783 | 0.005 | **PASS** |
| LG-rep-N600-pure | tech. report Sec. 4.5 | N = 600 nodes, pure: crossing time [L/c] | 1.9 | 1.8985 | 0.005 | **PASS** |
| LG-rep-shells | tech. report Sec. 4.5 | the crossing is in shell N-1 for every N except the pure N = 300 run, which crosses at the centre (shell 0) |  | [True, True, True, True, True, True] |  | **PASS** |
| LG-rep-nomask | tech. report Sec. 4.5 | without the weight the run fails at the first step (t = 0.001, shell 0, both EOS) |  | [('neu', 0.001), ('pure', 0.001)] |  | **PASS** |
| LG-rep-amp | tech. report Sec. 4.5 | central-density oscillation of the pure re-creation: maximum |rho_c/rho_c0 - 1| | 0.2 | 0.20543 | 0.01 | **PASS** |
| LG-rep-erratic | tech. report Sec. 4.5 | rho_c/rho_c0 of the pure re-creation at t = 7 ... 14.5 lies between 0.9 and 1.2 and is not monotone |  | 0.897 ... 1.205 |  | **PASS** |
| LG-rep-conv10 | tech. report Sec. 4.5 | converged central density at t = 10: rises by [%] | 1.8 | 1.75566 | 0.05 | **PASS** |
| LG-rep-conv20 | tech. report Sec. 4.5 | ... at t = 20 [%] (Richardson limit of pure N = 400, 800, 1600) | 4.7 | 4.73533 | 0.05 | **PASS** |
| LG-rep-mono | tech. report Sec. 4.5 | the converged central density (inviscid, N = 1600) rises monotonically up to t = 50 |  | [0.0, 0.0175, 0.0474, 0.0879, 0.1437, 0.2202] |  | **PASS** |
| LG-eul-ms-neu | tech. report Sec. 5.2 | crossing of the neutronised re-creation in ms ("approximately 12-80 ms") | 12 | 11.79 | 0.5 | **PASS** |
| LG-eul-ms-pure | tech. report Sec. 5.2 | crossing of the pure re-creation in ms ("approximately 12-80 ms") | 80 | 78.99 | 1.5 | **PASS** |
| LG-claim-ratio | tech. report abstract / Not established | ratio 9.1/14.2 of the originally quoted crossing times | 0.64 | 0.640845 | 0.005 | **PASS** |
| LG-claim-pct | tech. report Not established | "35 % sooner" = 100 (1 - 9.1/14.2) [%] | 35 | 35.9155 | 1 | **PASS** |
| LG-claim-f4 | tech. report summary / Not established | converged pure first-loss time over the originally quoted 14.2 ("factor 4-7") | 4 | 4.17295 | 0.5 | **PASS** |
| LG-claim-f7 | tech. report summary / Not established | converged neutronised first-loss time over the originally quoted 9.1 | 7 | 7.31161 | 0.5 | **PASS** |
| LG-abs-300x | tech. report abstract | neutronised central density reaches 300x at t = [L/c] | 103.2 | 103.176 | 0.05 | **PASS** |
| LG-abs-59 | tech. report abstract | pure first loss of ordering t approx [L/c] | 59 | 59.2559 | 0.5 | **PASS** |
| LG-abs-66 | tech. report abstract / summary | neutronised first loss of ordering t approx [L/c] | 66 | 66.5357 | 1 | **PASS** |
| LG-est-spectrum | tech. report Established | linear spectrum agrees with the independent solver to (largest relative difference of the fundamental mode, all four cases) | 3e-05 | 3.41826e-05 | 5e-06 | **PASS** |
| LG-log-pass | tech. report Sec. 4 (listed tests) | logs of lagrangian_L: 23 PASS lines (verify_ms_equations 1, test_a 3, test_b 12, repro_spline 2, test_i 4, test_j 1) and no FAIL line |  | 23 PASS, 0 FAIL [] |  | **PASS** |
| LG-log-exit | tech. report Sec. 4 | timing.log: all 17 steps of run_all.sh exited with 0 |  | 17 steps |  | **PASS** |
| LG-dust-hom-err | tech. report Sec. 4.3 (e) | homogeneous ball: max relative error of R(t) at N = 400 |  | 1.7786e-07 <= 1.8e-7 (+ half unit 5e-09) |  | **PASS** |
| LG-dust-inh-err | tech. report Sec. 4.3 (e) | inhomogeneous ball: max relative error of R(t) at N = 400 |  | 1.4408e-07 <= 1.4e-7 (+ half unit 5e-09) |  | **PASS** |
| LG-dust-order | tech. report Sec. 4.3 (e) | observed order 2.00 in R and U (all grid pairs N = 25 ... 400, both cases; within 0.01) |  | orders 1.9946 ... 2.0000 |  | **PASS** |
| LG-dust-m-hom | tech. report Sec. 4.3 (e) | m_j constant to (homogeneous, N = 400) |  | 1.415e-07 <= 1.4e-7 (+ half unit 5e-09) |  | **PASS** |
| LG-dust-iter3 | tech. report Sec. 4.3 (e) | with three iterations: error 3.0e-6 at 1600 steps against 1.8e-7 converged; order 1.1 |  | fixed3 2.97e-06, converged 1.78e-07, order 1.08 |  | **PASS** |
| LG-dg2-pure | tech. report Sec. 4.1 | max |Gamma^2(t) - Gamma^2(0)| up to t = 50, pure | 0.00027 | 0.000274522 | 5e-06 | **PASS** |
| LG-dg2-neu | tech. report Sec. 4.1 | same, neutronised | 0.00025 | 0.000246642 | 5e-06 | **PASS** |
| LG-dg2-U | tech. report Sec. 4.1 | largest |U|/c in those runs | 0.004 | 0.00373408 | 0.0005 | **PASS** |
| LG-dg2-2mR-pure | tech. report Sec. 4.1 | largest 2m/R, pure | 0.0045 | 0.00446856 | 5e-05 | **PASS** |
| LG-dg2-2mR-neu | tech. report Sec. 4.1 | largest 2m/R, neutronised | 0.0039 | 0.00386136 | 5e-05 | **PASS** |
| LG-dg2-frac | tech. report Sec. 4.1 | the variation is about 6 % of the largest 2m/R (5.5 % ... 6.6 %) |  | fractions 0.061, 0.064 |  | **PASS** |
| LG-node-counts | tech. report Sec. 4.5 | the weight switches the pressure force off for 3 (pure) / 4 (neutronised) outer nodes (repro_spline.log) |  | see logs/repro_spline.log |  | **PASS** |
| LG-ref-sums | tech. report Sec. 4 (unmodified copies of the report routines) | ref/: sha256 of the 11 copies equal SHA256SUMS.copy, and each copy is byte-identical to the script of the same name in this repository (ancillary/ or ancillary/linear_perturbation/) |  | 11 files; differing: [] |  | **PASS** |
| LG-ref-log | tech. report Sec. 4 | logs/ref_checksums.log: 11 lines, all ": OK" |  | 11 lines |  | **PASS** |

## Generators run in this invocation

| step | state | wall [s] | CPU [s] | estimate (one laptop core) | note |
|---|---|---|---|---|---|
| sympy_eos_checks | OK (last recorded, not run now) | 15.0 | 2.0 | 3 s |  |
| sympy_weak_form | OK (last recorded, not run now) | 7.5 | 1.0 | 10 s |  |
| sympy_fold_scaling | OK (last recorded, not run now) | 25.1 | 3.4 | 10 s |  |
| verify_exact_poisson | OK (last recorded, not run now) | 8.3 | 1.1 | 1 s |  |
| analytic_check2 | OK (last recorded, not run now) | 12.2 | 1.7 | 10 s |  |
| analytic_check_uniform | OK (last recorded, not run now) | 13.0 | 1.8 | 10 s |  |
| verify_pure_mode | OK (last recorded, not run now) | 7.4 | 1.0 | 1 s |  |
| background | OK (last recorded, not run now) | 6.7 | 0.9 | 1 s |  |
| eos_numbers | OK (last recorded, not run now) | 7.5 | 1.0 | 2 s |  |
| tov_solve | OK (last recorded, not run now) | 8.4 | 2.0 | 2 s |  |
| misc_checks | OK (last recorded, not run now) | 153.0 | 18.3 | 16 s |  |
| gr_rayleigh_ritz | OK (last recorded, not run now) | 10.4 | 2.9 | 2 s |  |
| rr_at | OK (last recorded, not run now) | 15.3 | 4.3 | 5 s |  |
| weak_strong_identity | OK (last recorded, not run now) | 17.6 | 4.4 | 5 s |  |
| validate_galerkin_sg | OK (last recorded, not run now) | 96.0 | 24.1 | 5 s |  |
| fig1_pekeris | OK (last recorded, not run now) | 111.5 | 30.7 | 22 s |  |
| figs | OK (last recorded, not run now) | 1767.1 | 153.0 | 2.2 min |  |
| figs_friedman | OK (last recorded, not run now) | 11.0 | 3.0 | 3 s |  |
| turning_point_first | OK (last recorded, not run now) | 1864.6 | 161.4 | 3 min | regression test of the first implementation; the shipped log code/turning_point_table.log is the output of the three sections run together (--full) |
| turning_point_table | OK (last recorded, not run now) | 551.8 | 345.2 | 10-25 min |  |
| ap_check_eos | OK (last recorded, not run now) | 2.2 | 0.6 | 1 s |  |
| ap_table | OK (last recorded, not run now) | 2335.1 | 205.0 | 1.5 min |  |
| ap_conv1 | OK (last recorded, not run now) | 172.2 | 39.6 | 1 min |  |
| ap_nwcmp | OK (last recorded, not run now) | 1472.8 | 127.6 | 1 min |  |
| ap_mscan | OK (last recorded, not run now) | 40.1 | 8.7 | 20 s |  |
| ap_myRR | OK (last recorded, not run now) | 41.7 | 8.6 | 20 s |  |
| ap_fem | OK (last recorded, not run now) | 142.5 | 27.4 | 1 min |  |
| ap_resid | OK (last recorded, not run now) | 9.5 | 1.8 | 10 s |  |
| ap_run_shoot | OK (last recorded, not run now) | 89.1 | 17.6 | 20 s |  |
| ap_projdiag | OK (last recorded, not run now) | 12.3 | 2.4 | 5 s |  |
| ap_projres | OK (last recorded, not run now) | 39.0 | 7.6 | 1 min |  |
| ap_projdens | OK (last recorded, not run now) | 94.3 | 18.7 | 2 min |  |
| homologous_calibrated | OK (last recorded, not run now) | 14.5 | 2.8 | 4 s |  |
| exact_energy_model | OK (last recorded, not run now) | 58.4 | 11.5 | 10 s |  |
| integrate_collapse | OK (last recorded, not run now) | 810.2 | 77.7 | 2 min |  |
| final_single | OK (last recorded, not run now) | 4999.9 | 402.8 | 6 min |  |
| fig8_collapse_timescale | OK (last recorded, not run now) | 255.1 | 51.9 | 46 s | reads the shipped final_single.json |
| twomode_final | OK (last recorded, not run now) | 5479.3 | 435.0 | 5 min |  |
| integrate_twomode | OK (last recorded, not run now) | 5530.9 | 437.2 | 5 min |  |
| survey_spot | OK (last recorded, not run now) | 507.6 | 109.7 | 1-2 min | two rows of the composition survey recomputed (ideal gas at offset 2.0; Fe56 at the fiducial density) and compared with survey4.json |
| survey4 | OK (last recorded, not run now) | 6512.8 | 1435.0 | 16 min |  |
| fold_step1 | OK (last recorded, not run now) | 183.1 | 43.4 | 37 s |  |
| fold_step2 | OK (last recorded, not run now) | 23.1 | 2.7 | 2 s |  |
| fold_dyn | OK (last recorded, not run now) | 243.2 | 42.5 | 34 s |  |
| fold_dyn2 | OK (last recorded, not run now) | 124.7 | 98.9 | 1.7 min | the original script of the archive (restored); the shipped fold_dyn2.json is kept and compared with the regenerated copy (checks FT-json-regen, PI-old-estimator) |
| pi_tritronquee | OK (last recorded, not run now) | 66.8 | 7.9 | 7 s |  |
| pi_checks | OK (last recorded, not run now) | 59.1 | 7.0 | 30 s |  |
| figs_fold | OK (last recorded, not run now) | 105.2 | 12.6 | 11 s |  |
| damped_pi | OK (last recorded, not run now) | 604.3 | 72.0 | 57 s |  |
| damp_phys | OK (last recorded, not run now) | 22.1 | 2.6 | 2 s |  |
| fold_numbers | OK (last recorded, not run now) | 1.5 | 0.2 | 1 s |  |
| fold_comp | OK (last recorded, not run now) | 1157.1 | 171.4 | 2.3 min |  |
| salp_crit | OK (last recorded, not run now) | 512.2 | 61.1 | 54 s |  |
| lag_eigs_quick | OK (last recorded, not run now) | 207.8 | 28.2 | 1 min | N = 100, 200 only (the shipped lag_eigs.log has N = 100 ... 800) |
| lag_eigs | OK (last recorded, not run now) | 1224.4 | 270.5 | 5 min |  |
| lag_jac | OK (last recorded, not run now) | 32.8 | 4.0 | 4 s |  |
| lag_test1 | OK (last recorded, not run now) | 154.7 | 15.4 | 12 s |  |
| lag_checks_resid | OK (last recorded, not run now) | 107.0 | 12.8 | 2 min |  |
| lag_checks_stages | OK (last recorded, not run now) | 613.9 | 73.3 | 20 s |  |
| lag_checks_shells_N200 | OK (last recorded, not run now) | 518.5 | 61.9 | 1.5 min |  |
| lag_checks_velnorm_N200_base | OK (last recorded, not run now) | 187.4 | 22.4 | 1 min |  |
| lag_checks_velnorm_N200_div | OK (last recorded, not run now) | 186.9 | 22.3 | 1 min |  |
| lag_checks_velnorm_N200_mul | OK (last recorded, not run now) | 181.7 | 21.7 | 1 min |  |
| fold_inhom_200 | OK (last recorded, not run now) | 535.9 | 63.9 | 51 s |  |
| fold_inhom_100 | OK (last recorded, not run now) | 520.4 | 62.5 | 15 s |  |
| fold_inhom_400 | OK (last recorded, not run now) | 612.5 | 107.3 | 2 min |  |
| fold_inhom_800 | OK (last recorded, not run now) | 2419.4 | 219.0 | 4 min |  |
| fold_inhom_1600 | OK (last recorded, not run now) | 3022.9 | 543.6 | 10 min |  |
| fold_to_BH_200 | OK (last recorded, not run now) | 385.7 | 46.0 | 36 s |  |
| fold_to_BH_400 | OK (last recorded, not run now) | 546.3 | 94.2 | 1.5 min |  |
| fold_to_BH_800 | OK (last recorded, not run now) | 2571.3 | 225.5 | 4 min |  |
| fold_to_BH_1600 | OK (last recorded, not run now) | 3116.3 | 585.7 | 12 min |  |
| fold_trapped_200 | OK (last recorded, not run now) | 426.5 | 50.4 | 1 min |  |
| fold_trapped_200_cfl0.2 | OK (last recorded, not run now) | 582.5 | 102.1 | 2 min |  |
| fold_trapped_200_cfl1.2 | OK (last recorded, not run now) | 103.9 | 22.2 | 30 s |  |
| fold_trapped_400 | OK (last recorded, not run now) | 349.4 | 96.1 | 2 min |  |
| fold_trapped_800 | OK (last recorded, not run now) | 2410.6 | 257.8 | 5 min |  |
| fold_trapped_1600 | OK (last recorded, not run now) | 4148.3 | 680.3 | 12 min |  |
| lag_checks_shells_N400 | OK (last recorded, not run now) | 360.6 | 104.8 | 4 min |  |
| tr_eos | OK (last recorded, not run now) | 6.2 | 0.8 | 1 s |  |
| tr_initial_profile | OK (last recorded, not run now) | 167.4 | 20.1 | 20 s |  |
| tr_neutronization | OK (last recorded, not run now) | 6.6 | 0.8 | 2 s |  |
| tr_eos_neutronized | OK (last recorded, not run now) | 7.9 | 1.0 | 5 s |  |
| tr_initial_profile_neu | OK (last recorded, not run now) | 100.6 | 12.0 | 20 s |  |
| tr_collapse2 | OK (last recorded, not run now) | 38.6 | 4.6 | 10 s |  |
| tr_evolve2 | OK (last recorded, not run now) | 7.7 | 0.9 | 10 s | right-hand-side demonstration only (pure-EOS version of the Lagrangian code; evolve_neu.py was not retained) |
| tr_apply_to_real_star | OK (last recorded, not run now) | 40.4 | 5.4 | 20 s |  |
| tr_sanity_check_real | OK (last recorded, not run now) | 27.4 | 3.4 | 20 s |  |
| tr_stable_crosscheck | OK (last recorded, not run now) | 58.6 | 6.9 | 20 s |  |
| tr_eigenfunction_shape | OK (last recorded, not run now) | 48.0 | 6.4 | 20 s |  |
| tr_euler_convergence | OK (last recorded, not run now) | 821.5 | 47.6 | 1 min CPU |  |
| tr_control_test | OK (last recorded, not run now) | 1301.8 | 74.9 | 3-4 min CPU |  |
| tr_study1 | OK (last recorded, not run now) | 766.4 | 44.6 | 1 min CPU |  |
| tr_study2 | OK (last recorded, not run now) | 3577.8 | 397.6 | 15-20 min CPU |  |
| tr_lagrangian_envelope_pure | OK (last recorded, not run now) | 2155.3 | 124.0 | 2-3 min CPU | attempted re-creation of the lost shell-crossing driver; does NOT reproduce the report (see its docstring); no check depends on it |
| tr_lagrangian_envelope_neu | OK (last recorded, not run now) | 379.6 | 21.9 | 1 min CPU | as above |
| tr_get_eigenfunction | OK (last recorded, not run now) | 3.2 | 3.2 | 20 s |  |
| tr_check_guard_equivalence_strong800 | OK (last recorded, not run now) | 31.4 | 28.5 | 30 s | guarded and original seeded module give identical logged states up to the stop (eps = 1e-2, N = 800) |
| tr_check_guard_equivalence_1600 | OK (last recorded, not run now) | 253.2 | 169.9 | 6 min CPU | same for eps = 1e-4, N = 1600 (Sec. 7.2 of the report) |
| tr_run_seeded2_guarded | OK (last recorded, not run now) | 80.7 | 68.6 | 4 min CPU |  |
| tr_run_pure_mode | OK (last recorded, not run now) | 1576.0 | 91.0 | ? |  |
| tr_run_seeded_strong_guarded | OK (last recorded, not run now) | 21.1 | 18.1 | 2 min CPU |  |
| tr_lagrangian_L_extra_checks | OK (last recorded, not run now) | 40.5 | 40.2 | 40 s | first-step kick for N = 100 ... 1600, central 1 - a, zones of the fixed-step core crossing (statements of Sec. 4 that no log of lagrangian_L contains) |

Wall times on a shared machine are upper bounds; CPU times are the cost of the step itself.
