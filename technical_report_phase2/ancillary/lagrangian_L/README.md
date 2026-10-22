# Convergence-verified Lagrangian (Misner-Sharp) evolution of the near-critical white-dwarf configuration

Directory: `technical_report_phase2/ancillary/lagrangian_L/`
Purpose: convergence-verified re-do of the Lagrangian (Misner-Sharp) evolution of Section 4 of the technical report. The report's own EOS and profile routines are used only through hash-verified, unmodified copies in `ref/` (`logs/ref_checksums.log`).
## How to run

```
bash run_all.sh        # python3 + numpy + scipy only; regenerates every log, results/*.json and results/RESULTS.md
```
Wall time of the last complete run on the 2 shared cores: **1177 s (19.6 min)** (`logs/wall_time.log`, per-step times in `logs/timing.log`;
two job streams, `nice`d).  Every number quoted below and in Section 4 of the technical report (`technical_report_phase2/main.tex`) is in `results/RESULTS.md` or in a log.

## Layout

| path | content |
|---|---|
| `derivations.md` | the exact system, derivation of the momentum equation, the discrepancies with Sec. 2 of the report and with `evolve2.py` |
| `src/lag.py` | the scheme (staggered, conservative; kick-drift-kick leapfrog; optional von Neumann-Richtmyer viscosity; free surface) |
| `src/eosmodel.py`, `structure.py`, `initdata.py`, `harness.py` | EOS tables built from the report's routines, TOV solver (enthalpy form), initial-data families, drivers (fixed and CFL-adaptive step) |
| `src/linear_gr.py`, `report_galerkin.py` | independent GR Sturm-Liouville FEM solver; the report's Galerkin solver called as in `apply_to_real_star.py` |
| `src/verify_ms_equations.py` | exact 2-jet automatic differentiation of the Einstein tensor (no CAS): identities I1-I6 (field equations, lapse equation, energy equation) + negative controls |
| `src/test_a_static.py` ... `test_j_gamma_size.py` | tests (a) static, (b) linear, (c) convergence, (d) vacuum edge, (e) conservation, (f) distance from equilibrium, (g) T=150, (h) viscous collapse metrics, (i) exact Oppenheimer-Snyder / Lemaitre-Tolman-Bondi dust solutions, (j) size of the time variation of Gamma^2 |
| `src/repro_spline.py` | exact replica of `lagrangian_envelope.py` and its diagnosis |
| `src/make_tables.py` -> `results/RESULTS.md` | the results table |

## What was done

1. **System** (`derivations.md`): metric `-a^2 dt^2 + X^2 dr^2 + R^2 dOmega^2`, `U = a^-1 dR/dt`, `Gamma^2 = 1+U^2-2m/R`; `m' = 4 pi e R^2 R'`, `dm/dt = -4 pi a p R^2 U`,
   `dGamma/dt = a' U/X`, `a'/a = -p'/(e+p)`, `dU/dt = -a[(m+4 pi R^3 p)/R^2 + Gamma^2 p'/((e+p)R')]`.  The first three identities, the algebra of the last, the lapse equation and the energy equation were verified on
   the Einstein tensor / the divergence of the perfect-fluid stress tensor by exact jets (identities I1-I6, residuals <= 1.5e-13 except I6 3.6e-15; negative controls detected at 4e-4 and 3e-3).
2. **Scheme**: fixed baryon content per zone, radii and velocities on interfaces, thermodynamics at zone centres, `m_j = sum e V` exactly, `Gamma` by fixed-point iteration of the coupled rho_N-m-Gamma relation to max|dGamma| < 1e-13 (a fixed count of 3 iterations, used in the first version of this study, leaves an O(dt) error in strongly relativistic motion: see test (i); it did not change any result of the weak-field runs, differences <= 1e-10 in the quoted times),
   leapfrog with CFL 0.4 on `c_s^2 = Gamma_1 p/(e+p)`, free surface.  Barotropic closure `e = rho_report`, `p = P(e)`; `rho_N` from the first law.
3. **Tests** (all PASS lines are in `logs/`, summarised at the end of `run_all_console.log`).

## Findings

### PROVEN (by the tests; numbers are in `results/RESULTS.md`)
* The scheme is second order: pointwise-sampled TOV data deviate from static by 2.2e-4 ... 8.1e-7 (N = 50...800, observed order 2.0-2.1; rho_c order 1.7-1.8); the discretely balanced
  equilibrium stays static to rounding (radii <= 1e-14 L, central density <= 1.2e-11 at t <= 30 and <= 5.7e-11 at t <= 150, the latter including the growth of rounding noise in the unstable neutronised mode).
* **Exact nonlinear solutions** (test (i), `results/test_i_dust.json`): for pressureless dust the scheme reproduces the Oppenheimer-Snyder cycloid (homogeneous ball, 2M/R0 = 0.6) and the Lemaitre-Tolman-Bondi cycloids (inhomogeneous ball, 2M/R0 = 0.5) up to 0.8 of the collapse time with max relative error in R(t) of 1.8e-7 and 1.4e-7 at N = 400 (4 steps per zone) and observed order 2.00 in both R and U, with m_j conserved to 1.4e-7 / 4e-8.  With 3 fixed-point iterations of the Gamma coupling the same runs converge only at order 1.1 (error 3.0e-6): the production setting iterates to 1e-13.
* Size of the Gamma variation (test (j)): up to t = 50, max |Gamma^2(t) - Gamma^2(0)| = 2.7e-4 (pure), 2.5e-4 (neutronised) with |U| <= 0.004 c, i.e. about 6 % of the largest 2m/R (4.5e-3, 3.9e-3).  It is not negligible for a nearly neutrally stable star, but it is not 'comparable to the whole relativistic correction'.
* Linear spectrum of the scheme = independent GR Sturm-Liouville FEM to 2e-5...3.4e-5 (rho_c=5e12) and 6e-6 (1e11); Newtonian limit = the report's Galerkin solver to 2e-5...8e-5; time-domain
  frequency = eigenvalue to 1.3e-6; seeded growth rate of the unstable mode 0.02417-0.02419 c/L vs 0.024189 (<= 1e-3).
* Exact equilibria (rho_c = 5e12): pure Chandrasekhar `omega^2 = +9.273 s^-2` (stable; TOV M = 1.4101 Msun, R = 1599.3 km); neutronised `omega^2 = -20.144 s^-2`, **e-folding 41.3 L/c = 0.223 s**
  (TOV M = 1.1901 Msun, R = 1613.7 km).  The report's `-18.68` is the Newtonian limit (we get -18.721); GR correction +7.6 %.  **Neither exact equilibrium ever crosses** (R' = 1 to rounding, t <= 150).
* The report's Newtonian initial profile is **not** an equilibrium of the Misner-Sharp equations: initial acceleration -0.0046 g (pure), -0.0040 g (neutronised), N-independent; N-independent
  transient (+1.76 %/+1.59 % in rho_c at t=10, +8.8 % at t=30).  Because the star is nearly neutrally stable, the report's central density is 55 % below (pure) / 29 % above (neutronised) that of the nearest GR
  equilibrium with the same zone masses.  Every earlier Lagrangian run therefore evolved a large-amplitude forced response.
* **Time to first R' <= 0 for those data (inviscid, N = 200...1600)**: pure 59.5, 59.5, 59.3, 59.0 L/c; neutronised 66.9, 66.4, 66.4, 66.1 L/c, always in the outermost zones (< 3.1e-6 of the mass
  outside).  Bulk quantities converge at second order; the crossing time does not converge (range 58.97-59.55 and 66.14-66.93; halving dt changes it by <= 0.26 %).  Ratio t_neu/t_pure = 1.12
  (1.116-1.124 at equal N; 1.11-1.14 over the full ranges): the neutronised crossing is **later**, not 35 % sooner.
* With the viscous scheme (and a step limited by the current CFL bound) the neutronised central density runs away: rho_c/rho_c0 = 2, 5, 10, 100, 300 at t = 76.216, 92.867, 97.746, 102.618, 103.176 L/c
  (Richardson limits, observed order 2.0-2.1; N = 400, 800, 1600; independent of dt/2).  The pure star oscillates: maximum 5.44x at t = 189 (half the linear period 193 L/c), back to 1.63x at t = 300,
  no crossing.  Times to equal compression: ratio t_neu/t_pure = 0.88, 0.80, 0.69, 0.56 for 1.5, 2, 3, 5x - it depends on the threshold.
* A fixed time step chosen from the initial CFL bound produces a spurious core "crossing" at t = 98.1, 102.3, 103.4 for dt x 1, 1/2, 1/4 (time-step artefact; not used).
* **The re-creation `lagrangian_envelope.py` is a discretisation artefact**: an exact matrix replica reproduces the original RHS (3e-10) and the quoted crossings to the digit (pure 14.783 in shell 0 = the
  inner boundary of the 1 % hole; neutronised 2.188 in shell 299).  Diagnosis: node accelerations of -0.77 g / -0.35 g (node 0) and -1.0 g (outer node) in the initial data; the outermost pair free-falls its
  spacing in 2.25 / 2.46 L/c; crossing time and shell change with N (2.73/2.19/2.23 and 2.48/14.78/1.90); without the edge weight it fails at step 1; the "+-20 %" rho_c oscillation is erratic grid noise
  (converged rho_c rises monotonically by 1.8 % at t=10).
* The report's "spurious pressure-floor kick" is reproduced: floors <= 1e-10 p_c change nothing; 1e-8 and 1e-6 p_c give first-step kicks of 9-11 and ~1e3 g dt and destroy the run.
* Mass drift <= 2e-11 over 50 L/c; Newtonian energy drift <= 5e-8 |W|; discrete energy-equation defect O(h^2) (orders 2.00-2.12); Gamma definition to 1e-16.
* **Errors in the report/code** (all with evidence): (i) Sec. 2 "R'/X independent of tau" is false (identity I3); (ii) `evolve2.py` uses `dU/dt = -a^2[...]` with `dR/dt = aU` (factor a too many); (iii) `evolve2.py` freezes
  `Gamma` at its initial value; (iv) the Newtonian initial data are not an equilibrium; (v) the stable control (`stable_crosscheck.py`, rho_c = 1e11) uses `rmax = 3e6` m, i.e. a truncated star (omega^2 = 4.65 quoted, 1.756 for the full star).

### CONSISTENT (numerical evidence only)
* The edge crossing is a vacuum-edge shock (sound speed -> 0): it is removed by viscosity while the bulk changes by <= 2e-5.
* The collapse of the neutronised star is a nonlinear large-amplitude process (d ln rho_c/dt rises from 0.012 to 0.053 vs the linear 0.0242).

### NOT ESTABLISHED
* Any crossing time as a property of the star: it is set by the initial perturbation, and the edge crossing time does not converge.  The original 14.2 / 9.1 (lost driver) are not reproduced by anything converged
  (the converged first R' <= 0 is at ~59 and ~66, i.e. 4-7 times later) and their origin is unknown.
* The ratio t_neu/t_pure: it is 1.12 for the criterion of the claim (first R' <= 0), 0.56-0.88 for equal central compression, and undefined for core collapse (the pure star does not collapse before t = 300).
* Pure-star behaviour beyond t = 300; evolution past rho_c ~ 1.5e15 kg/m^3 (the neutronisation EOS is not valid there); dependence on the viscosity coefficients beyond the two tested; non-radial or non-barotropic dynamics.
* Convergence of anything involving the vacuum edge at the 1e-5 mass level.

## Verdict on "35 % sooner"
Not supported.  For the claim as stated (first shell crossing, R' <= 0, outer shell) the converged scheme gives ratio **1.12 +- 0.02** (neutronised later), and the crossing is a grid-scale vacuum-edge event whose time does not converge.
The only robust effect of neutronisation is physical and different: the exact neutronised equilibrium is linearly unstable (41.3 L/c) and the pure one is not, and with viscosity only the neutronised core collapses.
The quantity "35 % sooner" is not a well-defined property of the star: with equal-compression criteria the ratio ranges over 0.56-0.88 (the quoted 0.64 lies inside this range, but for no distinguished threshold).

## Limitations / honesty notes
* The independent GR Sturm-Liouville solver and the scheme were written by the same author (me); their agreement, the agreement with the report's own Galerkin solver (Newtonian limit) and, independent of the author, the agreement with the exact dust solutions of test (i) are the checks. The pressure terms are checked by (a), (b), (f), (g) but have no exact nonlinear solution to compare with.
* The lapse equation (`a'/a = -p'/(e+p)`) and the energy equation are machine-checked as identities I5, I6 of `verify_ms_equations.py` (divergence of the perfect-fluid stress tensor with the Christoffel symbols of the metric).  The viscosity is a standard von Neumann-Richtmyer form; results beyond the edge crossing depend on it (the bulk does not).
* Time normalisation: `t` is the proper time of the outermost shell (lapse 1 at the surface); to a distant observer multiply by 1/sqrt(1-2M/R) = 1.0013.
