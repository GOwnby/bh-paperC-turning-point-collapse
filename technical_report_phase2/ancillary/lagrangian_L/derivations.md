# Misner-Sharp system: statement, derivation of the momentum equation, and discrepancies with Section 2 of the report

Units G = c = 1, time unit L/c, length unit L. Metric (comoving, label r):
`ds^2 = -a^2 dt^2 + X^2 dr^2 + R^2 dOmega^2`, fluid at rest in (t, r), `T^mu_nu = diag(-e, p, p, p)`.
Definitions: `U = (1/a) dR/dt` (derivative with respect to the proper time of a fluid element), `Gamma = R'/X`,
Misner-Sharp mass `m` by `Gamma^2 = 1 + U^2 - 2m/R`.

## The system (Misner & Sharp 1964, Phys. Rev. 136, B571; May & White 1966, Phys. Rev. 141, 1232)

| # | equation | origin |
|---|----------|--------|
| 1 | `dR/dt = a U` | definition of U |
| 2 | `m' = 4 pi e R^2 R'` | `G^t_t = -8 pi e` (identity I1) |
| 3 | `dm/dt = -4 pi p R^2 dR/dt = -4 pi a p R^2 U` | `G^r_r = 8 pi p` (identity I2) |
| 4 | `dGamma/dt = a' (dR/dt) / (a X) = a' U / X` | `G_tr = 0` (identity I3) |
| 5 | `a'/a = -p'/(e+p)` | radial component of `nabla_mu T^{mu r} = 0` |
| 6 | `dU/dt = -a [ (m + 4 pi R^3 p)/R^2 + Gamma^2 p' /((e+p) R') ]` | consequence of 1-5 (below) |
| 7 | `p = P(e)` | closure (the report's EOS) |

Derivation of 6 (pure algebra, checked numerically as identity I4): differentiate `Gamma^2 = 1 + U^2 - 2m/R`,
`Gamma Gamma_dot = U U_dot - m_dot/R + m R_dot/R^2`; insert 3 (`m_dot = -4 pi p R^2 a U`), 4 (`Gamma_dot = a' U/X`) and `R_dot = aU`:
`U_dot = Gamma a'/X - a (m/R^2 + 4 pi p R)`; with 5 and `Gamma/X = Gamma^2/R'`:
`U_dot = -a [ (m + 4 pi R^3 p)/R^2 + Gamma^2 p'/((e+p) R') ]`.
The factor in front is **a** (one power).

Barotropic fluid: with `rho_N` the conserved number-density function (`d e = (e+p) d rho_N / rho_N`) and `h = (e+p)/rho_N`:
`dp/(e+p) = d ln h`, hence Eq. 5 integrates to `a h = const` (hydrostatic equilibrium: redshifted enthalpy constant), and the
rest-mass content of a shell is `dM_b = 4 pi R^2 rho_N dR/Gamma` (this is what the scheme keeps fixed).
Static limit: `dH/dr = -(m + 4 pi r^3 p)/(r^2 (1-2m/r))` with `H = ln h` (TOV in enthalpy form).

Machine verification of 2, 3, 4 and of the algebra behind 6: `src/verify_ms_equations.py` (exact second-order Taylor/jet arithmetic for
the Einstein tensor of the above metric at random points, no computer algebra system); `logs/verify_ms_equations.log`.
The Euler/lapse equation 5 is the r-component of `nabla_mu T^mu_nu = 0`: for `T^mu_nu = diag(-e,p,p,p)` one has `nabla_mu T^mu_r = p' + sum_mu Gamma^mu_{mu r}(p - T^mu_mu) = p' + Gamma^t_{tr}(e+p)` with `Gamma^t_{tr} = a'/a`, and the t-component gives `edot = -(e+p)(Xdot/X + 2 Rdot/R)`, the continuity equation behind the conserved `dM_b`.  Both are checked to rounding by exact jets in `src/verify_ms_equations.py` (identities I5, I6, `logs/verify_ms_equations.log`).

## Discrepancies with the report (Section 2) and with ancillary/evolve2.py

1. **Eq. (constraint) of the report, "`R'/X` is independent of tau along each shell", is false in general.** Identity I3 gives
   `G_tr = -(2X/R) (dGamma/dt - a' (dR/dt)/(a X))`, so `G_tr = 0` (the momentum constraint) implies `dGamma/dt = a' (dR/dt)/(aX)`, which vanishes only if
   `a' = 0` (dust or uniform pressure) or the shell is momentarily at rest. The argument given in the report ("follows from `T^tau_r = 0` ... does not reference the
   lapse") overlooks that `G_tr` contains the lapse gradient. `collapse2.py`/`evolve2.py` implement the false statement: `X = R'/Gamma_shell` with `Gamma_shell`
   frozen at its initial value (`X` in the pressure-gradient term), instead of `Gamma^2 = 1 + U^2 - 2m/R` at the current time.
   Size (`src/test_j_gamma_size.py`, `results/test_j_gamma_size.json`): up to t = 50 L/c on the report's data (N = 400) `max |Gamma^2(t) - Gamma^2(0)|` = 2.7e-4 (pure) and 2.5e-4 (neutronised), with `|U| <= 0.004 c`; the largest `2m/R` on the same interfaces is 4.5e-3 and 3.9e-3, so the variation is about 6 % of the relativistic correction (larger later in the collapse).
2. **Factor `a` versus `a^2` in the Euler equation.** The report writes `dU/dtau = -a [...]` and `dm/dtau = -4 pi p R^2 U` (both consistent only if `U` is the
   *coordinate-time* derivative of R), but `evolve2.py` uses `dR/dt = a U`, `dm/dt = -4 pi a p R^2 U` (so `U` = proper-time derivative) together with
   `dU/dt = -a^2 [...]`. With `U` the proper-time derivative the correct factor is `a` (derivation above). The code therefore multiplies the acceleration by an extra
   `a = 1/h <= 1` (0.3-0.35 % at the centre for rho_c = 5e12): internally inconsistent, small. (Effect measured in `logs/repro_spline.log`, variant "a^1 instead of a^2".)
3. **Time variable.** With the lapse normalised to 1 at the free surface (`a = 1/h`, `h -> 1` in vacuum) the time `t` is the proper time of the outermost shell,
   not the time of a distant observer; the two differ by `sqrt(1-2M/R) = 0.9987`. `L/c` is the unit throughout (`L` = the report's Newtonian `R_star`, 1601.9 km pure,
   1615.8 km neutronised). Mode frequencies from the Schwarzschild-time linear equations are converted with `omega^2_here = omega^2_inf/(1-2M/R)`.
4. **Initial data.** The report's Newtonian-hydrostatic profile is not a solution of the Misner-Sharp hydrostatic equations: the Newtonian force balance
   `p' = -rho m/r^2` differs from `p' = -(e+p)(m + 4 pi r^3 p)/(r^2 (1-2m/r))` by a relative `~0.4-0.5 %` (all three corrections have the same sign: more gravity), see
   `logs/test_f_forced.log` (initial acceleration -0.0046 g pure, -0.0040 g neutronised, independent of N, versus < 1e-4 g for pointwise TOV data) and `logs/test_a_static.log`. For a configuration that is marginally stable this is not negligible (see the results).
