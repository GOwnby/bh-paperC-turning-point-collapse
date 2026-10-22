# Code for Paper C

Scripts, logs and data of "Delayed Collapse of a Degenerate Core at Its Maximum Mass: General-Relativistic Instability and the Fold-Passage Law"
(doi:10.5281/zenodo.23455581). The table "result in the paper -> script -> output file", with the exact commands and run times, is in the
`README.md` of the repository root. Run every script from this directory (`cd code`). Python 3 with NumPy >= 2.0, SciPy, Matplotlib (`../requirements.txt`).

## Layout

- Libraries: `background.py` (Chandrasekhar equation of state, optional Salpeter lattice term `Zlat`, switch `use_eps` for the thermodynamically
  consistent energy density), `tov_solve.py` (TOV), `metric.py` (metric functions), `gr_rayleigh_ritz.py` (Chandrasekhar's pulsation equation, Rayleigh-Ritz),
  `homologous_calibrated.py` and `twomode.py` (calibrated nonlinear models), `exact_energy_model.py`, `lag_gr.py` (Misner-Sharp Lagrangian scheme),
  `survey.py` (composition survey driver used by `survey4.py`).
- Drivers: one script per result of the paper; see the root `README.md` for the mapping. The docstring at the top of each newer script states what it computes, its
  usage and its run time.
- Results: `*.log` are the screen output of the commands in the root README, `*.json` hold numerical results (for example `fold_step1.json` with the fold
  coefficients, read by the fold scripts), `*.npy` are the Lagrangian time series (`fold_inhom_N<N>_y0.1.npy`, `fold_BH_N<N>.npy`).
- `audit_pulsation/`: independent re-implementation of the pulsation eigenvalue (own equation of state `eos.py`, shooting on two formulations `pulse.py`/`pulseB.py`,
  finite elements `fem.py`, Newtonian check `newton.py`). Run its scripts from inside that directory (`cd audit_pulsation`); `proj*.py` import the modules of `code/`.

- `symbolic/`: sympy/mpmath scripts behind every "verified symbolically" statement of the paper (`sympy_eos_checks.py`, `sympy_weak_form.py`, `sympy_fold_scaling.py`), each with its log.
- Verification helpers: `eos_numbers.py`, `fold_numbers.py`, `pi_checks.py`, `rr_at.py`, `survey_spot.py`, `lag_checks.py` recompute individual quoted numbers; they are run by `../verify.sh`
  and compared with the manuscript values in `../verification/REPORT.md`.

## Notes on the shipped outputs

- `survey4.json`/`survey4.log` are the data of the composition survey table of the paper. The `e_fold` field of these files is the n = 1 (single trial function) value; the paper's 1/sigma column uses the five-function value.
- `fold_dyn2.py` and `fold_dyn2.json` are the ORIGINAL files of the archive (Zenodo v1.0.1), restored unchanged. The 14 rows of the `out` field (source of Table 3) are reproduced
  bit for bit by `fold_dyn2.py`. **Known issue:** the `PI` field of `fold_dyn2.json` (values 3.3157, 3.3454, 3.4659, 3.4534 for the start points s_0 = -20, -30, -45, -60) comes from an early pole
  estimator with an initial-slope sign error (the PI block of `fold_dyn2.py` starts the solution on the slow manifold with the wrong sign of the first-order slope) and is NOT the value used in the paper
  (the script as shipped gives 3.3206, 3.3503, 3.4708, 3.4583, see `fold_dyn2.log`; the shipped json was written by a slightly earlier version of the block). The paper's s_* = 3.41167 is produced
  by `pi_tritronquee.py`, which starts on the slow manifold with its first adiabatic correction. `verify.sh` regenerates the json into `verification/logs/` and leaves the shipped file as it is.
- `survey3.json` and `survey3.log` are superseded: they differ from `survey4.*` (for the ideal gas M_max = 1.41625 instead of 1.41661 and sigma^2 = 0 at 2.3132e13 instead of 2.3107e13 kg/m^3) and are quoted nowhere in the paper
  (see `survey3_SUPERSEDED.txt`); the table of the paper is `survey4.*`.
- Three figure generators are re-creations written for the verification, not the original scripts (the originals were not retained): `fig1_pekeris.py`, `fig8_collapse_timescale.py`,
  `figs_friedman.py`. The original figure files are in `../figures_original/`; see `../PROVENANCE.md`. `figs.py` and `figs_fold.py` are the original generators (paths changed only);
  the figures they write are pixel-identical to the originals.
- `sol_t.npy`, `sol_y.npy` (written by `integrate_collapse.py` and `final_single.py`) are scratch trajectories and are not tracked.
- Hand-off interface used by Paper D: `fold_inhom.handoff`, `fold_dyn.py`, `homologous_calibrated.py`, `lag_gr.py`, `metric.py`, `tov_solve.py`, `background*.py`
  and `fold_step1.json` must stay in this directory with these names.
