# Code for Paper C (merged from former "C1" + C2)
Union of the former "C1" (linear instability: Newtonian/GR pulsation solvers, Rayleigh-Ritz, turning-point test) and
C2 (nonlinear timescale model, Painleve-I fold passage, dissipation, Lagrangian test) code. See the old per-folder READMEs in
../../archive_split_C1_C2_C3/ for script-to-figure mapping. Python 3, numpy/scipy.

## Known issue in saved outputs
The `PI` field in `fold_dyn2.json` (values 3.3157, 3.3454, 3.4659, 3.4534) comes from an early pole estimator with an initial-slope sign error
and is NOT the value used in the paper. The paper's s_* = 3.41167 is produced by `pi_tritronquee.py` (agrees with the Boutroux
tritronquee data, 6^{1/5} x* with x* = 2.384169). Use `pi_tritronquee.py`, not that field.
