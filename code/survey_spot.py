"""Spot check of two rows of the composition survey (Table 2), recomputed from scratch with the same code path as survey4.py:
the ideal gas at offset 2.0 and Fe56 at the fiducial density rho_c = 5e13 kg/m^3.  The JSON lines are compared by verify.py with survey4.json."""
import json
import background as bg
from survey import run
from scipy.optimize import minimize_scalar
from metric import build_metric
for name, Z, mue, ratio in (('ideal-C(Z=0,mue=2)', 0, 2.0, 2.0), ('Fe56', 26, 56/26, None)):
    bg.Zlat = float(Z); bg.mu_e = mue
    r = minimize_scalar(lambda lr: -build_metric(10**lr)['M'], bounds=(12.8, 14.2), method='bounded', options={'xatol': 1e-4})
    rmax = 10**r.x
    if ratio is None: ratio = 5e13/rmax
    bg.Zlat = float(Z); bg.mu_e = mue
    o = run(mue, ratio, rmax); o.update(name=name, Z=Z, M_max=-r.fun/bg.M_sun)
    print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in o.items()}), flush=True)
