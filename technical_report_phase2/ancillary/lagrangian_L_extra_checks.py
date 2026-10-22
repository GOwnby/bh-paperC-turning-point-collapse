"""Three statements of Sec. 4 of the technical report that are not in any log of ancillary/lagrangian_L (written for the verification, round 2 of the audit).

   (1) intrinsic first-step kick of an unfloored surface zone, 'about 0.57 g*dt for every N' (lagrangian_L/src/test_d_edge.py runs N = 400 only);
   (2) 1 - a = 1 - 1/h at the centre, '0.3-0.35 %' (lagrangian_L/derivations.md; a*h = const with h = 1 at the surface);
   (3) the fixed-step 'core crossing' of the neutronised star (test_h_collapse.py, N = 400, factors 1, 1/2, 1/4): zones with R' <= 0 at the crossing step.

It uses src/ of lagrangian_L unmodified and the same set-up as the tests above.  Usage (from this directory):  python3 lagrangian_L_extra_checks.py lagrangian_L
Run time about 1 min.  Output: lagrangian_L_extra_checks.log (checked by verify.py: LG-kick-allN, LG-fixed-zone)."""
import sys, os, numpy as np
sys.dont_write_bytecode = True
root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else 'lagrangian_L'); sys.path.insert(0, os.path.join(root, 'src')); os.chdir(os.path.join(root, 'src'))
from harness import *
for kind in ('pure', 'neu'):
    S = Setup(kind)
    print(f'{kind}: central enthalpy H_c = ln h_c = {S.tov["H_c"]:.5f}  ->  1 - 1/h_c = {100*(1-np.exp(-S.tov["H_c"])):.3f} %')
    for N in (100, 200, 400, 800, 1600):
        d = S.newton_data(N)
        star2 = S.make(dict(d)); star2.start()
        dt = star2.dt_cfl(star2.geo, 0.4)
        g = star2.geo['m'][-1] / star2.R[-1]**2
        star2.step(dt)
        print(f'   N = {N:5d}: first-step kick = {np.max(np.abs(star2.U[-5:])) / (g * dt):.4f} g dt')
print()
S = Setup('neu')
for ds in (1.0, 0.5, 0.25):
    d = S.newton_data(400)
    star = S.make(d, cq2=2.0)
    r = evolve(star, 120.0, nout=240, dt_scale=ds)
    gaps = np.diff(star.R)
    print(f'neutronised N = 400, fixed dt x {ds}: status {r["status"]}, t_cross {r["t_cross"]:.3f}, zone of minimum gap {int(np.argmin(gaps))}, zones with R\' <= 0 at that step: {np.where(gaps <= 0)[0].tolist()}')
