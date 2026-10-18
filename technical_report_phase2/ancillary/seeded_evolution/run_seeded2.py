import importlib
import euler_module_seeded as EMS
importlib.reload(EMS)
import numpy as np, pickle

eps = 1e-4
results = {}
for N in [800, 1600]:
    t, rho, v, p, log, status, wall = EMS.run(N, eps, t_final=3.0, wall_budget=150)
    results[N] = log
    print(f"N={N}: status={status}, t_reached={t:.3f}s, wall={wall:.1f}s")
    for i in range(0, len(log), max(1,len(log)//12)):
        tt, rc, vc = log[i]
        print(f"    t={tt:.3f}s  rho_c/rho_c0={rc:.6f}  v_c={vc:+.3e} m/s")

pickle.dump(results, open('seeded_run2.pkl','wb'))
