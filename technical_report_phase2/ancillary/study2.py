import os
import euler_module as EM
import pickle, time

results = {}
for N in [3200, 6400]:
    t0=time.time()
    t, rho, v, p, log, status, wall, rho_c0 = EM.run(N, 0.40, wall_budget=float(os.environ.get('WALL_BUDGET', 260)))
    print(f"N={N:5d}  status={status}  t={t*1000:.2f}ms  rho_c/rho0={log[-1][1]:.6f}  mass_ratio={log[-1][2]:.6f}  wall={wall:.1f}s", flush=True)
    results[N] = dict(t=t, log=log, status=status, wall=wall, rho=rho, v=v, p=p)
    pickle.dump(results, open('study2_results.pkl','wb'))
    if time.time()-t0 > float(os.environ.get('WALL_BUDGET', 260)): break
