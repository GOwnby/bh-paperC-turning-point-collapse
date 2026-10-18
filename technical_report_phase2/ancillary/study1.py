import euler_module as EM

for N in [800, 1600]:
    t, rho, v, p, log, status, wall, rho_c0 = EM.run(N, 0.40, wall_budget=90)
    print(f"N={N:5d}  status={status}  t={t*1000:.2f}ms  rho_c/rho0={rho[EM.make_run(N)['g']] if False else log[-1][1]:.5f}  mass_ratio={log[-1][2]:.6f}  wall={wall:.1f}s")
