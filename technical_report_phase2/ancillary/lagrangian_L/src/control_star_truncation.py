"""Diagnosis: the report's rho_c = 1e11 'stable control' (linear_perturbation/stable_crosscheck.py, euler_module.py/control_test.py) builds the
star with rmax = 3e6 m.  The true star is larger, so the integration stops at rmax with non-zero density: a TRUNCATED star with an artificial
hard edge (the same failure the report describes and claims fixed for its polytrope builder, Sec. 'two bugs')."""
import numpy as np, os, io, contextlib, sys
import report_galerkin as RG, eosmodel as M
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'ref'))
with contextlib.redirect_stdout(io.StringIO()):
    import initial_profile_neu as IPN, initial_profile as IP
for kind, IPm in (('neu', IPN), ('pure', IP)):
    print(f'--- {kind}, rho_c = 1e11 kg/m^3')
    for rmax, label in ((3e6, "report's rmax = 3e6 m"), (2.5e7, 'rmax = 2.5e7 m (surface reached)')):
        sol = IPm.build_star(1e11, rmax=rmax, n=int(rmax / 250) if rmax > 5e6 else 8000)
        rho_end = sol.y[1, -1]
        print(f'  {label}: stopped at R = {sol.t[-1]/1e3:.1f} km, M = {sol.y[0,-1]/M.MSUN:.4f} Msun, rho at the last point = {rho_end:.3e} kg/m^3 ({rho_end/1e11:.2e} of rho_c), surface event reached: {sol.status==1}')
        w, R, Mt = RG.galerkin_modes(kind, 1e11, 1000, rmax=(None if rmax == 3e6 else 2.5e7))
        print(f'      report Galerkin (N=1000): omega^2 = {np.array2string(w[:3], precision=4)} s^-2   (fundamental period {2*np.pi/np.sqrt(w[0]):.3f} s)')
