"""Equivalence test of the two versions of the seeded Eulerian module (added for the verification; no counterpart in the original archive).

euler_module_seeded.py           original module: stops only on NaN / negative density
euler_module_seeded_guarded.py   same code plus one more stop condition, max|v| > 1e12 m/s (status BLOWUP)

The guarded run is made first; then the original module is run for the same (N, epsilon) with a wall-clock budget of 3 x (wall time of the guarded run) + 10 s.
The test is that the logged states (t, rho_c/rho_c0, v_c) of the guarded run are IDENTICAL (exact floating-point equality) to the first entries of the
original run, i.e. that the added condition changes nothing before the stop.  The part of the original run after the stop is not compared (it is the
meaningless state reached after the run-away, see README.txt).
Usage (from this directory, after get_eigenfunction.py, PYTHONPATH=..:../linear_perturbation):  python3 check_guard_equivalence.py N epsilon t_final
"""
import sys, importlib
N = int(sys.argv[1]); eps = float(sys.argv[2]); tf = float(sys.argv[3])
G = importlib.import_module('euler_module_seeded_guarded')
t, rho, v, p, lg, sg, wg = G.run(N, eps, t_final=tf, wall_budget=1e9)
print(f"guarded:  N={N}, eps={eps:g}: status={sg}, t_reached={t:.4f} s, {len(lg)} steps, wall={wg:.1f} s")
O = importlib.import_module('euler_module_seeded')
budget = 3*wg + 10
t2, rho2, v2, p2, lo, so, wo = O.run(N, eps, t_final=tf, wall_budget=budget)
print(f"original: N={N}, eps={eps:g}: status={so}, t_reached={t2:.4f} s, {len(lo)} steps, wall={wo:.1f} s (budget {budget:.0f} s)")
n = len(lg)
same = len(lo) >= n and lo[:n] == lg
print(f"logged states of the guarded run identical to the first {n} of the original run: {same}")
print(f"last state of the guarded run: t={lg[-1][0]:.6f} s, rho_c/rho_c0={lg[-1][1]:.6g}, v_c={lg[-1][2]:.4g} m/s")
if len(lo) > n:
    print(f"state of the original run after the stop (not compared): step {n+1}: t={lo[n][0]:.6f} s, rho_c/rho_c0={lo[n][1]:.6g}; last step: rho_c/rho_c0={lo[-1][1]:.6g}, v_c={lo[-1][2]:.4g} m/s")
sys.exit(0 if same else 1)
