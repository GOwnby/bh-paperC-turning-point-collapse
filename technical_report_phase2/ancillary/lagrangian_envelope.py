"""RE-CREATION (not the original code) of the Lagrangian shell-crossing runs of the technical report, Sec. 4 ('Lagrangian evolution: shell crossing in the
envelope, not the core').  The original driver and the neutronized-EOS variants evolve_neu.py / collapse_neu.py were not retained (see README.txt).

What is re-created, and how:
  * initial data: exactly collapse2.py (rho_c = 5e12 kg/m^3, 300 comoving shells from 1 % of R to R, Newtonian hydrostatic profile), with the
    pure Chandrasekhar EOS ('pure', initial_profile.py) or the neutronized EOS ('neu', initial_profile_neu.py / eos_neutronized.py);
  * right-hand side: exactly evolve2.py (Misner-Sharp, Eq. of Sec. 2 of the report) with the pressure P(rho) of the chosen EOS;
  * time stepping: classical fourth-order Runge-Kutta with a fixed step dt (the report mentions RK4 sub-stages); time unit L/c with L = R_star;
  * stop at the first time min_i R'(r_i) < 0 (shell crossing; R' by the same cubic-spline derivative as the right-hand side).
It prints the crossing time in units of L/c and in ms, the shell where R' first turns negative, and the change of the central density.
usage: python3 lagrangian_envelope.py pure|neu [dt]       (the values recorded originally, 14.2 (pure) and 9.1 (neutronized) in these units, are no longer quoted as results in the report)
Run from this directory.

RESULT OF THIS RE-CREATION: it does NOT reproduce the numbers of the report.  With dt = 1e-3 (the largest stable step; dt = 4e-3 is unstable and 'crosses' at
t ~ 1.9 at the surface) it finds, for the pure EOS, R' < 0 at t = 14.78 (79 ms), at the innermost shell (not an outer shell), after the central density has
oscillated by +-20 % (report: 14.2, an outer shell, central density within 1-2 %), and for the neutronized EOS R' < 0 at t = 2.19 (11.8 ms) at the outermost
shell (report: 9.1, i.e. 49 ms).  The original driver, its step size and the stopping rule were not retained, so the report (Sec. 4) states no crossing time, shell or neutronization effect from the Lagrangian runs; this file is kept only as the record of the attempt, and nothing else uses it.
"""
import sys, time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
import eos as E

kind = sys.argv[1] if len(sys.argv) > 1 else 'pure'
dt = float(sys.argv[2]) if len(sys.argv) > 2 else 2e-3
if kind == 'neu':
    import initial_profile_neu as IP, eos_neutronized as EOS
    P_of_rho = EOS.P_of_rho
else:
    import initial_profile as IP
    P_of_rho = IP.P_of_rho

rho_c = 5e12
sol = IP.build_star(rho_c, rmax=2e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 0)
r_surf = r_raw[-1]
N = 300
r_grid = np.linspace(0.01*r_surf, r_surf, N)
m_grid = np.interp(r_grid, r_raw, m_raw)
L = r_surf
r_lab = r_grid/L
m0 = (m_grid*E.G/E.c**2)/L
Gamma_shell = np.sqrt(np.maximum(1 - 2*m0/r_lab, 1e-10))
Gscale = E.G*L**2/E.c**2
print(f'[{kind}] Star: M={m_raw[-1]/E.Msun:.4f} Msun, R={r_surf/1e3:.1f} km; time unit L/c = {L/E.c*1e3:.3f} ms; dt = {dt} (L/c)')

def p_dimless(rho_d):
    return Gscale*(P_of_rho(np.maximum(rho_d, 0)/Gscale)/E.c**2)
def derivs(y, r): return CubicSpline(r, y)(r, 1)
def weight_of(rho, p_d):
    denom = rho + p_d; threshold = 1e-6*np.max(rho)
    return denom, np.maximum(denom, 1e-300), denom**2/(denom**2 + threshold**2)
def rhs(t, y):
    R, m, U = y[0:N], y[N:2*N], y[2*N:3*N]
    Rp = derivs(R, r_lab); mp = derivs(m, r_lab)
    rho = np.maximum(mp/(4*np.pi*np.maximum(R, 1e-12)**2*np.where(np.abs(Rp) > 1e-12, Rp, 1e-12)), 0)
    p_d = p_dimless(rho)
    dp_dr = derivs(p_d, r_lab)
    denom, denom_safe, w = weight_of(rho, p_d)
    integrand = w*(-dp_dr/denom_safe)
    a = np.clip(np.exp(cumulative_trapezoid(integrand[::-1], r_lab[::-1], initial=0.0)[::-1]), 1e-6, 1e2)
    X = Rp/np.maximum(Gamma_shell, 1e-12)
    pgrad = w*(Gamma_shell/denom_safe)*(dp_dr/np.where(np.abs(X) > 1e-12, X, 1e-12))
    return np.concatenate([a*U, -4*np.pi*a*p_d*R**2*U, -a**2*((m + 4*np.pi*R**3*p_d)/np.maximum(R, 1e-12)**2 + pgrad)]), rho, Rp

y = np.concatenate([r_lab.copy(), m0.copy(), np.zeros(N)])
_, rho0, Rp0 = rhs(0.0, y)
t = 0.0; t0 = time.time(); nstep = 0; status = 'tmax'
while t < 30.0:
    k1, rho, Rp = rhs(t, y)
    if np.min(Rp) < 0 or not np.all(np.isfinite(y)):
        status = 'shell crossing (R\' < 0)'; break
    k2 = rhs(t + dt/2, y + dt/2*k1)[0]; k3 = rhs(t + dt/2, y + dt/2*k2)[0]; k4 = rhs(t + dt, y + dt*k3)[0]
    y = y + dt/6*(k1 + 2*k2 + 2*k3 + k4); t += dt; nstep += 1
    if nstep % 500 == 0:
        print(f'  t={t:8.3f}  rho_c/rho_c0={rho[0]/rho0[0]:.5f}  min R\'={np.min(Rp):.4f}', flush=True)
i = int(np.argmin(Rp))
print(f'[{kind}] {status} at t = {t:.3f} (L/c) = {t*L/E.c*1e3:.2f} ms; shell {i} of {N} (r/R = {r_lab[i]:.3f}, rho/rho_c = {rho[i]/rho0[0]:.3e});'
      f' central density changed by {100*(rho[0]/rho0[0]-1):+.2f} %; wall {time.time()-t0:.0f} s')
