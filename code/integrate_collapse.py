import numpy as np
from scipy.integrate import solve_ivp
from homologous_calibrated import CalibratedHomologousModel
import background as bg

G, c = bg.G, bg.c

m = CalibratedHomologousModel(5e13)
f_cal, k_N, k_F = m.calibrate_to_RR()
sigma = np.sqrt(m.target_k)
print(f"f_cal={f_cal:.4f}, sigma={sigma:.4f} s^-1, e-folding time 1/sigma={1/sigma:.4f} s")

# Horizon-proximity guard: stop if the lapse factor at the stellar surface
# (lam*R) approaches zero (collapse through the would-be horizon).
comp_max = np.max(2*G*m.Mr_bg[1:]/(m.r[1:]*c**2))   # max_r 2GM_r/(r c^2): lapse factor first vanishes at this interior shell
def horizon_factor(lam):
    return 1 - comp_max/lam

print(f"Background compactness at lambda=1: 1-2GM/(Rc^2) = {horizon_factor(1.0):.6f}")

def rhs(t, y):
    lam, lamdot = y
    if lam <= 1e-6:
        return [0.0, 0.0]
    return [lamdot, m.lambda_ddot(lam)]

lam_crit = comp_max
print(f"Critical lambda = max_r 2GM_r/(r c^2) (interior shell r~0.42R, where the frozen-metric lapse factor first vanishes) = {lam_crit:.6e}")

def event_horizon(t, y):
    return horizon_factor(y[0]) - 0.02   # stop with a safety margin before the true singularity
event_horizon.terminal = True
event_horizon.direction = -1

def event_small(t, y):
    return y[0] - 1e-3
event_small.terminal = True
event_small.direction = -1

# Start from a tiny inward density perturbation, zero velocity -- the
# generic initial condition for a star nudged off its (unstable) critical
# equilibrium; the growing mode dominates after a few e-folding times
# regardless of the exact split between growing/decaying components.
delta0 = 1e-6
y0 = [1.0-delta0, 0.0]

t_span = (0.0, 3.0)
sol = solve_ivp(rhs, t_span, y0, method='Radau', events=[event_horizon, event_small],
                 max_step=0.02, rtol=1e-9, atol=1e-13, dense_output=True)

print("\nIntegration status:", sol.status, sol.message)
print("t_events (horizon, small):", sol.t_events)
print(f"Final t = {sol.t[-1]:.6f} s,  lambda_final = {sol.y[0][-1]:.6e}")

# Report lambda(t) at a series of checkpoints
import numpy as np
checkpoints = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
print("\n t(s)    lambda        1-lambda      central_density/rho_c")
for tc in checkpoints:
    if tc <= sol.t[-1]:
        lam_c = sol.sol(tc)[0]
        rho_ratio = 1.0/lam_c**3 if lam_c > 0 else np.inf
        print(f"{tc:5.2f}  {lam_c:.6e}  {1-lam_c:.6e}  {rho_ratio:.4e}")

np.save('sol_t.npy', sol.t)
np.save('sol_y.npy', sol.y)
