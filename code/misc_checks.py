"""Small numerical statements of Paper C that are not produced by the larger scripts (runtime ~2 min):
 1. Newtonian Chandrasekhar-mass limit of the electron-gas EOS used here (Sec. 2.1: 1.456 Msun, mu_e=2);
 2. Newtonian Rayleigh-Ritz sigma^2 at the test configuration (positive: Newtonian solver reports 'stable', Sec. 3.3);
 3. eigenfunction shape xi_GR(r)/r (center/surface = 1.133) and the position of lambda_crit (r = 0.42 R), Sec. 5.1-5.2;
 4. flatness of M(rho_c) near the maximum (relative mass deficit at +-25 %, +50 %), Sec. 3.6 (ideal gas and Salpeter 12C);
 5. equipartition rms velocity of the fold mode (Sec. 5.5, Limits (ii))."""
import numpy as np, json, os
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
import background as bg
from metric import build_metric
import gr_rayleigh_ritz as rr
from homologous_calibrated import CalibratedHomologousModel
HERE = os.path.dirname(os.path.abspath(__file__)); G, c = bg.G, bg.c
bg.mu_e = 2.0; bg.Zlat = 0.0
# 1. Newtonian limit
def newt_mass(rc):
    xc = bg.x_of_rho(rc)
    rhs = lambda r, y: [-G*bg.rho_of_x(y[0])*y[1]/r**2/(bg.A_P*bg.dfdx(y[0])), 4*np.pi*r**2*bg.rho_of_x(y[0])]
    ev = lambda r, y: y[0] - 1e-4*xc; ev.terminal = True
    s = solve_ivp(rhs, [1.0, 1e9], [xc, 4/3*np.pi*rc], events=ev, rtol=1e-10, atol=[1e-14, 1e-3], method='DOP853')
    return s.y_events[0][0][1]/bg.M_sun
for rc in (5e12, 1e14, 1e16, 1e18, 1e20): print(f"Newtonian hydrostatic mass, rho_c={rc:.0e}: {newt_mass(rc):.5f} Msun")
# 2. Newtonian RR (c -> infinity limit of the same variational form), as in figs.py
src = open(os.path.join(HERE, 'figs.py')).read()
exec(src[:src.index("import os\nout=")]); exec(src[src.index('def newtonian_s2'):src.index('rcs=[')])
print("Newtonian RR sigma^2 at 5e13 (compactness):", newtonian_s2(5e13))
# 3. eigenfunction shape and lambda_crit location
d = build_metric(5e13); A, B, xis, r = rr.build_RR_matrices(d, 5); w, v = eigh(A, B)
xi = sum(v[i, 0]*xis[i] for i in range(5)); q = xi[1:]/r[1:]
print("xi_GR/r center/surface =", q[0]/q[-1])
m = CalibratedHomologousModel(5e13); comp = 2*G*m.Mr_bg[1:]/(m.r[1:]*c**2); j = np.argmax(comp)
print("lambda_crit =", comp.max(), " at r/R =", m.r[1:][j]/m.R, " surface compactness =", comp[-1], " T_dyn =", np.sqrt(m.R**3/(G*m.M)))
# 4. flatness of the mass curve
fc = json.load(open(os.path.join(HERE, 'fold_comp.json')))
for name, Z in (('ideal(Z=0)', 0.0), ('C12', 6.0)):
    bg.Zlat = Z; rm = fc[name]['rho_max']; Mm = fc[name]['Mmax']*bg.M_sun
    print(name, {y: float((Mm - build_metric(rm*(1 + y))['M'])/Mm) for y in (-0.25, 0.25, 0.5)})
bg.Zlat = 0.0
# 5. equipartition rms velocity of the fold mode
kB = 1.380649e-23
mf = CalibratedHomologousModel(json.load(open(os.path.join(HERE, 'fold_step1.json')))['rho_max'])   # fold configuration
print("fold configuration: model moment of inertia I =", mf.I, " sqrt(kB T/I) at 1e8 K =", np.sqrt(kB*1e8/mf.I))
