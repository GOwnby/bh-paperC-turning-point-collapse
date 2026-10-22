import numpy as np, json, sys
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from homologous_calibrated import CalibratedHomologousModel
G, c = bg.G, bg.c

def rho_max(mue):
    bg.mu_e = mue
    f = lambda lr: -build_metric(10**lr)['M']
    res = minimize_scalar(f, bounds=(13.2, 14.0), method='bounded', options={'xatol': 1e-4})
    return 10**res.x, -res.fun/bg.M_sun

def run(mue, ratio, rmax=None, d0=1e-6):
    bg.mu_e = mue
    if rmax is None: rmax, Mmax = rho_max(mue)
    else: Mmax=None
    rc = ratio*rmax
    bgd = build_metric(rc)
    out = dict(mu_e=mue, ratio=ratio, rho_max=rmax, M_max=Mmax, rho_c=rc, M=bgd['M']/bg.M_sun, R=bgd['R'], C=2*G*bgd['M']/(bgd['R']*c**2))
    for n in (1, 5):
        A, B, _, _ = build_RR_matrices(bgd, n); out[f's2_n{n}'] = eigh(A, B)[0][0]
    m = CalibratedHomologousModel(rc)
    k_target = -out['s2_n1']
    fcal, kN, kF = m.calibrate(target_k=k_target)
    out.update(f_cal=fcal, k_N=kN, k_F=kF, e_fold=1/np.sqrt(k_target))
    lam_h = lambda lam: 1-2*G*m.M/(lam*m.R*c**2)
    ev = lambda t, y: y[0]-0.15; ev.terminal = True
    smax = 40/np.sqrt(k_target)
    s = solve_ivp(lambda t, y: [y[1], m.lambda_ddot(y[0])], (0, smax), [1-d0, 0], method='Radau', events=ev, max_step=0.02, rtol=1e-9, atol=1e-13)
    lam = s.y[0]; th = s.t[np.argmax(lam < 0.5)]; out.update(t_half=th, t_015=s.t[-1], dt_plunge=s.t[-1]-th, Tdyn=np.sqrt(bgd['R']**3/(G*bgd['M'])))
    return out

if __name__ == "__main__":
    res = []
    rm2, M2 = rho_max(2.0); print("mu_e=2: rho_max", rm2, "Mmax", M2, flush=True)
    for ratio in ():
        r = run(2.0, ratio, rm2); res.append(r); print(json.dumps(r), flush=True)
    rmF, MF = rho_max(56/26); print("mu_e=56/26: rho_max", rmF, "Mmax", MF, flush=True)
    for ratio in (1.43,):
        r = run(56/26, ratio, rmF); r['M_max']=MF; res.append(r); print(json.dumps(r), flush=True)
    json.dump(res, open('survey.json', 'w'))
