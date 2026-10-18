import numpy as np, json
import background as bg
from survey import run
from scipy.optimize import minimize_scalar, brentq
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
cases = {'ideal-C(Z=0,mue=2)':(0,2.0),'He4':(2,2.0),'C12':(6,2.0),'O16':(8,2.0),'Fe56':(26,56/26)}
res=[]
for name,(Z,mue) in cases.items():
    bg.Zlat=float(Z); bg.mu_e=mue
    r=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(12.8,14.2),method='bounded',options={'xatol':1e-4})
    rmax=10**r.x
    s2=lambda rc: eigh(*build_RR_matrices(build_metric(rc),5)[:2])[0][0]
    rzero=brentq(s2,0.5*rmax,1.5*rmax,xtol=1e9)
    print(name,"M peak rho",rmax,"M",-r.fun/bg.M_sun,"sigma2=0 rho",rzero,flush=True)
    for ratio in (1.1,1.5,2.0,3.0,5e13/rmax):
        bg.Zlat=float(Z); bg.mu_e=mue
        o=run(mue,ratio,rmax); o.update(name=name,Z=Z,M_max=-r.fun/bg.M_sun,rho_zero=rzero); res.append(o)
        print(json.dumps({k:(round(v,5) if isinstance(v,float) else v) for k,v in o.items()}),flush=True)
json.dump(res,open('survey4.json','w'))
