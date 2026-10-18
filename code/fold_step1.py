import numpy as np, json, background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
bg.mu_e=2.0
s2=lambda rc,n=5: eigh(*build_RR_matrices(build_metric(rc),n)[:2])[0][0]
f=lambda lr:-build_metric(10**lr)['M']
res=minimize_scalar(f,bounds=(13.30,13.45),method='bounded',options={'xatol':1e-8})
rm=10**res.x; Mm=-res.fun
print('rho_max',rm,'Mmax',Mm/bg.M_sun,flush=True)
out=dict(rho_max=rm,Mmax_kg=Mm,Mmax=Mm/bg.M_sun,pts=[])
for y in (0.002,0.005,0.01,0.02):
    rec={'y':y}
    for sgn in (+1,-1):
        d=build_metric(rm*(1+sgn*y)); rec[f's2_{sgn}']=s2(rm*(1+sgn*y)); rec[f'dM_{sgn}']=(Mm-d['M'])/bg.M_sun
    rec['kappa']=-(rec['s2_1']-rec['s2_-1'])/(2*y)   # sigma^2 = -kappa*y + ...
    rec['kappa2']=(rec['s2_1']+rec['s2_-1'])/(2*y*y)
    rec['a']=0.5*(rec['dM_1']+rec['dM_-1'])/y**2
    print(rec,flush=True); out['pts'].append(rec)
json.dump(out,open('fold_step1.json','w'))
