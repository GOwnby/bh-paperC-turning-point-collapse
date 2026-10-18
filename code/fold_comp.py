import numpy as np, json, background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
cases={'ideal(Z=0)':(0,2.0,None),'He4':(2,2.0,1.39e11),'C12':(6,2.0,3.97e10),'O16':(8,2.0,1.94e10),'Fe56':(26,56/26,1.18e9)}
s2=lambda rc,n=5: eigh(*build_RR_matrices(build_metric(rc),n)[:2])[0][0]
out={}
for name,(Z,mue,thr) in cases.items():
    bg.Zlat=float(Z); bg.mu_e=mue
    res=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(13.2,13.6),method='bounded',options={'xatol':1e-8})
    rm=10**res.x; Mm=-res.fun; rec=dict(Z=Z,mu_e=mue,rho_max=rm,Mmax=Mm/bg.M_sun,thr_gcc=thr,pts=[])
    for y in (0.005,0.01):
        e={}
        for sg in (1,-1):
            d=build_metric(rm*(1+sg*y)); e[sg]=(s2(rm*(1+sg*y)),(Mm-d['M'])/bg.M_sun)
        rec['pts'].append(dict(y=y,kappa=-(e[1][0]-e[-1][0])/(2*y),a=0.5*(e[1][1]+e[-1][1])/y**2,s2p=e[1][0],s2m=e[-1][0]))
    rec['kappa']=np.mean([p['kappa'] for p in rec['pts']]); rec['a']=np.mean([p['a'] for p in rec['pts']])
    out[name]=rec; print(name,{k:v for k,v in rec.items() if k!='pts'},flush=True)
json.dump(out,open('fold_comp.json','w'))
