import sys,time
sys.dont_write_bytecode=True
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import numpy as np, background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import json
bg.mu_e=2.0; bg.Zlat=0.0
t=time.time()
res=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(13.30,13.45),method='bounded',options={'xatol':1e-9})
rm=10**res.x; Mm=-res.fun
print('rho_max',rm,'Mmax',Mm/bg.M_sun,time.time()-t,flush=True)
def s2(rc,n=5): return eigh(*build_RR_matrices(build_metric(rc),n)[:2])[0][0]
for n in (5,):
  for y in (0.005,):
    e={}
    for sg in (1,-1):
        d=build_metric(rm*(1+sg*y)); e[sg]=(s2(rm*(1+sg*y),n),(Mm-d['M'])/bg.M_sun)
    print(n,y,'s2+-',e[1][0],e[-1][0],'kappa',-(e[1][0]-e[-1][0])/(2*y),'a',0.5*(e[1][1]+e[-1][1])/y**2,flush=True)
