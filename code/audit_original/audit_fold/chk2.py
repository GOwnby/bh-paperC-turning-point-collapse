import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import numpy as np, background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
which=sys.argv[1]
cfg={'C12':(6,2.0),'Fe56':(26,56/26),'ideal':(0,2.0)}[which]
bg.Zlat=float(cfg[0]); bg.mu_e=cfg[1]
res=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(13.2,13.6),method='bounded',options={'xatol':1e-9})
rm=10**res.x; Mm=-res.fun
print(which,'rho_max',rm/1e3,'g/cc','Mmax',Mm/bg.M_sun,flush=True)
ns=(3,5,8) if which=='ideal' else (5,)
for n in ns:
  y=0.005; e={}
  for sg in (1,-1):
    d=build_metric(rm*(1+sg*y)); e[sg]=(eigh(*build_RR_matrices(d,n)[:2])[0][0],(Mm-d['M'])/bg.M_sun)
  print(which,'n',n,'kappa',-(e[1][0]-e[-1][0])/(2*y),'a',0.5*(e[1][1]+e[-1][1])/y**2,flush=True)
