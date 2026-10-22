import sys, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
from lag_jac import jac_eigs
for p in (1.0,3.0):
    v=[]
    for N in (100,200,400):
        L=Lag(5e13,N=N,p=p,pert=None); ev,_=jac_eigs(L); v.append(ev.real.max())
    d=np.diff(v); print('p=%.1f growth'%p,np.round(v,5),'ratio',d[0]/d[1],'-> limit %.5f'%(v[-1]+d[-1]/(d[0]/d[1]-1)),flush=True)
