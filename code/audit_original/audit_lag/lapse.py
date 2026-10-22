import sys, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
from lag_jac import jac_eigs
G,c=bg.G,bg.c
class LagUnif(Lag):  # variant: spatially uniform lapse (violates alpha*h=const) -- own copy, scratch only
    def f(self,t,y,subtract=True,bracket=False):
        # call base f but rescale alpha via monkeypatching thermo-h: emulate by multiplying result by h-ratio
        out=super().f(t,y,subtract,bracket)
        if bracket: return out
        N=self.N; r,U,m=self.unpack(y); Gam,rho,P,e=self.zone(r,U,m); h=(e+P)/rho
        hz=h/c**2
        hj=np.empty(N+1); hj[1:N]=2/(1/hz[1:]+1/hz[:-1]); hj[N]=1  # inverse of the arithmetic mean of 1/h
        am=np.empty(N+1); am[1:N]=0.5*(1/hz[1:]+1/hz[:-1]); am[N]=1.0
        fac=1/am[1:]  # undo the c^2/h factor at interfaces
        o=out.copy(); o[:N]*=fac; o[N:2*N]*=fac; o[2*N:]*=fac
        return o
N=int(sys.argv[1])
for name,cls,fa in (('proper',Lag,None),('alpha_s=1 (h-profile kept)',Lag,1.0),('uniform lapse alpha_s',LagUnif,None)):
    L=cls(5e13,N=N,pert=None)
    if fa is not None: L.freeze_alpha=fa
    # recompute S for frozen alpha is irrelevant (bracket only)
    ev,V=jac_eigs(L); print(name,'N=%d growth %.5f'%(N,ev.real.max()),flush=True)
