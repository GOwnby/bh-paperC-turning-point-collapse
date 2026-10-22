# sensitivity ("mutation") test: does the linear validation discriminate wrong GR terms?  own subclass; module untouched
import sys, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
from lag_jac import jac_eigs
G,c=bg.G,bg.c
class Mut(Lag):
    flag=None
    def f(self,t,y,subtract=True,bracket=False):
        N=self.N; r,U,m=self.unpack(y); Gam,rho,P,e=self.zone(r,U,m); enth=e+P; h=enth/rho
        Ms=m[-1]; alpha_s=(1-2*G*Ms/(r[-1]*c**2))/Gam[-1]
        a_zone=alpha_s*c**2/h
        alpha=np.empty(N+1); alpha[1:N]=0.5*(a_zone[1:]+a_zone[:-1]); alpha[N]=alpha_s; alpha[0]=a_zone[0]
        Pj=np.zeros(N+1); Pj[1:N]=0.5*(P[1:]+P[:-1]); Pj[N]=0.0; Pj[0]=P[0]
        dmu=self.dmu; Mnode=np.empty(N); Mnode[:N-1]=0.5*(dmu[:-1]+dmu[1:]); Mnode[N-1]=0.5*dmu[-1]
        hz=h*dmu; hnode=np.empty(N); hnode[:N-1]=0.5*(hz[:-1]+hz[1:])/Mnode[:N-1]; hnode[N-1]=h[-1]
        dP=np.empty(N); dP[:N-1]=P[1:]-P[:-1]; dP[N-1]=-P[-1]
        gp=0.0 if self.flag=='noPterm' else 1.0
        Gf=Gam[1:]**(1 if self.flag=='Gamma1' else 1)   # baseline
        if self.flag=='Gam0': Gf=1.0
        acc=-Gf*c**2*4*np.pi*r[1:]**2*dP/(hnode*Mnode)-G*(m[1:]+gp*4*np.pi*r[1:]**3*Pj[1:]/c**2)/r[1:]**2
        if subtract: acc=acc-self.S
        if bracket: return acc
        a1=alpha[1:] if self.flag!='noalpha' else 1.0
        dm=-a1*4*np.pi*r[1:]**2*Pj[1:]*U[1:]/c**2
        if self.flag=='nodm': dm=0*dm
        return np.concatenate([a1*U[1:], a1*acc, dm])
N=int(sys.argv[1])
for flag in (None,'noPterm','Gam0','noalpha','nodm'):
    Mut.flag=flag
    L=Mut(5e13,N=N,pert=None); ev,_=jac_eigs(L); g=ev.real.max()
    L2=Mut(2.3106026e13,N=N,pert=None); ev2,_=jac_eigs(L2); idx=np.argsort(np.abs(ev2)); nz=[i for i in idx if abs(ev2[i])>1e-3]
    print('mutation=%-8s N=%d growth(5e13)=%.5f   lowest sigma2 at rho_max=%.4f'%(flag,N,g,-(ev2[nz[0]]**2).real),flush=True)
