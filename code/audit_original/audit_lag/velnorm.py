import sys, os, numpy as np, io, contextlib
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
G,c=bg.G,bg.c
N=int(sys.argv[1]); mode=sys.argv[2]
th,state=handoff(sc,yh=0.1); lam,ld=state
L=Lag(rm,N=N,pert=(lam,ld))
r=L.r0*lam
# lapse at the unperturbed equilibrium nodes
r0=L.r0.copy(); U0=np.zeros_like(r0); m0=L.solve_m(r0,U0)
Gam,rho,P,e=L.zone(r0,U0,m0); h=(e+P)/rho
alpha_s=(1-2*G*m0[-1]/(r0[-1]*c**2))/Gam[-1]
az=alpha_s*c**2/h; aj=np.empty(N+1); aj[1:N]=0.5*(az[1:]+az[:-1]); aj[N]=alpha_s; aj[0]=az[0]
print('alpha_s %.5f alpha_center %.5f'%(alpha_s,aj[1]))
if mode=='div': U=ld*L.r0/aj
elif mode=='mul': U=ld*L.r0*aj
else: U=ld*L.r0
U[0]=0
m=L.solve_m(r,U); L.y0=np.concatenate([r[1:],U[1:],m[1:]])
jh=np.searchsorted(L.mu,0.5*L.mu_tot)
y=L.y0.copy(); t=0.0; lp=lam; t5=None
while True:
    dt=L.dt_cfl(y,0.8); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt
    rr=L.unpack(y)[0]; l=rr[jh]/L.r0[jh]
    if l<0.5: t5=t-dt*(0.5-l)/(lp-l); break
    lp=l
print('N=%d mode=%s t(0.5)=%.5f'%(N,mode,t5))
