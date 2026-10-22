import sys, os, numpy as np, io, contextlib, time
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
G,c=bg.G,bg.c
N=int(sys.argv[1]); th,state=handoff(sc,yh=0.1)
for cfl in [float(x) for x in sys.argv[2:]]:
    L=Lag(rm,N=N,pert=(state[0],state[1])); jh=np.searchsorted(L.mu,0.5*L.mu_tot)
    y=L.y0.copy(); t=0.0; lp=1.0; t5=None; gmax=0; ok=True; k=0; w0=time.time()
    while t<2.1:
        dt=L.dt_cfl(y,cfl); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
        y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
        r,U,m=L.unpack(y)
        if not np.isfinite(y).all() or np.min(np.diff(r))<=0: ok=False; break
        lam=r[jh]/L.r0[jh]
        if t5 is None and lam<0.5: t5=t-dt*(0.5-lam)/(lp-lam)
        lp=lam
        if k%50==0:
            Gam=np.sqrt(np.maximum(1+U[1:]**2/c**2-2*G*m[1:]/(r[1:]*c**2),0)); gmax=max(gmax,Gam.max())
        if np.max(2*G*m[1:]/(r[1:]*c**2))>=1.0: break
    print('N=%d cfl=%.2f ok=%s steps=%d t_end=%.5f t(0.5)=%s max Gamma(sampled)=%.4f wall %.0fs'%(N,cfl,ok,k,t,t5,gmax,time.time()-w0),flush=True)
