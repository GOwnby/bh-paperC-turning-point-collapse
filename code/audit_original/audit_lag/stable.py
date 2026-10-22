import sys, os, numpy as np, time
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
G,c=bg.G,bg.c
N=int(sys.argv[1]); T=float(sys.argv[2]); amp=1e-3
L=Lag(1e13,N=N,pert=(1-amp,0.0)); jh=np.searchsorted(L.mu,0.5*L.mu_tot)
y=L.y0.copy(); t=0.0; k=0; ts=[];ls=[];cons=[]
while t<T:
    dt=L.dt_cfl(y,0.4); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
    if k%20==0:
        r,U,m=L.unpack(y); Gam,rho,P,e=L.zone(r,U,m); V=4*np.pi/3*(r[1:]**3-r[:-1]**3)
        ts.append(t); ls.append(r[jh]/L.r0[jh]); cons.append(abs(m[-1]-np.sum(e*V)/c**2)/m[-1])
ts=np.array(ts); ls=np.array(ls)-1
# omega from upward zero crossings of lam_h-<lam>
x=ls-ls.mean(); zc=[ts[i]-x[i]*(ts[i+1]-ts[i])/(x[i+1]-x[i]) for i in range(len(x)-1) if x[i]<0<=x[i+1]]
per=np.diff(zc); print('N=%d steps=%d period(s)=%s -> omega=%s  (Jacobian N=200: 2.84199; RR 2.80176)'%(N,k,per,2*np.pi/per))
print(' amplitude start %.3e, max %.3e, min %.3e; constraint |m_N - sum eV/c^2|/M: start %.1e max %.1e end %.1e'%(ls[0],ls.max(),ls.min(),cons[0],max(cons),cons[-1]))
