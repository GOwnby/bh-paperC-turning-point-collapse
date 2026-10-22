import sys, numpy as np, time
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c; Msun=bg.M_sun
th,state=handoff(sc,yh=0.1)
L=Lag(rm,N=N,pert=(state[0],state[1]))
y=L.y0.copy(); t=0.0; thr=[0.9,0.97,0.99,1.0]; done=set()
def snap(tag):
    r,U,m=L.unpack(y); Gam,rho,P,e=L.zone(r,U,m)
    comp=2*G*m[1:]/(r[1:]*c**2); j=np.argmax(comp)+1
    fr=[0.05,0.2,0.5,0.8,0.95,0.99]
    rr=[np.interp(f*L.mu_tot,L.mu,np.where(L.r0>0,r/np.where(L.r0>0,L.r0,1),0)) for f in fr]
    # parabolic-refined max of compactness vs node mass
    print(N,tag,'t=%.5f'%t,'maxcomp=%.4f'%comp.max(),'node',j,'m/M=%.4f'%(m[j]/m[-1]),'rho_c/rho_max=%.3e'%(rho[0]/rm),'radii',np.round(rr,4),'surf %.4f'%(r[-1]/L.r0[-1]),flush=True)
while True:
    r,U,m=L.unpack(y); mc=(2*G*m[1:]/(r[1:]*c**2)).max()
    for th_ in thr:
        if mc>=th_ and th_ not in done: done.add(th_); snap('comp>=%.2f'%th_)
    if mc>=1.0 or t>40: break
    dt=L.dt_cfl(y)
    k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt
