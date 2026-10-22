import sys, os, numpy as np, io, contextlib
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c
th,state=handoff(sc,yh=0.1); L=Lag(rm,N=N,pert=(state[0],state[1]))
y=L.y0.copy(); t=0.0
while True:
    dt=L.dt_cfl(y); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; r,U,m=L.unpack(y)
    if np.max(2*G*m[1:]/(r[1:]*c**2))>=1.0: break
f=[0.05,0.2,0.5,0.8,0.95,0.99,1.0]
print('N',N,'shell r/r0 at AH:',[(x,round(float(np.interp(x*L.mu_tot,L.mu,r)/np.interp(x*L.mu_tot,L.mu,L.r0)),4)) for x in f])
Gam,rho,P,e=L.zone(r,U,m); print(' rho_c/rho_c0 (zone0)=%.3e ; central rest-mass density in g/cc %.2e'%(rho[0]/L.rho_c,rho[0]*1e-3))
print(' minGam2',(1+U[1:]**2/c**2-2*G*m[1:]/(r[1:]*c**2)).min())
