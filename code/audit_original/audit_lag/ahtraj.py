import sys, os, json, time, numpy as np
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c
th,state=handoff(sc,yh=0.1)
L=Lag(rm,N=N,pert=(state[0],state[1]))
y=L.y0.copy(); t=0.0; rows=[]
while True:
    r,U,m=L.unpack(y); comp=2*G*m[1:]/(r[1:]*c**2); jm=int(np.argmax(comp))
    if comp.max()>0.03: rows.append([t,comp.max(),jm+1,m[jm+1]/m[-1],L.mu[jm+1]/L.mu_tot,r[jm+1]/L.r0[jm+1],U[jm+1]/c, m[1]/m[-1]])
    if comp.max()>=1.02: break
    dt=L.dt_cfl(y); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt
rows=np.array(rows); np.save('/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/audit_lag/ahtraj_N%d.npy'%N,rows)
print('N',N,'rows',len(rows))
for th_ in (0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.8,0.9,0.95,0.98,1.0):
    i=np.argmax(rows[:,1]>=th_); r_=rows[i]
    print('comp>=%.2f t=%.5f comp=%.4f jm=%d m/M=%.4f mu/mu_tot=%.4f r/r0=%.5f U/c=%.3f  m(j=1)/M=%.2e'%(th_,*r_))
