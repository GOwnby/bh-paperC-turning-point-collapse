import sys, os, json, time, numpy as np
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm, model_from
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c
th,state=handoff(sc,yh=0.1)
L=Lag(rm,N=N,pert=(state[0],state[1]))
fr=[0.02,0.1,0.25,0.5,0.75,0.9,1.0]
idx=[int(np.argmin(abs(L.mu/L.mu_tot-f))) for f in fr]
print('mass-fraction nodes',[(f,i,round(L.mu[i]/L.mu_tot,4)) for f,i in zip(fr,idx)])
y=L.y0.copy(); t=0.0; prev=None; cross={}
def lams(y):
    r=L.unpack(y)[0]; return np.array([r[i]/L.r0[i] for i in idx])
prev=lams(y)
while True:
    dt=L.dt_cfl(y); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    yn=y+dt/6*(k1+2*k2+2*k3+k4); cur=lams(yn)
    for v in (0.5,0.15):
        for k,f in enumerate(fr):
            if (v,f) not in cross and cur[k]<v<=prev[k]: cross[(v,f)]=t+dt*(prev[k]-v)/(prev[k]-cur[k])
    prev=cur; y=yn; t+=dt
    r,U,m=L.unpack(y)
    if np.max(2*G*m[1:]/(r[1:]*c**2))>=1.0 or t>10: break
print('N',N,'t_AH',t)
for f in fr:
    a=cross.get((0.5,f)); b=cross.get((0.15,f))
    print('  mass frac %.2f: t(0.5)=%.4f t(0.15)=%s plunge=%s'%(f,a,'%.4f'%b if b else None,'%.4f'%(b-a) if b else None))
# homologous
for scl in (sc,1.0):
    t5,t15=model_from(state,scl); print('  homologous scale %.3f t(0.5)=%.4f t(0.15)=%.4f plunge=%.4f'%(scl,t5,t15,t15-t5))
