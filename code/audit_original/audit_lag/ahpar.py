import sys, os, numpy as np, io, contextlib
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c
th,state=handoff(sc,yh=0.1); L=Lag(rm,N=N,pert=(state[0],state[1]))
def peak(y):
    r,U,m=L.unpack(y); comp=2*G*m[1:]/(r[1:]*c**2); j=int(np.argmax(comp))
    if j==0 or j==N-1: return comp[j],m[j+1]/m[-1],j+1
    x=np.array([m[j],m[j+1],m[j+2]])/m[-1]; yv=comp[j-1:j+2]
    a,b,c0=np.polyfit(x,yv,2)
    if a>=0: return comp[j],m[j+1]/m[-1],j+1
    xs=-b/(2*a); xs=min(max(xs,x[0]),x[2]); return np.polyval([a,b,c0],xs),xs,j+1
y=L.y0.copy(); t=0.0; hist=[]
while True:
    pk=peak(y); hist.append((t,)+pk)
    if pk[0]>=1.02: break
    dt=L.dt_cfl(y); k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt
h=np.array(hist); i=np.argmax(h[:,1]>=1.0)
f=(1-h[i-1,1])/(h[i,1]-h[i-1,1]); mm=h[i-1,2]+f*(h[i,2]-h[i-1,2]); tt=h[i-1,0]+f*(h[i,0]-h[i-1,0])
print('N=%d: parabolic-peak AH at t=%.5f  m_AH/M=%.4f  (prev step comp*=%.4f m=%.4f, next comp*=%.4f m=%.4f)'%(N,tt,mm,h[i-1,1],h[i-1,2],h[i,1],h[i,2]))
print('  last 6 steps (t, comp*, m*/M, node):'); 
for row in h[-6:]: print('   %.7f %.4f %.4f %d'%tuple(row))
