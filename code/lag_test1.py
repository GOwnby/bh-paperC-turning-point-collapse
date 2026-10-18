import numpy as np, time, sys
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
N=int(sys.argv[1]); rc=float(sys.argv[2]); T=float(sys.argv[3]); amp=float(sys.argv[4]) if len(sys.argv)>4 else 1e-4
L=Lag(rc,N=N,pert=(1-amp,0.0))
d=build_RR_matrices(L.bgd,5); s2=eigh(d[0],d[1])[0][0]; print('RR sigma^2',s2,flush=True)
y=L.y0.copy(); t=0; k=0; rows=[]; t0=time.time(); nxt=0
while t<T:
    dt=L.dt_cfl(y)
    k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
    if t>=nxt:
        r,U,m=L.unpack(y); jh=np.searchsorted(L.mu,0.5*L.mu_tot); jq=np.searchsorted(L.mu,0.25*L.mu_tot)
        rows.append([t,r[jh]/L.r0[jh],r[jq]/L.r0[jq],U[jh]]); nxt+=T/400
        if len(rows)%40==0: print(k,t,dt,rows[-1],time.time()-t0,flush=True)
np.save(f'lag_t1_{N}_{rc:.1e}.npy',np.array(rows))
