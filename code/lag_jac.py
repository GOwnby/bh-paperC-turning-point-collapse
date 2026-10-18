import numpy as np, background as bg; bg.mu_e=2.0
from lag_gr import Lag
def jac_eigs(L):
    N=L.N; y0=L.y0.copy(); n=len(y0); J=np.zeros((n,n))
    for k in range(n):
        h=1e-7*max(abs(y0[k]),1e-30) if (k<N or k>=2*N) else 1e-3
        yp=y0.copy(); yp[k]+=h; ym=y0.copy(); ym[k]-=h
        J[:,k]=(L.f(0,yp)-L.f(0,ym))/(2*h)
    return np.linalg.eig(J)
if __name__=="__main__":
    L=Lag(1e13,N=100,pert=None); ev,V=jac_eigs(L); k=np.argmax(ev.real); v=V[:,k]; N=100
    print(ev[k]); 
    for nm,seg in (('r',v[:N]),('U',v[N:2*N]),('m',v[2*N:])): print(nm,np.abs(seg[-6:])/np.max(np.abs(seg)))
    Gam,rho,P,e=L.zone(*L.unpack(L.y0)); print('rho last',rho[-4:],'P last',P[-4:])
    r=L.unpack(L.y0)[0]; print('dr last',np.diff(r)[-4:])
