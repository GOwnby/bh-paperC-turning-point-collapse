# P1 finite-element generalized eigenproblem on the Sturm-Liouville functional
#   int (Pk zeta'^2 - Q zeta^2) dr = w2 int W zeta^2 dr ,  zeta(0)=0 essential, surface natural.
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulse import *
import scipy.linalg as sla
def fem(rc,N=400,q=2.0,x_end=1e-4,ng=6):
    sol,R,yR,xc=background(rc,x_end=x_end)
    M=yR[0]; sc=sigma2_scale(R,M,yR[1])
    s=np.linspace(0,1,N+1)
    r=R*(1-(1-s)**q)       # clusters nodes near surface
    r[0]=0.0
    xg,wg=np.polynomial.legendre.leggauss(ng)
    K=np.zeros((N+1,N+1)); Mm=np.zeros((N+1,N+1))
    for i in range(N):
        a,b=r[i],r[i+1]; h=b-a
        rr=0.5*(a+b)+0.5*h*xg; ww=0.5*h*wg
        rr=np.maximum(rr,1e-12)
        Y=sol.sol(np.minimum(rr,sol.t[-1]))
        Pk,W,Q=coeffsA(rr,Y[0],Y[1],Y[2])
        phi0=(b-rr)/h; phi1=(rr-a)/h; d0=-1/h; d1=1/h
        loc=[(phi0,d0),(phi1,d1)]
        for p in range(2):
            for qq in range(2):
                K[i+p,i+qq]+=np.sum(ww*(Pk*loc[p][1]*loc[qq][1]-Q*loc[p][0]*loc[qq][0]))
                Mm[i+p,i+qq]+=np.sum(ww*W*loc[p][0]*loc[qq][0])
    K=K[1:,1:]; Mm=Mm[1:,1:]   # drop node 0 (zeta(0)=0)
    ev=sla.eigh(K,Mm,eigvals_only=True,subset_by_index=[0,2])
    return ev*sc
if __name__=="__main__":
    rc=float(sys.argv[1])
    for N in [100,200,400,800,1600]:
        print(N,fem(rc,N=N),flush=True)
