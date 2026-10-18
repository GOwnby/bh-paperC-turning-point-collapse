import sys; sys.path.insert(0,'.')
from pulse import *
import scipy.linalg as sla
def rr(rc,nb=5,x_end=1e-3,npan=400,ng=10,Rcut=None,basis='proj'):
    sol,R,yR,xc=background(rc,x_end=x_end); sc=sigma2_scale(R,yR[0],yR[1])
    Rb=R if Rcut is None else Rcut
    # panels uniform in s with r=R(1-(1-s)^2)
    s=np.linspace(0,1,npan+1); edges=Rb*(1-(1-s)**2) if Rcut is None else None
    if Rcut is not None:
        edges=np.linspace(0,Rcut,npan+1)
    xg,wg=np.polynomial.legendre.leggauss(ng)
    A=np.zeros((nb,nb)); B=np.zeros((nb,nb))
    for i in range(npan):
        a,b=edges[i],edges[i+1]; h=b-a
        r=0.5*(a+b)+0.5*h*xg; w=0.5*h*wg
        r=np.maximum(r,1e-12)
        Y=sol.sol(np.minimum(r,sol.t[-1]))
        m,nu,x=Y[0],Y[1],Y[2]
        n1=nu1(r,m,pgeo(x))
        Pk,W,Q=coeffsA(r,m,nu,x)
        xi=[];dxi=[]
        u=(r/Rb)
        for k in range(nb):
            if k==0: f=r; df=np.ones_like(r)
            else:
                f=r*(1-u**2)**k; df=(1-u**2)**k+r*k*(1-u**2)**(k-1)*(-2*r/Rb**2)
            xi.append(f);dxi.append(df)
        zeta=[np.exp(-nu)*r**2*f for f in xi]
        dz=[np.exp(-nu)*(r**2*df+(2*r-n1*r**2)*f) for f,df in zip(xi,dxi)]
        for p in range(nb):
            for q in range(nb):
                A[p,q]+=np.sum(w*(Pk*dz[p]*dz[q]-Q*zeta[p]*zeta[q]))
                B[p,q]+=np.sum(w*W*zeta[p]*zeta[q])
    ev=sla.eigh(A,B,eigvals_only=True)
    return ev[0]*sc, R
if __name__=="__main__":
    for nb in [1,2,3,4,5,6,8,10]:
        print(nb, rr(5e13,nb=nb),flush=True)
