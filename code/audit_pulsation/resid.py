import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulseB import *
rc=5e13
sol0,R,yR,xc=background(rc); sc=sigma2_scale(R,yR[0],yR[1])
s2=-28.555947372879654; w2=s2/sc
F,R,yR,sol=shootA(rc,w2,full=True)
fB=rhsB(w2)
def XdP(r):
    m,nu,x,z,eta=sol.sol(r)
    e=egeo(x);Pg=pgeo(x);G1=Gamma1(x)
    e2l=1/(1-2*m/r); el=np.sqrt(e2l)
    X=np.exp(nu)*z/r**2
    dP=-np.exp(-np.log(el)-2*nu)*eta-G1*Pg*4*np.pi*r*(e+Pg)*e2l*X
    return X,dP
for fr in [0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95]:
    r=fr*R; h=1e-5*R
    m,nu,x,z,eta=sol.sol(r)
    X,dP=XdP(r)
    Xp=(XdP(r+h)[0]-XdP(r-h)[0])/(2*h); dPp=(XdP(r+h)[1]-XdP(r-h)[1])/(2*h)
    d=fB(r,[m,nu,x,X,dP])
    print(fr, "Xp res",(d[3]-Xp)/abs(Xp), " dPp: B=",d[4]," A=",dPp, " rel diff",(d[4]-dPp)/abs(dPp))
