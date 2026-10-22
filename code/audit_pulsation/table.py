import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulseB import *
def root(shoot,rc,a=-300,b=300):
    sol,R,yR,xc=background(rc); sc=sigma2_scale(R,yR[0],yR[1])
    f=lambda s2: shoot(rc,s2/sc)
    # find bracket by scanning from -1e5 upward for + to -
    g=[-3e3,-3e3,-1e3,-300,-100,-30,-10,-3,-1,0,1,3,10,30,100,300,1e3,3e3,1e4]
    vals=[f(s) for s in g]
    s0=np.sign(vals[0]); vals=[v*s0 for v in vals]
    for i in range(len(g)-1):
        if vals[i]>0 and vals[i+1]<=0:
            return brentq(lambda s:s0*f(s),g[i],g[i+1],xtol=1e-11,rtol=1e-14), R*L/1e3, yR[0]*L*c**2/G/Msun, sc
    return None
print("rho_c   sigma2_A   sigma2_B   R_km  M/Msun")
for rc in [1e13,2e13,2.31e13,3e13,5e13]:
    a=root(shootA,rc); b=root(shootB,rc)
    # node count of zeta for A at root
    F,R,yR,sol=shootA(rc,a[0]/a[3],full=True)
    zs=sol.y[3]; nodes=int(np.sum(zs[1:]*zs[:-1]<0))
    print("%.4g  %.8g  %.8g  %.3f  %.5f  nodes(zeta)=%d"%(rc,a[0],b[0],a[1],a[2],nodes),flush=True)
# zero crossing in rho_c
def g(lr):
    return shootA(10**lr,0.0)
lz=brentq(g,13.30,13.45,xtol=1e-13)
print("zero crossing rho_c (A):",10**lz)
def gB(lr): return shootB(10**lr,0.0)
lzB=brentq(gB,13.30,13.45,xtol=1e-13)
print("zero crossing rho_c (B):",10**lzB)
