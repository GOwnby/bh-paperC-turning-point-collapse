import sys; sys.path.insert(0,'.')
from pulseB import *
def rootA(rc,g=(-3e3,-300,-30,-3,0,3,30,300,3e3,1e4)):
    sol,R,yR,xc=background(rc,r0=1e-6); sc=sigma2_scale(R,yR[0],yR[1])
    f=lambda s: shootA(rc,s/sc)
    v=[f(s) for s in g]; s0=np.sign(v[0])
    for i in range(len(g)-1):
        if s0*v[i]>0 and s0*v[i+1]<=0:
            return brentq(lambda s:s0*f(s),g[i],g[i+1],xtol=1e-12), 2*yR[0]/R
for rc in [1e8,1e9,1e10]:
    print(rc,rootA(rc))
