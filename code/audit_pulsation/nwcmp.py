import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulseB import *
import newton
def rootgen(f,scale):
    g=[scale*u for u in (-30,-10,-3,-1,-0.3,-0.1,0,0.1,0.3,1,3,10,30)]
    v=[f(s) for s in g]; s0=np.sign(v[0])
    for i in range(len(g)-1):
        if s0*v[i]>0 and s0*v[i+1]<=0:
            return brentq(lambda s:s0*f(s),g[i],g[i+1],xtol=1e-13,rtol=1e-12)
for rc in [1e8,1e9,1e10,1e11]:
    sc0=G*rc
    sol,R,yR,xc=background(rc,x_end=1e-3); sc=sigma2_scale(R,yR[0],yR[1])
    a=rootgen(lambda s:shootA(rc,s/sc),sc0)
    b=rootgen(lambda s:shootB(rc,s/sc),sc0)
    n=rootgen(lambda s:newton.run(rc,s)[0],sc0)
    print("rho_c %.1e  compactness 2M/R=%.3e  GR_A=%.8g GR_B=%.8g Newton=%.8g  A/N-1=%.3e"%(rc,2*yR[0]/R,a,b,n,a/n-1),flush=True)
