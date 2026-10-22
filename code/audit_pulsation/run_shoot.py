import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulse import *
import numpy as np
def sig2_A(rc,x_end=1e-3,r0=1e-6,rtol=1e-12,lo=None,hi=None):
    sol,R,yR,xc=background(rc,x_end=x_end,r0=r0,rtol=rtol)
    M=yR[0]; sc=sigma2_scale(R,M,yR[1])
    f=lambda s2: shootA(rc,s2/sc,x_end=x_end,r0=r0,rtol=rtol)
    # scan upward from very negative to find first sign change from + to -
    grid=-np.logspace(5,-2,40)  # negative
    grid=np.concatenate([grid,[0.0],np.logspace(-2,5,60)])
    prev=None
    for g in grid:
        v=f(g)
        if prev is not None and prev[1]>0 and v<=0:
            s2=brentq(f,prev[0],g,xtol=1e-10,rtol=1e-13)
            return s2,R*L/1e3,M*L*c**2/G/Msun
        prev=(g,v)
    return None
if __name__=="__main__":
    for rc in [float(a) for a in sys.argv[1:]]:
        print(rc, sig2_A(rc))
