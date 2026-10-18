import sys; sys.path.insert(0,'.')
from run_shoot import *
rc=5e13
def bracket(rc,x_end=1e-3,r0=1e-6,rtol=1e-12,a=-100,b=100):
    sol,R,yR,xc=background(rc,x_end=x_end,r0=r0,rtol=rtol)
    sc=sigma2_scale(R,yR[0],yR[1])
    f=lambda s2: shootA(rc,s2/sc,x_end=x_end,r0=r0,rtol=rtol)
    return brentq(f,a,b,xtol=1e-10,rtol=1e-13)
for r0 in [1e-4,1e-5,1e-7]:
    print("r0",r0,bracket(rc,r0=r0),flush=True)
for rt in [1e-9,1e-10,1e-13]:
    print("rtol",rt,bracket(rc,rtol=rt),flush=True)
print("x_end 1e-5",bracket(rc,x_end=1e-5),flush=True)
