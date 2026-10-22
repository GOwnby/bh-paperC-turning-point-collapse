import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
def run(S,thr,rtol,meth='DOP853',start='x'):
    ev=lambda s,u:u[0]-thr; ev.terminal=True
    x=np.sqrt(S); u0=-x-1/(8*S**2); v0=1/(2*x)
    sol=solve_ivp(f,(-S,60),[u0,v0],method=meth,events=ev,rtol=rtol,atol=rtol,max_step=0.02)
    s=sol.t[-1]; u=sol.y[0,-1]
    return s-np.sqrt(6/u), s, u
for thr in (1e3,1e4,1e5,1e6,1e7,1e8):
    for rtol in (1e-10,1e-13):
        print(thr,rtol,run(60,thr,rtol))
print('Radau',run(60,1e6,1e-12,'Radau'))
print('DOP853 S=60 w/ x start only', )
