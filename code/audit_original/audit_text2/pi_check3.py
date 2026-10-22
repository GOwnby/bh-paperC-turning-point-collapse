import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
S=160; x=np.sqrt(S)
for thr in (1e4,1e6,1e8):
    ev=lambda s,u:u[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[-x-1/(8*S**2),1/(2*x)],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02)
    s=sol.t[-1]; u=sol.y[0,-1]
    print(thr,'s*(correct, +)=%.6f'%(s+np.sqrt(6/u)),' s*(code, -)=%.6f'%(s-np.sqrt(6/u)))
# Laurent-corrected: u=6/x^2 - s*/10 x^2 - x^3/6
# verify via u*(s*-s)^2 at s*=3.41167
sstar=3.411676
for thr in (1e2,1e3,1e4):
    ev=lambda s,u:u[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[-x-1/(8*S**2),1/(2*x)],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02)
    s=sol.t[-1]; u=sol.y[0,-1]; print(thr,u*(sstar-s)**2)
# also nonadiabatic starts check
for S in (40,80,320):
    x=np.sqrt(S);ev=lambda s,u:u[0]-1e6; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[-x-1/(8*S**2),1/(2*x)],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02)
    print(S,sol.t[-1]+np.sqrt(6/sol.y[0,-1]))
