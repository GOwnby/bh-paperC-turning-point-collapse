import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
ev=lambda s,u:u[0]-1e7; ev.terminal=True
for S in (40,80,160,320):
    x=np.sqrt(S); u0=-x-1/(8*S**2); v0=1/(2*x)   # slow-manifold start (first-order adiabatic correction)
    sol=solve_ivp(f,(-S,60),[u0,v0],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02)
    s=sol.t[-1]; u=sol.y[0,-1]
    # also max |u + sqrt(-s)| over the slow part (oscillation residual)
    m=sol.t<-5; res=np.max(np.abs(sol.y[0][m]+np.sqrt(-sol.t[m])))
    print(S,'s*=',s+np.sqrt(6/u),'resid on slow manifold',res,flush=True)
# crossing at u_h
S=160; x=np.sqrt(S); sol=solve_ivp(f,(-S,60),[-x-1/(8*S**2),1/(2*x)],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02,dense_output=True)
for uh in (1e2,1e3,1e4,1e5):
    ss=np.linspace(0,sol.t[-1],200001); uu=sol.sol(ss)[0]; i=np.argmax(uu>uh); print('u=',uh,'s=',ss[i])
