import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
S=100;x=np.sqrt(S)
for A in (-1,-0.3,-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3,1):
    ev=lambda s,u:u[0]-1e7; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[-x-1/(8*S**2)+A,1/(2*x)],method='DOP853',events=ev,rtol=1e-11,atol=1e-12,max_step=0.05,dense_output=True)
    s=sol.t[-1]+np.sqrt(6/sol.y[0,-1])
    m=(sol.t>-6)&(sol.t<-2); amp=np.max(np.abs(sol.y[0][m]+np.sqrt(-sol.t[m])))-0  # includes slow offset
    print('A(u at s=-100)=',A,'s*=%.4f'%s,' dev from slow manifold in -6<s<-2: %.3f'%amp,flush=True)
