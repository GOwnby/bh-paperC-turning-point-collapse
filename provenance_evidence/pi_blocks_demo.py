import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
print("OLD block (original fold_dyn2.py):")
ev=lambda s,u:u[0]-1e6; ev.terminal=True
for S in (20,30,45,60):
    sol=solve_ivp(f,(-S,60),[-np.sqrt(S),-0.5/np.sqrt(S)],method='DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.02)
    s_end=sol.t[-1]; u_end=sol.y[0,-1]; print(S, s_end+np.sqrt(6/u_end))
print("NEW block (modified fold_dyn2.py):")
ev=lambda s,u:u[0]-1e7; ev.terminal=True
for S in (40,80,160,320):
    x=np.sqrt(S); u0=-x-1/(8*S**2); v0=1/(2*x)
    sol=solve_ivp(f,(-S,60),[u0,v0],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02)
    s_end=sol.t[-1]; u_end=sol.y[0,-1]; print(S, s_end+np.sqrt(6/u_end))
