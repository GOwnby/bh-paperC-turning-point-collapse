import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
ev=lambda s,u:u[0]-1e8; ev.terminal=True
for S in (40,80,160,320):
    x=np.sqrt(S); u0=-x-1/(8*S**2); v0=1/(2*x)
    sol=solve_ivp(f,(-S,60),[u0,v0],method='DOP853',events=ev,rtol=1e-13,atol=1e-13,max_step=0.02,dense_output=True)
    s=sol.t[-1]; u=sol.y[0,-1]
    print(S,'s*=%.6f'%(s-np.sqrt(6/u)))
    if S==160:
        ss=np.linspace(0,sol.t[-1],400001); uu=sol.sol(ss)[0]
        # u at which? y=0.01 for Mdot=1e-7
        kap=18.733;a=0.004027;yr=365.25*86400
        md=1e-7/yr; beta=(4*a/(kap**2*md))**0.2; alpha=2/(kap*beta**2)
        print('beta',beta,'alpha',alpha,'u(y=0.01)=',0.01/alpha)
        uh=0.01/alpha; i=np.argmax(uu>uh); print('s at u_h',ss[i], 'resc model', 359.348/beta, 'diff %', (359.348/beta/ss[i]-1)*100)
        # with beta using kappa? y def: y=rho/rho_max-1 ; model uses (lam^-3 -1)
