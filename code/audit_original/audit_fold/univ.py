import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
S=100;x=np.sqrt(S)
kap=18.733;a=0.004027;yr=3.15576e7;Md=1e-7/yr
be=(4*a/(kap**2*Md))**0.2; al=2/(kap*be**2); uh=0.01/al
print('beta',be,'alpha',al,'u at y=1e-2',uh)
sol=solve_ivp(f,(-S,3.5),[-x-1/(8*S**2),1/(2*x)],method='DOP853',rtol=1e-12,atol=1e-13,max_step=0.01,dense_output=True)
ss=np.linspace(2.5,3.41,900001);uu=sol.sol(ss)[0];i=np.argmax(uu>uh);print('universal s at u=uh:',ss[i],'; model(rescaled) 3.33666-3.33685; raw 3.4373?')
# in terms of lambda=0.5 equivalent? y=lam^-3-1 =7 -> u=7/al
for yv in (7.0,):
    u7=yv/al; i=np.argmax(uu>u7); print('normal-form s at y=7 (lam=0.5)',ss[i],'-> t=',ss[i]*be,' (table 366.4); u=',u7)
