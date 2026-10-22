import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ai_zeros
a1=ai_zeros(1)[0][0]
def pole(nu,thr=1e8):
    S=max(60,40+60*nu**2)
    x=np.sqrt(S); u0=-x; up0=1/(2*x); upp0=(1/4)*x**-3
    d=(nu*up0+upp0)/(2*u0); u=u0+d; v=up0
    f=lambda s,y:[y[1],y[0]**2+s-nu*y[1]]
    ev=lambda s,y:y[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,200),[u,v],method='Radau' if nu>5 else 'DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.02 if nu<5 else 0.2)
    s=sol.t[-1]; uu=sol.y[0,-1]; xx=np.sqrt(6/uu)
    # u=6/x^2-(6nu/5)/x : solve for x given u
    # 6/x^2 - 1.2 nu/x = uu -> uu x^2 + 1.2 nu x - 6 =0
    xr=(-1.2*nu+np.sqrt((1.2*nu)**2+24*uu))/(2*uu)
    return s+xr, s+xx
for nu in (0,0.03,0.1,0.3,1,3,10):
    s,s2=pole(nu); asy=-a1*nu**(2/3)
    print(nu,'%.5f %.5f'%(s,s2),'airy %.4f'%asy,'ratio %.3f'%(s/asy if nu>0 else 0), 'old(table) shift', flush=True)
