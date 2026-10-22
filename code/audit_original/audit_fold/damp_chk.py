import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ai_zeros
a1=ai_zeros(1)[0][0]
def pole(nu,S,start,method,rt,thr=1e8):
    x=np.sqrt(S); u0=-x
    if start=='first':
        d=(1/(4*x**3)+nu/(2*x))/(2*u0); u=u0+d; v=1/(2*x)
    else: u=u0; v=0.0
    f=lambda s,y:[y[1],y[0]**2+s-nu*y[1]]
    ev=lambda s,y:y[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,400),[u,v],method=method,events=ev,rtol=rt,atol=rt*1e-3,max_step=0.1)
    return sol.t[-1]+np.sqrt(6/sol.y[0,-1])
print('nu  first(S1)  first(S2)  rest   asym ratio')
for nu in (0,0.1,0.3,1,3,10):
    S1=max(100,60*nu**2+50); S2=3*S1
    p1=pole(nu,S1,'first','LSODA',1e-10); p2=pole(nu,S2,'first','Radau',1e-10)
    pr=pole(nu,S1+300,'rest','LSODA',1e-10) if nu>=0.3 else np.nan
    print(nu,'%.5f %.5f %.5f'%(p1,p2,pr),'airy %.4f'%(-a1*nu**(2/3)),'ratio %.4f'%(p1/(-a1*nu**(2/3)) if nu else 0),flush=True)
