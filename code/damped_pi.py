import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.special import ai_zeros
a1=ai_zeros(1)[0][0]; print('first Airy zero',a1)
def pole(nu,S=None):
    S=S or max(60,40+60*nu**2)    # start far enough that the adiabatic start is accurate
    x=np.sqrt(S); u0=-x; up0=1/(2*x); upp0=(1/4)*x**-3
    d=(nu*up0+upp0)/(2*u0); u=u0+d; v=up0     # first-order adiabatic correction
    f=lambda s,y:[y[1],y[0]**2+s-nu*y[1]]
    ev=lambda s,y:y[0]-1e7; ev.terminal=True
    sol=solve_ivp(f,(-S,200),[u,v],method='Radau' if nu>5 else 'DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.02 if nu<5 else 0.2)
    return sol.t[-1]+np.sqrt(6/sol.y[0,-1])
out=[]
for nu in (0,0.1,0.3,1,3,10):
    s=pole(nu); asy=-a1*nu**(2/3) if nu>0 else np.nan
    out.append((nu,s,asy)); print(nu,s,asy,s/asy if nu>0 else '',flush=True)
json.dump(out,open('damped_pi.json','w'))
