import numpy as np
from scipy.integrate import solve_ivp
S=100; x=np.sqrt(S)
f=lambda s,u:[u[1],u[0]**2+s]
for m,rt in (('DOP853',1e-13),('DOP853',1e-10),('LSODA',1e-11),('Radau',1e-12)):
    sol=solve_ivp(f,(-S,3.45),[-x-1/(8*S**2),1/(2*x)],method=m,rtol=rt,atol=1e-14,max_step=0.02,dense_output=True)
    out=[]
    for ut in (1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8):
        ss=np.linspace(2.0,sol.t[-1],2000001); uu=sol.sol(ss)[0]; i=np.argmax(uu>ut)
        if uu[i]>ut and i>0: out.append((ut,ss[i]-np.sqrt(6/uu[i])))
    print(m,rt,'tend',sol.t[-1],[ (u,round(v,6)) for u,v in out],flush=True)
