import numpy as np
from scipy.integrate import solve_ivp
f=lambda s,u:[u[1],u[0]**2+s]
B=6**-0.2; A=6*B**2
print('B',B,'A',A,' u(0) expected from y(0)=-0.18755430834:',A*-0.1875543083404949,' u\'(0) expected', A*B*0.3049055602612289)
for S,m,rt,thr in ((60,'DOP853',1e-13,1e8),(200,'DOP853',1e-12,1e9),(600,'DOP853',1e-12,1e7),(30,'DOP853',1e-13,1e8)):
    x=np.sqrt(S)
    ev=lambda s,u:u[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,3.5),[-x-1/(8*S**2),1/(2*x)],method=m,rtol=rt,atol=1e-14,max_step=0.02,events=ev,dense_output=True)
    s=sol.t[-1];u=sol.y[0,-1]
    print(S,thr,'s*=',s+np.sqrt(6/u),'(author formula s-sqrt(6/u)=',s-np.sqrt(6/u),') u(0)=',sol.sol(0.0)[0],"u'(0)=",sol.sol(0.0)[1],flush=True)
    print('   x-pole in y-variable',(s+np.sqrt(6/u))*B)
