import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import airy, ai_zeros
print('Ai zero1',ai_zeros(1)[0][0], 'Ai(a1)',airy(ai_zeros(1)[0][0])[0])
def pole_pi(S,method,rtol,thr=1e6,start='first'):
    x=np.sqrt(S)
    if start=='first': u0=-x-1/(8*S**2); v0=1/(2*x)+ (-4*(-1/(8*S**2))*0) 
    f=lambda s,u:[u[1],u[0]**2+s]
    ev=lambda s,u:u[0]-thr; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[u0,v0],method=method,events=ev,rtol=rtol,atol=rtol*1e-2,max_step=0.05)
    s=sol.t[-1]; u=sol.y[0,-1]
    return s-np.sqrt(6/u)
for S,m,r in ((25,'LSODA',1e-10),(50,'Radau',1e-11),(100,'LSODA',1e-10),(500,'DOP853',1e-10),(1000,'LSODA',1e-9)):
    print(S,m,r,pole_pi(S,m,r),flush=True)
# threshold dependence
for thr in (1e3,1e4,1e5,1e7,1e8):
    print('thr',thr,pole_pi(100,'LSODA',1e-11,thr))
# raw oscillation sensitivity: rest start
def pole_rest(S):
    x=np.sqrt(S); f=lambda s,u:[u[1],u[0]**2+s]
    ev=lambda s,u:u[0]-1e6; ev.terminal=True
    sol=solve_ivp(f,(-S,60),[-x,0],method='DOP853',events=ev,rtol=1e-11,atol=1e-12,max_step=0.05)
    return sol.t[-1]-np.sqrt(6/sol.y[0,-1])
print('rest start',[ (S,pole_rest(S)) for S in (25,100,400)])
# Riccati/Airy
a1=ai_zeros(1)[0][0]
for T in (200,):
    tau0=-T; w0=-np.sqrt(T)+0.0
    # adiabatic: w=-sqrt(-tau)-1/(4 tau)?? derive: w'=w^2+tau; w=-x+d, x=sqrt(-tau): w'=1/(2x) ; 1/(2x)=x^2... => -2xd = 1/(2x) -> d=-1/(4x^2)
    x=np.sqrt(T); w0=-x-1/(4*x**2)
    f=lambda t,w:w**2+t
    ev=lambda t,w:w[0]-1e7; ev.terminal=True
    sol=solve_ivp(f,(-T,10),[w0],method='LSODA',events=ev,rtol=1e-12,atol=1e-12,max_step=0.05)
    print('Riccati pole',sol.t[-1]-1/np.sqrt(sol.y[0,-1]) if False else sol.t[-1]+1/sol.y[0,-1],' -a1=',-a1)
# check w=Ai'(-tau)/Ai(-tau) vs -sqrt(-tau)
for tau in (-50,-10,-2,0,1,2):
    Ai,Aip,_,_=airy(-tau); w=Aip/Ai  # w=-phi'/phi with phi=Ai(-tau): phi'=-Ai'(-tau)
    print(tau,'w',w,'-sqrt(-tau)',-np.sqrt(-tau) if tau<0 else None)
# numeric check of w'=w^2+tau
tau=0.7; h=1e-5
W=lambda t:(lambda A:A[1]/A[0])(airy(-t))
print('w\'',(W(tau+h)-W(tau-h))/(2*h),'w^2+tau',W(tau)**2+tau)
