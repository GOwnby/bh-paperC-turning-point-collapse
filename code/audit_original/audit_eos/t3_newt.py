import numpy as np
from scipy.integrate import solve_ivp
from t2_tov import *
hbar=h/2/np.pi
omega3=2.01824
MCh=np.sqrt(3*np.pi)/2*omega3*(hbar*c/G)**1.5/(mue*mu)**2
print("Newtonian Chandrasekhar mass (n=3 Lane-Emden, mu_e=2, m_u): %.5f Msun(code)  %.5f Msun(GM)  -> 5.83/mu_e^2 form: %.4f"%(MCh/Msun_code,MCh/Msun_GM, MCh/Msun_GM*mue**2))
# Newtonian structure with full Chandrasekhar EOS
def newt(rhoc):
    xc=(rhoc/(mue*mu*nfac))**(1/3); Pc=Pf(xc)
    def rhs(r,y):
        x,m=y; rho=mue*mu*nfac*x**3
        return [-G*rho*m/r**2/dPdx(x), 4*np.pi*r*r*rho]
    ev=lambda r,y:Pf(y[0])-1e-10*Pc; ev.terminal=True
    s=solve_ivp(rhs,[1,5e8],[xc,4/3*np.pi*rhoc],events=ev,rtol=1e-12,atol=[1e-14,1e-6],method='DOP853')
    return s.y_events[0][0][1]/Msun_GM, s.t_events[0][0]
for rc in (1e11,1e12,1e13,2.3e13,5e13,1e14,1e15,1e17):
    print(" Newtonian rho_c=%.1e : M=%.4f Msun  R=%.1f km"%(rc,*[newt(rc)[0],newt(rc)[1]/1e3]))
