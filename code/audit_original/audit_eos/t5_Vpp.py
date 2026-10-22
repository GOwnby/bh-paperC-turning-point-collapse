import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
from exactV import Star
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
st=Star(5e13)
h=2e-3
Vpp=(st.M_of_lam(1+h)-2*st.M_of_lam(1.0)+st.M_of_lam(1-h))/h**2*c**2
# Richardson
h2=4e-3
Vpp2=(st.M_of_lam(1+h2)-2*st.M_of_lam(1.0)+st.M_of_lam(1-h2))/h2**2*c**2
Vpp_R=(4*Vpp-Vpp2)/3
bgd=build_metric(5e13)
r=bgd['r'];a=bgd['a'];b=bgd['b'];eps=bgd['eps'];P=bgd['P'];rho=bgd['rho']
I_code=np.trapezoid(4*np.pi*r**4*rho,r)
I_rel=np.trapezoid(4*np.pi*np.exp(-a+3*b)*(eps+P/c**2)*r**4,r)
I_eps=np.trapezoid(4*np.pi*(eps+P/c**2)*r**4,r)
I_epsonly=np.trapezoid(4*np.pi*eps*r**4,r)
A,B,_,_=build_RR_matrices(bgd,1)
print("RR n=1: sigma^2=%.5f ; A/B"%(A[0,0]/B[0,0]), " 4 pi B = %.5e  (I_rel=%.5e)"%(4*np.pi*B[0,0],I_rel))
print("exact V'' = %.5e J (Richardson %.5e);  I_code=%.5e  I_eps+P=%.5e I_rel(e^{-a+3b})=%.5e"%(Vpp,Vpp_R,I_code,I_eps,I_rel))
print("k_exact using I_rel: %.4f s^-2   using I_code: %.4f ;   RR target 27.810"%(-Vpp_R/I_rel,-Vpp_R/I_code))
print("I_rel/I_code = %.5f ; (eps+P)/rho factor: %.5f ; e^{-a+3b} factor: %.5f"%(I_rel/I_code,I_eps/I_code,I_rel/I_eps))
print("also with 4*pi*A: A*4pi*? : -V''/(4 pi A_RR) etc: sigma^2*4piB=",A[0,0]/B[0,0]*4*np.pi*B[0,0], " vs V''=",Vpp_R)
