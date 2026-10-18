import numpy as np, json, background as bg
from homologous_calibrated import CalibratedHomologousModel
bg.mu_e=2.0
d=json.load(open('fold_step1.json')); kap=18.733; a=0.004027; yr=3.15576e7
m=CalibratedHomologousModel(d['rho_max'])
G1W=m.Gamma1avg*m.W_P
print('I',m.I,'W_P',m.W_P,'Gamma1avg',m.Gamma1avg,'9*G1avg*W_P/I [s^-2] =',9*G1W/m.I)
c9=9*G1W/m.I   # nu = c9 * Delta_r * tau  (omega tau << 1)
for Md in (1e-9,1e-7,1e-5):
    beta=(4*a/(kap**2*Md/yr))**0.2
    print(Md,'beta',beta,'Delta_r*tau needed for nu*beta=1: %.3e s'%(1/(c9*beta)))
# thermal: ion/electron pressure ratio at T
kB=1.380649e-23; T=1e8
rho=d['rho_max']; ne=rho/(2*bg.m_u) if hasattr(bg,'m_u') else rho/(2*1.66053907e-27)
P=bg.P_of_rho(rho); print('P_e(rho_max)',P)
for A in (4,12,16,56):
    ni=rho/(A*1.66053907e-27); print('A',A,'P_ion/P_e at 1e8 K',ni*kB*T/P)
# Fermi energy / kT
hbar=1.054571817e-34; c=2.99792458e8; pF=hbar*(3*np.pi**2*ne)**(1/3); EF=np.sqrt((pF*c)**2+(9.109e-31*c**2)**2)-9.109e-31*c**2
print('EF [MeV]',EF/1.602e-13,'kT/EF',kB*T/EF,'(kT/EF)^2',(kB*T/EF)**2)
