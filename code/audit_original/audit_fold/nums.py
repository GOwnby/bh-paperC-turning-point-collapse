import sys,json
sys.dont_write_bytecode=True
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import numpy as np, background as bg
from metric import build_metric
bg.mu_e=2.0; bg.Zlat=0.0
rm=23106017651.354734e3
d=build_metric(rm); r=d['r'];rho=d['rho'];P=d['P']
G1=np.array([bg.Gamma1_of_rho(q) if q>1 else 5/3 for q in rho])
I=np.trapezoid(4*np.pi*r**4*rho,r); W=np.trapezoid(4*np.pi*r**2*P,r); G1a=np.trapezoid(G1*P*4*np.pi*r**2,r)/W
print('R [km]',d['R']/1e3,'M',d['M']/bg.M_sun,'I',I,'W_P',W,'G1avg',G1a,'9G1W/I',9*G1a*W/I)
print('2GM/Rc2',2*bg.G*d['M']/(d['R']*bg.c**2),'Pc/(rho c2)',P[0]/(rho[0]*bg.c**2))
yr=3.15576e7; kap=18.733; a=0.004027
c9=9*G1a*W/I
for Md in (1e-9,1e-8,1e-7,1e-6,1e-5):
    be=(4*a/(kap**2*Md/yr))**0.2; al=2/(kap*be**2)
    print('Mdot %g beta %.2f alpha %.3e 3.4117beta %.1f 3.4101beta %.1f  1/(c9 beta)=%.2e'%(Md,be,al,3.41167*be,3.41012*be,1/(c9*be)))
# examples
for Md,nu in ((1e-7,10),(1e-5,100)):
    be=(4*a/(kap**2*Md/yr))**0.2; nt=nu*be
    t=2.33811*be*nt**(2/3); dM=Md/yr*t
    dM2=2.33811*(Md/yr)**(2/3)*(4*a/kap**2)**(1/3)*nu**(2/3)
    print('Mdot',Md,'nu',nu,'beta',be,'nutilde',nt,'t',t,'dM',dM,dM2,'(per yr conv) dM [Msun]')
# ion pressure etc
kB=1.380649e-23;T=1e8;mu=1.66053906660e-27;me=9.1093837015e-31;c=2.99792458e8;hbar=1.054571817e-34
rho0=rm; ne=rho0/(2*mu); pF=hbar*(3*np.pi**2*ne)**(1/3)
EF=np.sqrt((pF*c)**2+(me*c**2)**2)-me*c**2
print('EF MeV',EF/1.602176634e-13,'pF c MeV',pF*c/1.602176634e-13,'(kT/EF)^2',(kB*T/EF)**2,'(kT/(pFc))^2',(kB*T/(pF*c))**2)
Pe=bg.P_of_rho(rho0)
for A in (4,12,16,56): print('A',A,kB*T*rho0/(A*mu)/Pe)
print('equip: sqrt(kT/I)',np.sqrt(kB*T/I),' *beta(1e-7)',np.sqrt(kB*T/I)*107.69,' sqrt(kT*beta^2/I)',np.sqrt(kB*T*107.69**2/I))
# Composition table
cj=json.load(open('/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve/fold_comp.json'))
Md=1e-7/yr; base=None
for k,v in cj.items():
    be=(4*v['a']/(v['kappa']**2*Md))**0.2; al=2/(v['kappa']*be**2); t=3.41167*be; tt=3.41012*be
    if base is None: base=(v['kappa'],v['a'],be,t,al)
    print(k,'rho %.3e Mmax %.4f kappa %.2f a %.3f beta %.1f t(3.41167)=%.1f t(3.41012)=%.1f  dk %.1f%% da %.1f%% dt %.1f%% dalpha %.1f%%'%(v['rho_max']/1e3,v['Mmax'],v['kappa'],v['a']*1e3,be,t,tt,100*(v['kappa']/base[0]-1),100*(v['a']/base[1]-1),100*(t/base[3]-1),100*(al/base[4]-1)))
