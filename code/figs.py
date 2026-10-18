import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices, make_basis
from tov_solve import rho_of_P
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
out='/mnt/user-data/outputs/consolidated/paperC_core_instability/'
# ---- Newtonian structure + Newtonian RR (c->inf limit of same variational form) ----
def newtonian_s2(rc,nb=5,n=12000):
    Pc=bg.P_of_rho(rc)
    def rhs(r,y):
        P,M=y
        if P<=1e-10*Pc: return [0,0]
        rho=rho_of_P(P,rc); return [-G*rho*M/r**2, 4*np.pi*r**2*rho]
    ev=lambda r,y:y[0]-1e-10*Pc; ev.terminal=True
    r0=1.0; s=solve_ivp(rhs,[r0,5e7],[Pc,4/3*np.pi*r0**3*rc],events=ev,rtol=1e-10,atol=[Pc*1e-14,1e-3],dense_output=True,max_step=2e4)
    R=s.t_events[0][0]; r=np.linspace(r0,R,n); P=np.clip(s.sol(r)[0],0,None); M=s.sol(r)[1]
    rho=np.array([rho_of_P(p,rc) for p in P]); dPdr=-G*rho*M/r**2
    Gam=np.array([bg.Gamma1_of_rho(x) if x>1 else 5/3 for x in rho])
    fns=make_basis(R,nb); xis=[f(r) for f in fns]; zs=[r**2*x for x in xis]; dz=[np.gradient(z,r) for z in zs]
    Pi=Gam*P/r**2; Q=4*dPdr*r; S=rho*r**2
    A=np.array([[np.trapezoid(Pi*dz[i]*dz[j]+Q*xis[i]*xis[j],r) for j in range(nb)] for i in range(nb)])
    B=np.array([[np.trapezoid(S*xis[i]*xis[j],r) for j in range(nb)] for i in range(nb)])
    return eigh(A,B)[0][0], 2*G*M[-1]/(R*c**2)
rcs=[1e8,3e8,1e9,3e9,1e10,3e10,1e11,3e11,1e12,3e12,1e13]
comp=[];rat=[]
for rc in rcs:
    sN,C=newtonian_s2(rc)
    d=build_metric(rc); A,B,_,_=build_RR_matrices(d,5); sG=eigh(A,B)[0][0]
    comp.append(2*G*d['M']/(d['R']*c**2)); rat.append(sG/sN); print(rc,comp[-1],sN,sG,sG/sN,flush=True)
# ---- mass curve ----
rr=np.logspace(12.6,14.6,50); Ms=np.array([build_metric(x)['M']/bg.M_sun for x in rr])
rpk=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(13,14),method='bounded',options={'xatol':1e-4}); rpk=10**rpk.x
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].semilogx(comp,rat,'o-',color='#1f5fb4'); ax[0].axhline(1,ls='--',c='gray'); ax[0].axhline(0,ls=':',c='gray')
ax[0].set_xlabel(r'compactness $2GM/(Rc^2)$'); ax[0].set_ylabel(r'$\sigma^2_{\rm GR}/\sigma^2_{\rm Newtonian}$'); ax[0].set_title('Convergence to the Newtonian limit'); ax[0].grid(alpha=.3)
ax[1].semilogx(rr[rr<=rpk],Ms[rr<=rpk],c='#1f5fb4',label='stable'); ax[1].semilogx(rr[rr>=rpk],Ms[rr>=rpk],c='#b02a18',label='unstable')
ax[1].axhline(1.4176,ls='--',c='gray'); ax[1].text(1.2e13,1.4180,'Chandrasekhar-Tooper (1964): 1.4176',fontsize=8,color='gray')
ax[1].axvline(rpk,ls=':',c='k'); ax[1].plot([5e13],[build_metric(5e13)['M']/bg.M_sun],'*',ms=14,c='#7b1fa2',label='test configuration')
ax[1].set_xlabel(r'central density $\rho_c$ [kg/m$^3$]'); ax[1].set_ylabel(r'mass $M$ [$M_\odot$]'); ax[1].set_title('TOV + Chandrasekhar EOS')
Mpk=build_metric(rpk)['M']/bg.M_sun; ax[1].text(0.97,0.06,f'$M_{{max}}$={Mpk:.4f} $M_\\odot$ at {rpk:.3e}',transform=ax[1].transAxes,ha='right',fontsize=8); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(out+'fig5_gr_validation.png',dpi=200); plt.close()
print("peak",rpk,Ms.max())
# ---- RR convergence ----
d=build_metric(5e13); vals=[eigh(*build_RR_matrices(d,n)[:2])[0][0] for n in range(1,6)]
plt.figure(figsize=(5,3.6)); plt.plot(range(1,6),vals,'o-'); plt.xlabel('number of basis functions'); plt.ylabel(r'$\sigma^2$ [s$^{-2}$]'); plt.grid(alpha=.3); plt.tight_layout(); plt.axhline(-28.556,ls=':',c='gray'); plt.savefig(out+'fig_rayleigh_ritz_convergence.png',dpi=200); plt.close(); print(vals)
