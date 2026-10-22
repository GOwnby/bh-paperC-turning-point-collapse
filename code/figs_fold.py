import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from fold_dyn import *
sc=1.5*kap/(-coef[0])
def traj(Mdot,S=30.0):
    md=Mdot/yr; beta=(4*a/(kap**2*md))**0.2; cc=-coef[0]*sc; kk=md/(9*a)
    fz=lambda t: cc*(-md*t)/(9*a); t0=-S*beta; e0=-np.sqrt(-kk*t0); ed=0.5*kk/np.sqrt(-kk*t0)
    rhs=lambda t,y:[y[1],sc*F(y[0])+fz(t)]
    ev=lambda t,y:(y[0]**-3-1)-0.3; ev.terminal=True
    s=solve_ivp(rhs,(t0,t0+60*beta),[1-e0,-ed],method='DOP853',events=ev,rtol=1e-11,atol=1e-14,max_step=beta/40)
    return s,beta,md
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
f=lambda s,u:[u[1],u[0]**2+s]; ev=lambda s,u:u[0]-1e4; ev.terminal=True
S=40; pi=solve_ivp(f,(-S,5),[-np.sqrt(S)-1/(8*S**2),0.5/np.sqrt(S)],method='DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.01)
for Md,col in ((1e-9,'C0'),(1e-7,'C1'),(1e-5,'C2')):
    s,beta,md=traj(Md); y=s.y[0]**-3-1; al=2/(kap*beta**2)
    ax[0].plot(s.t/beta,y/al,color=col,label=f'$\\dot M$={Md:.0e} $M_\\odot$/yr ($\\alpha$={al:.1e})')
ax[0].plot(pi.t,pi.y[0],'k--',lw=1,label='Painlevé-I (universal)')
ss=np.linspace(-12,0,200); ax[0].plot(ss,-np.sqrt(-ss),'gray',lw=0.8,label='quasi-static $-\\sqrt{-s}$')
ax[0].set_xlim(-12,3.6); ax[0].set_ylim(-5,14); ax[0].set_xlabel('$s=t/\\beta$  ($s=0$: $M=M_{\\max}$)'); ax[0].set_ylabel('$u=y/\\alpha$'); ax[0].legend(fontsize=7,loc='upper left'); ax[0].set_title('Fold passage in scaled variables')
Mds=np.logspace(-9,-5,9); tp=[]; 
for Md in Mds:
    md=Md/yr; beta=(4*a/(kap**2*md))**0.2; tp.append(3.41167*beta)
ax[1].loglog(Mds,tp,'k-',label='$3.412\\,\\beta(\\dot M)\\propto\\dot M^{-1/5}$')
d=json.load(open('fold_dyn2.json'))
pts=[(r['Mdot'],r['t_half']) for r in d['out'] if abs(r['scale']-sc)<1e-6 and r['S']==45]
ax[1].loglog(*zip(*pts),'ro',label='homologous model (rescaled), $\\lambda=0.5$')
ax[1].axhline(0.187,color='b',ls=':',label='$1/\\sigma$ at test config (0.187 s)')
ax[1].set_xlabel('accretion rate $\\dot M$ [$M_\\odot$/yr]'); ax[1].set_ylabel('time from $M=M_{\\max}$ to collapse [s]'); ax[1].legend(fontsize=7); ax[1].set_title('Collapse delay vs accretion rate')
plt.tight_layout(); import os; plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','fig_fold_passage.png'),dpi=160)
