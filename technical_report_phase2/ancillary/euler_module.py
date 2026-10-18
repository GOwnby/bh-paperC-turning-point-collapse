"""Shared Eulerian hydro setup, factored out of euler_convergence.py so each
resolution/duration run can import it directly without re-pasting code."""
import numpy as np, time
import eos_neutronized as EN
import eos as E
import initial_profile_neu as IPN
from euler_hydro import EulerGrid, minmod
from scipy.interpolate import interp1d

rho_c = 5e12
_sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = _sol.t, _sol.y[0], np.maximum(_sol.y[1], 1.0)
R_star = r_raw[-1]
Mstar_total = m_raw[-1]
rho_interp = interp1d(r_raw, rho_raw, bounds_error=False, fill_value=1.0)

def make_run(N):
    grid = EulerGrid(1.6*R_star, N, rho_interp, ghost=2)
    dr = grid.dr; g = grid.ghost
    def apply_bc(rho,v,p):
        rho[:g]=rho[g:2*g][::-1]; v[:g]=-v[g:2*g][::-1]; p[:g]=p[g:2*g][::-1]
        rho[-g:]=rho[-g-1]; v[-g:]=v[-g-1]; p[-g:]=p[-g-1]
    def compute_mass(r,rho): return np.concatenate([[0.0], np.cumsum(4*np.pi*r[:-1]**2*rho[:-1]*dr)])
    def gravsrc(r,rho,v):
        m = compute_mass(r,rho)
        comp = np.clip(2*E.G*m/(r*E.c**2),0,0.5)
        return -E.G*m/np.maximum(r,1.0)**2/np.maximum(1-comp,0.5), m
    def flux(rho,v,p,eint):
        return np.array([rho*v, rho*v**2+p, v*(rho*(eint+0.5*v**2)+p)])
    def hll(rhoL,vL,pL,rhoR,vR,pR):
        eiL=grid.eint(rhoL,pL); eiR=grid.eint(rhoR,pR)
        csL=grid.sound_speed(rhoL,pL); csR=grid.sound_speed(rhoR,pR)
        SL=np.minimum(vL-csL,vR-csR); SR=np.maximum(vL+csL,vR+csR)
        UL=np.array([rhoL,rhoL*vL,rhoL*(eiL+0.5*vL**2)]); UR=np.array([rhoR,rhoR*vR,rhoR*(eiR+0.5*vR**2)])
        FL=flux(rhoL,vL,pL,eiL); FR=flux(rhoR,vR,pR,eiR)
        Fh=(SR*FL-SL*FR+SL*SR*(UR-UL))/np.maximum(SR-SL,1e-10)
        return np.where(SL>=0,FL,np.where(SR<=0,FR,Fh)), np.maximum(np.abs(SL),np.abs(SR))
    def recon(q):
        dl=q[1:-1]-q[:-2]; dr_=q[2:]-q[1:-1]; s=minmod(dl,dr_)
        return q[1:-1]-0.5*s, q[1:-1]+0.5*s
    def rhs(rho,v,p):
        apply_bc(rho,v,p); r=grid.r
        rhoL,rhoR=recon(rho); vL,vR=recon(v); pL,pR=recon(p)
        F,smax = hll(rhoR[:-1],vR[:-1],pR[:-1], rhoL[1:],vL[1:],pL[1:])
        dU=np.zeros((3,len(r))); riface=0.5*(r[1:]+r[:-1]); area=riface**2
        idx0=g-1; n_int=F.shape[1]
        for k in range(3):
            fd = area[idx0+1:idx0+n_int]*F[k,1:] - area[idx0:idx0+n_int-1]*F[k,:-1]
            dU[k,idx0+1:idx0+n_int] = -fd/(r[idx0+1:idx0+n_int]**2*dr)
        gg,m = gravsrc(np.abs(r),rho,v)
        dU[1]+=rho*gg+2*p/np.maximum(np.abs(r),1.0); dU[2]+=rho*v*gg
        return dU, np.max(smax)
    def c2p(U):
        D,Sr,Et=U; rho=np.maximum(D,1.0); v=Sr/rho
        eint=np.maximum(Et/rho-0.5*v**2,1e-10)
        Gl=np.where(rho>=EN._rho_tab[0],EN._dlogPdlogrho(np.log(np.maximum(rho,EN._rho_tab[0]))),5/3)
        p=np.maximum(eint*rho*np.maximum(Gl-1,0.1),0.0)
        return rho,v,p
    def total_mass(rho):
        r = np.abs(grid.r[g:-g])
        return np.sum(4*np.pi*r**2*rho[g:-g]*dr)
    return dict(grid=grid, dr=dr, g=g, apply_bc=apply_bc, rhs=rhs, c2p=c2p, total_mass=total_mass)

def run(N, t_final, t_checkpoint=None, state=None, wall_budget=250):
    """Run from t=0 (or from `state` if resuming) to t_final, or until wall_budget seconds elapse.
    Returns (t_reached, rho, v, p, log, finished, wall_time)."""
    ctx = make_run(N)
    grid, dr, g = ctx['grid'], ctx['dr'], ctx['g']
    if state is None:
        rho, v, p = grid.rho.copy(), grid.v.copy(), grid.p.copy()
        t = 0.0
    else:
        rho, v, p, t = state
    rho_c0 = rho[g] if state is None else state[4]
    mass0 = ctx['total_mass'](rho)
    t0 = time.time(); log = []
    step = 0
    while t < t_final:
        dU, smax = ctx['rhs'](rho, v, p)
        dt = min(0.4*dr/max(smax, 1.0), t_final - t)
        U0 = np.array([rho, rho*v, rho*(grid.eint(rho,p)+0.5*v**2)])
        Uh = U0 + 0.5*dt*dU
        rh, vh, ph = ctx['c2p'](Uh); ctx['apply_bc'](rh, vh, ph)
        dU2, _ = ctx['rhs'](rh, vh, ph)
        Un = U0 + dt*dU2
        rho, v, p = ctx['c2p'](Un); ctx['apply_bc'](rho, v, p)
        t += dt; step += 1
        if step % max(1, int(0.002*t_final/dt)) == 0 or step==1:
            m_now = ctx['total_mass'](rho)
            log.append((t, rho[g]/rho_c0, m_now/mass0, np.max(np.abs(v))))
        if np.any(np.isnan(rho)) or np.any(rho < 0):
            return t, rho, v, p, log, 'BLOWUP', time.time()-t0, rho_c0
        if time.time()-t0 > wall_budget:
            return t, rho, v, p, log, 'WALL_LIMIT', time.time()-t0, rho_c0
    return t, rho, v, p, log, 'DONE', time.time()-t0, rho_c0
