import sys; sys.path.insert(0,'..'); sys.path.insert(0,'../linear_perturbation')
"""
Same Eulerian setup as euler_module.py, but the initial density is perturbed
by the validated fundamental eigenfunction at a controlled, explicit amplitude,
instead of relying on grid-truncation noise to accidentally seed whatever mode
a given resolution happens to excite.
"""
import numpy as np, time, pickle
import eos_neutronized as EN
import eos as E
from euler_hydro import EulerGrid, minmod
from scipy.interpolate import interp1d

d = pickle.load(open('eigenfunction_data.pkl','rb'))
r_eig, xi_eig = d['r'], d['xi']
R_star = d['R_star']; rho_c = d['rho_c']
r_raw, m_raw, rho_raw = d['r_raw'], d['m_raw'], d['rho_raw']

xi_of_r = interp1d(r_eig, xi_eig, bounds_error=False, fill_value=(0.0, xi_eig[-1]))
rho0_of_r = interp1d(r_raw, rho_raw, bounds_error=False, fill_value=1.0)

def perturbed_rho(r, epsilon):
    r = np.abs(r)
    rho0 = rho0_of_r(r)
    # Near r=0, xi ~ a1*r exactly (regularity), so xi/r -> a1 (finite) and
    # xi'+2xi/r -> 3*a1 (finite, the homologous-center limit). Interpolation
    # noise in xi_of_r near r=0 can corrupt the raw xi/r ratio there, so use
    # the analytic near-center limit directly below a small cutoff radius,
    # rather than trust 2*xi/r computed from interpolated values at tiny r.
    r_cut = R_star*1e-3
    h = R_star*1e-4
    xip = (xi_of_r(r+h) - xi_of_r(np.maximum(r-h,0)))/(h + np.minimum(r,h))
    a1 = xi_of_r(r_cut)/r_cut   # near-center slope, evaluated safely away from r=0
    divxi_far = xip + 2*xi_of_r(r)/np.maximum(r, r_cut)
    divxi = np.where(r > r_cut, divxi_far, 3*a1)
    delta_rho_over_rho0 = -divxi * R_star
    delta_rho_over_rho0 = np.clip(delta_rho_over_rho0, -50, 50)   # guard against any residual interpolation spike
    return np.maximum(rho0*(1 + epsilon*delta_rho_over_rho0), 1.0)

def make_run(N, epsilon):
    grid = EulerGrid(1.6*R_star, N, lambda r: perturbed_rho(r, epsilon), ghost=2)
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
        r = np.abs(grid.r[g:-g]); return np.sum(4*np.pi*r**2*rho[g:-g]*dr)
    return dict(grid=grid, dr=dr, g=g, apply_bc=apply_bc, rhs=rhs, c2p=c2p, total_mass=total_mass)

def run(N, epsilon, t_final, wall_budget=250):
    ctx = make_run(N, epsilon)
    grid, dr, g = ctx['grid'], ctx['dr'], ctx['g']
    rho, v, p = grid.rho.copy(), grid.v.copy(), grid.p.copy()
    rho_c0_true = rho_raw[0]   # UNPERTURBED central density (r_raw[0]~1m, the true center), for a consistent ratio baseline
    t = 0.0; t0 = time.time(); log = []
    while t < t_final:
        dU, smax = ctx['rhs'](rho, v, p)
        dt = min(0.4*dr/max(smax, 1.0), t_final - t)
        U0 = np.array([rho, rho*v, rho*(grid.eint(rho,p)+0.5*v**2)])
        Uh = U0 + 0.5*dt*dU
        rh, vh, ph = ctx['c2p'](Uh); ctx['apply_bc'](rh, vh, ph)
        dU2, _ = ctx['rhs'](rh, vh, ph)
        Un = U0 + dt*dU2
        rho, v, p = ctx['c2p'](Un); ctx['apply_bc'](rho, v, p)
        t += dt
        log.append((t, rho[g]/rho_c0_true, v[g]))
        if np.any(np.isnan(rho)) or np.any(rho < 0):
            return t, rho, v, p, log, 'BLOWUP', time.time()-t0
        if time.time()-t0 > wall_budget:
            return t, rho, v, p, log, 'WALL_LIMIT', time.time()-t0
    return t, rho, v, p, log, 'DONE', time.time()-t0

print("euler_module_seeded loaded OK")
