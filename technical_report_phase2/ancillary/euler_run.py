import numpy as np, time, pickle
import eos_neutronized as EN
import eos as E
import initial_profile_neu as IPN
from euler_hydro import EulerGrid, minmod

# Build the same initial profile as the Lagrangian run (rho_c=5e12) for direct comparison
rho_c = 5e12
sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1.0)
R_star = r_raw[-1]
Mstar = m_raw[-1]
print(f"Star: M={Mstar/E.Msun:.4f} Msun, R={R_star/1e3:.1f} km")

from scipy.interpolate import interp1d
rho_interp = interp1d(r_raw, rho_raw, bounds_error=False, fill_value=1.0)

N = 800
rmax = 1.6*R_star
grid = EulerGrid(rmax, N, rho_interp, ghost=2)
dr = grid.dr
g = grid.ghost

def apply_bc(rho, v, p):
    # reflecting at r=0 (v antisymmetric, rho/p symmetric); outflow at outer boundary
    rho[:g] = rho[g:2*g][::-1]; v[:g] = -v[g:2*g][::-1]; p[:g] = p[g:2*g][::-1]
    rho[-g:] = rho[-g-1]; v[-g:] = v[-g-1]; p[-g:] = p[-g-1]

def compute_mass(r, rho):
    # cumulative mass interior to each cell (cell-centered r, uniform dr)
    return np.concatenate([[0.0], np.cumsum(4*np.pi*r[:-1]**2*rho[:-1]*dr)])

def gravitational_source(r, rho, v):
    m = compute_mass(r, rho)
    # leading relativistic correction to the Newtonian potential gradient (weak-field):
    # g_eff = -Gm/r^2 * (1 + small correction from pressure/compactness); here we include
    # the standard leading term (1 - 2Gm/(rc^2))^{-1} factor from the TOV structure equation
    compactness = np.clip(2*E.G*m/(r*E.c**2), 0, 0.5)
    g_grav = -E.G*m/np.maximum(r,1.0)**2 / np.maximum(1-compactness, 0.5)
    return g_grav, m

def flux(rho, v, p, eint):
    F0 = rho*v
    F1 = rho*v**2 + p
    F2 = v*(rho*(eint+0.5*v**2)+p)
    return np.array([F0,F1,F2])

def hll_flux(rhoL,vL,pL, rhoR,vR,pR, grid):
    eintL = grid.eint(rhoL,pL); eintR = grid.eint(rhoR,pR)
    csL = grid.sound_speed(rhoL,pL); csR = grid.sound_speed(rhoR,pR)
    SL = np.minimum(vL-csL, vR-csR)
    SR = np.maximum(vL+csL, vR+csR)
    UL = np.array([rhoL, rhoL*vL, rhoL*(eintL+0.5*vL**2)])
    UR = np.array([rhoR, rhoR*vR, rhoR*(eintR+0.5*vR**2)])
    FL = flux(rhoL,vL,pL,eintL); FR = flux(rhoR,vR,pR,eintR)
    F_hll = (SR*FL - SL*FR + SL*SR*(UR-UL))/np.maximum(SR-SL,1e-10)
    F = np.where(SL>=0, FL, np.where(SR<=0, FR, F_hll))
    return F, np.maximum(np.abs(SL),np.abs(SR))

def reconstruct(q):
    dq_l = q[1:-1]-q[:-2]; dq_r = q[2:]-q[1:-1]
    slope = minmod(dq_l, dq_r)
    qL = q[1:-1] - 0.5*slope   # left state at interface i-1/2 -> actually build interface values
    qR = q[1:-1] + 0.5*slope
    return qL, qR   # qL,qR at cell centers 1..N-2 (reconstructed edge values within each cell)

def rhs(rho, v, p):
    apply_bc(rho, v, p)
    r = grid.r
    rhoL, rhoR = reconstruct(rho)
    vL, vR = reconstruct(v)
    pL, pR = reconstruct(p)
    # interface i+1/2 between cell i and i+1: use R-state of cell i, L-state of cell i+1
    F, smax = hll_flux(rhoR[:-1], vR[:-1], pR[:-1], rhoL[1:], vL[1:], pL[1:], grid)
    n_int = F.shape[1]
    dU = np.zeros((3, len(r)))
    r_iface = 0.5*(r[1:]+r[:-1])
    area = r_iface**2
    # flux divergence (finite volume, spherical): use interior interfaces only (aligned with reconstructed cells 1..N-2)
    idx0 = g-1
    for k in range(3):
        flux_diff = area[idx0+1:idx0+1+n_int-1]*F[k,1:] - area[idx0:idx0+n_int-1]*F[k,:-1]
        dU[k, idx0+1:idx0+n_int] = -flux_diff/(r[idx0+1:idx0+n_int]**2*dr)
    g_grav, m = gravitational_source(np.abs(r), rho, v)
    eint = grid.eint(rho,p)
    dU[1] += rho*g_grav + 2*p/np.maximum(np.abs(r),1.0)
    dU[2] += rho*v*g_grav
    return dU, np.max(smax)

def cons_to_prim(U, rho_old, p_old):
    D,Sr,Etot = U
    rho = np.maximum(D, 1.0)
    v = Sr/rho
    eint = np.maximum(Etot/rho - 0.5*v**2, 1e-10)
    Gamma_loc = np.where(rho>=EN._rho_tab[0], EN._dlogPdlogrho(np.log(np.maximum(rho,EN._rho_tab[0]))), 5/3)
    p = np.maximum(eint*rho*np.maximum(Gamma_loc-1,0.1), 0.0)
    return rho, v, p

rho, v, p = grid.rho.copy(), grid.v.copy(), grid.p.copy()
dU, smax = rhs(rho, v, p)
print("initial RHS test: max|dU|=", np.max(np.abs(dU)), " max wave speed=", smax, " dr=", dr)

print("\n=== Running the Eulerian evolution ===")
rho, v, p = grid.rho.copy(), grid.v.copy(), grid.p.copy()
t = 0.0
t_final = 0.10   # seconds -- comparable to where interesting Lagrangian dynamics occurred (t~9-15 in L/c~5.3ms units => ~0.05-0.08s)
CFL = 0.4
t0 = time.time()
step = 0
rho_c_init = rho[g]
log = []
while t < t_final:
    dU, smax = rhs(rho, v, p)
    dt = CFL*dr/max(smax, 1.0)
    dt = min(dt, t_final-t)
    U0 = grid.cons() if False else np.array([rho, rho*v, rho*(grid.eint(rho,p)+0.5*v**2)])
    # RK2 (midpoint)
    U_half = U0 + 0.5*dt*dU
    rho_h, v_h, p_h = cons_to_prim(U_half, rho, p)
    apply_bc(rho_h, v_h, p_h)
    dU2, smax2 = rhs(rho_h, v_h, p_h)
    U_new = U0 + dt*dU2
    rho, v, p = cons_to_prim(U_new, rho, p)
    apply_bc(rho, v, p)
    t += dt; step += 1
    if step % 50 == 0:
        r_phys = np.abs(grid.r)
        # find current "surface" (density drops below 1% of initial surface density proxy)
        surf_idx = np.where(rho[g:-g] > 1e6)[0]
        R_surf_now = r_phys[g+surf_idx[-1]] if len(surf_idx)>0 else np.nan
        log.append((t, rho[g]/rho_c_init, R_surf_now/R_star, np.max(np.abs(v))))
        if step % 500 == 0:
            print(f"step={step:5d} t={t*1000:.3f}ms  rho_c/rho0={rho[g]/rho_c_init:.5f}  R_surf/R0={R_surf_now/R_star:.4f}  max|v|={np.max(np.abs(v)):.3e}  wall={time.time()-t0:.1f}s", flush=True)
    if np.any(np.isnan(rho)) or np.any(np.isnan(v)) or np.any(rho<0):
        print(f"BLOWUP/negative density at step {step}, t={t*1000:.3f}ms"); break
    if time.time()-t0 > 240:
        print(f"wall-time budget reached at step {step}, t={t*1000:.2f}ms"); break

pickle.dump({'log':log, 't_final':t, 'rho_final':rho, 'v_final':v, 'p_final':p, 'r':grid.r, 'R_star':R_star, 'rho_c_init':rho_c_init},
            open('euler_run.pkl','wb'))
print(f"\ndone: reached t={t*1000:.3f} ms in {step} steps, wall={time.time()-t0:.1f}s")
