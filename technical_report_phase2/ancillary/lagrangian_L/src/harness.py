"""Common driver: build a star from one of the initial-data families, evolve with fixed dt, record diagnostics."""
import time, json, sys
import numpy as np
from initdata import Setup
from lag import Lag, relax_equilibrium, linear_modes, fd_jacobian

C_SI = 2.998e8


def make_star(kind, data, N, newt=False, rho_c=5e12, relax=None, **kw):
    """data in {'tov','newton'}; relax: None | 'db' (discrete-equilibrium Newton relaxation, labels reset)"""
    S = Setup(kind, rho_c=rho_c, newt=newt)
    if data == 'tov':
        d = S.tov_data(N)
    elif data == 'newton':
        d = S.newton_data(N, newt=newt)
    else:
        raise ValueError(data)
    star = S.make(d, **kw)
    if relax == 'db':
        relax_equilibrium(star)
        star.r_lab = star.R.copy()
    star.S = S
    return star


def mass_radii(star, geo, fr=(0.5, 0.9, 0.99, 0.999)):
    Mb = np.concatenate([[0.0], np.cumsum(star.dM)])
    f = Mb / Mb[-1]
    return [float(np.interp(x, f, star.R)) for x in fr]


def snapshot(star, t0state):
    geo = star.geometry(star.R, star.U); star.accel(geo)
    R = star.R
    dRl = np.diff(R) / np.diff(star.r_lab)
    Mb = np.concatenate([[0.0], np.cumsum(star.dM)]); f = Mb / Mb[-1]
    k = int(np.argmin(dRl))
    fout = 1.0 - f[k + 1]
    interior = f[1:] < 0.999            # zones whose outer interface encloses < 99.9 % of the mass
    return dict(t=star.t, rho_c=float(geo['e'][0]), Rpmin=float(dRl[k]), kmin=k, mass_out_at_min=float(fout),
                Rpmin_int=float(np.min(dRl[interior])), Rs=float(R[-1]), M=float(geo['m'][-1]),
                rl=mass_radii(star, geo), umax=float(np.max(np.abs(star.U))), emax=float(geo['e'].max()),
                dR_L2=float(np.sqrt(np.sum(star.dM * (R[1:] - t0state['R'][1:])**2) / star.dM.sum())),
                dR_bulk=float(np.max(np.abs(R[1:] - t0state['R'][1:])[f[1:] < 0.99])),
                q_max=float(geo['q'].max()))


def evolve(star, T, nout=10, C=0.4, dt_scale=1.0, stop_cross=True, verbose=False, max_wall=1e9):
    """fixed-dt leapfrog to T with nout equally spaced outputs; stops at ordering loss"""
    star.start()
    st0 = dict(R=star.R.copy())
    dt0 = star.dt_cfl(star.geo, C) * dt_scale
    nsteps = nout * int(np.ceil(T / dt0 / nout))
    dt = T / nsteps
    snaps = [snapshot(star, st0)]
    snaps[0]['rho_c0'] = snaps[0]['rho_c']
    status = 'tmax'; t_cross = None; t_wall = time.time()
    gprev = np.diff(star.R)
    minRp_hist = []
    for n in range(1, nsteps + 1):
        star.step(dt)
        gaps = np.diff(star.R)
        gmin = gaps.min()
        if gmin <= 0:
            kk = int(np.argmin(gaps))
            # linear interpolation of the time at which the gap of the crossing zone vanished
            # (use the gap of that zone at the previous step ~ stored in gap_prev_arr)
            t_cross = star.t - dt + dt * gprev[kk] / (gprev[kk] - gaps[kk])
            status = 'crossing'
            s = snapshot(star, st0); snaps.append(s)
            break
        gprev = gaps
        if n % (nsteps // nout) == 0:
            snaps.append(snapshot(star, st0))
        if time.time() - t_wall > max_wall:
            status = 'wall'; break
    return dict(snaps=snaps, status=status, t_cross=t_cross, nsteps=n, dt=dt, wall=time.time() - t_wall)


def evolve_adaptive(star, T, nout=10, C=0.4, dt_scale=1.0, stop_cross=True, max_wall=1e9, rho_stop=300.0):
    """as evolve(), but (i) the run stops with status 'collapse' when rho_c exceeds rho_stop x its initial value (the EOS tables end at 1e16 kg/m^3), (ii) the step is re-limited every step by the CFL bound of the CURRENT state, dt = C*dt_scale*dt_cfl(geo), cut to hit the output times exactly
    (needed when the sound speed grows by large factors during collapse)"""
    star.start()
    st0 = dict(R=star.R.copy())
    snaps = [snapshot(star, st0)]
    snaps[0]['rho_c0'] = snaps[0]['rho_c']
    status = 'tmax'; t_cross = None; t_wall = time.time(); nsteps = 0; e0 = snaps[0]['rho_c']
    gprev = np.diff(star.R)
    touts = np.linspace(0.0, T, nout + 1)[1:]
    done = False
    for tout in touts:
        while star.t < tout - 1e-12:
            dt = min(C * dt_scale * star.dt_cfl(star.geo, 1.0), tout - star.t)
            star.step(dt); nsteps += 1
            gaps = np.diff(star.R)
            if gaps.min() <= 0:
                kk = int(np.argmin(gaps))
                t_cross = star.t - dt + dt * gprev[kk] / (gprev[kk] - gaps[kk])
                status = 'crossing'; snaps.append(snapshot(star, st0)); done = True; break
            gprev = gaps
            if star.geo['e'][0] > rho_stop * e0:
                status = 'collapse'; snaps.append(snapshot(star, st0)); done = True; break
        if done: break
        snaps.append(snapshot(star, st0))
        if time.time() - t_wall > max_wall: status = 'wall'; break
    return dict(snaps=snaps, status=status, t_cross=t_cross, nsteps=nsteps, dt=T / max(nsteps, 1), wall=time.time() - t_wall)


def rel(snaps, key='rho_c'):
    r0 = snaps[0]['rho_c']
    return [s['rho_c'] / r0 - 1 for s in snaps]


def s2(x):  # seconds^-2 factor
    return None
