"""Conservative staggered-grid Lagrangian scheme for spherically symmetric Misner-Sharp collapse of a barotropic fluid.

CONTINUUM SYSTEM (G = c = 1; derived/verified in verify_ms_equations.py; Misner & Sharp 1964, May & White 1966):
   metric  ds^2 = -a^2 dt^2 + X^2 dr^2 + R^2 dOmega^2,  r a fixed Lagrangian (comoving) label, fluid at rest in r.
   U = (1/a) dR/dt  (proper-time derivative),  Gamma = R'/X,   m defined by  Gamma^2 = 1 + U^2 - 2m/R
   (1)  dR/dt = a U
   (2)  m' = 4 pi e R^2 R'            (Hamiltonian constraint)
   (3)  dm/dt = -4 pi a p R^2 U       (energy equation)
   (4)  a'/a = -p'/(e+p)              (Euler, radial component;  for a barotropic fluid:  a h = const)
   (5)  dU/dt = -a [ (m + 4 pi R^3 p)/R^2 + Gamma^2 p' / ((e+p) R') ]      (Euler, time component / momentum)
   (6)  dGamma/dt = (a'/(a X)) dR/dt  (so Gamma is NOT constant along a shell unless a' = 0)
   closure p = P(e)  (barotropic, report's EOS).  Rest-mass/number conservation: dM_b = 4 pi R^2 rho_N dR / Gamma.

SCHEME (second order):  interfaces j = 0..N carry R_j, U_j (R_0 = U_0 = 0 regular centre), zones k = j+1/2 carry
   the fixed number content dM_k, rho_N,k = Gamma_k dM_k / V_k  (V_k = 4 pi/3 (R_{k+1}^3 - R_k^3)),  e_k = e(rho_N,k),
   p_k = p(e_k).  m_j = sum_{k<j} e_k V_k  (discrete constraint (2), exact),  Gamma_j^2 = 1 + U_j^2 - 2 m_j/R_j,
   Gamma_k = (Gamma_j+Gamma_{j+1})/2 (the coupled rho_N-m-Gamma relation is solved by fixed-point iteration to max|dGamma| < 1e-13; a fixed count of 3 iterations leaves an O(dt) error in
   fast, strongly relativistic motion, see test_i_dust.py).
   Interface acceleration (momentum eq. (5) with  Gamma^2 p'/R' = Gamma (4 pi R^2 rho_N/Gamma...) written in mass form):
        dU_j/dt = -a_j [ (m_j + 4 pi R_j^3 pbar_j)/R_j^2 + Gamma_j 4 pi R_j^2 (1/h)_j (pq_{j+1/2} - pq_{j-1/2}) / dM^int_j ]
   a_j = 1/h_j (a h = const, normalised to 1 at the vacuum surface h=1), h = (e+p)/rho_N.
   pq = p + q, q = von Neumann-Richtmyer artificial viscosity (optional, off by default).
   Free surface: p = 0 outside the last zone.  Time stepping: kick-drift-kick leapfrog with a-trapezoid in the drift.
   CFL: dt = C min_k a_k (R_{k+1}-R_k)/(Gamma_k c_s,k),  c_s^2 = Gamma1 p/(e+p).
Newtonian mode (newt=True, c -> infinity):  Gamma = a = h = 1, e = rho_N = rho, force = -m/R^2 - 4 pi R^2 dp/dM.
"""
import numpy as np
from scipy.linalg import eig, solve

FOUR_PI = 4.0 * np.pi


class Lag:
    def __init__(self, eos, dM, R0, U0=None, newt=False, cq2=0.0, cl=0.0, p_floor=0.0, p_ext=0.0, n_iter=3, label=None, gtol=1e-13, n_iter_max=200):
        self.eos = eos
        self.newt = newt
        self.dM = np.asarray(dM, float).copy()
        self.N = len(self.dM)
        self.R = np.asarray(R0, float).copy()           # length N+1, R[0] = 0
        self.U = np.zeros(self.N + 1) if U0 is None else np.asarray(U0, float).copy()
        self.R[0] = 0.0; self.U[0] = 0.0
        self.r_lab = self.R.copy()                       # Lagrangian label r = initial areal radius
        self.cq2, self.cl = cq2, cl
        self.p_floor = p_floor
        self.p_ext = p_ext
        self.n_iter = n_iter              # number of fixed-point iterations when gtol = 0 (used only by test_e)
        self.gtol = gtol                  # gtol > 0: iterate the rho_N-m-Gamma coupling until max|dGamma| < gtol (cap n_iter_max)
        self.n_iter_max = n_iter_max
        self.dMpad = np.concatenate([self.dM, [0.0]])
        self.dMint = 0.5 * (self.dMpad[:-1] + self.dMpad[1:])
        self.t = 0.0
        self.Gint = np.ones(self.N + 1)
        self.geo = None
        self.nevals = 0

    # ------------------------------------------------------------------------------------------------------
    def geometry(self, R, U, Gguess=None):
        eos, N = self.eos, self.N
        R3 = R**3
        V = (FOUR_PI / 3.0) * (R3[1:] - R3[:-1])
        G = (self.Gint if Gguess is None else Gguess).copy()
        Rj = R[1:]
        if self.newt:
            rN = self.dM / V
            e, p, g1 = eos.thermo(rN)
            m = np.concatenate([[0.0], np.cumsum(self.dM)])
            G[:] = 1.0
        else:
            for _ in range(self.n_iter_max if self.gtol > 0 else self.n_iter):
                Gold = G.copy()
                Gz = 0.5 * (G[:-1] + G[1:])
                rN = self.dM * Gz / V
                e, p, g1 = eos.thermo(rN)
                m = np.concatenate([[0.0], np.cumsum(e * V)])
                G[0] = 1.0
                G[1:] = np.sqrt(np.maximum(1.0 + U[1:]**2 - 2.0 * m[1:] / Rj, 1e-12))
                if self.gtol > 0 and np.max(np.abs(G - Gold)) < self.gtol:
                    break
        if self.p_floor > 0:
            p = np.maximum(p, self.p_floor)
        self.nevals += 1
        return dict(R=R, U=U, V=V, G=G, rN=rN, e=e, p=p, g1=g1, m=m)

    def accel(self, geo):
        """dU/dt (length N, interfaces 1..N) and a_j (length N)"""
        R, U, G, rN, e, p, m = geo['R'], geo['U'], geo['G'], geo['rN'], geo['e'], geo['p'], geo['m']
        N = self.N
        Rj, mj, Gj = R[1:], m[1:], G[1:]
        # artificial viscosity (compression only), zone-centred
        q = np.zeros(N)
        if self.cq2 > 0 or self.cl > 0:
            dU = U[1:] - U[:-1]
            comp = np.minimum(dU, 0.0)
            cs = np.sqrt(np.maximum(geo['g1'] * p / np.maximum(e + p, 1e-300), 0.0))
            q = self.cq2 * rN * comp**2 + self.cl * rN * cs * (-comp)
        pq = p + q
        ppad = np.concatenate([p, [self.p_ext]])
        pqpad = np.concatenate([pq, [self.p_ext]])
        pbar = 0.5 * (ppad[:-1] + ppad[1:])
        dp = pqpad[1:] - pqpad[:-1]
        if self.newt:
            F = -(mj / Rj**2 + FOUR_PI * Rj**2 * dp / self.dMint)
            a = np.ones(N)
        else:
            h = (e + p) / rN
            hpad = np.concatenate([h, [1.0]])
            ihpad = 1.0 / hpad
            hj = 0.5 * (hpad[:-1] + hpad[1:])
            ihj = 0.5 * (ihpad[:-1] + ihpad[1:])
            a = 1.0 / hj
            F = -a * ((mj + FOUR_PI * Rj**3 * pbar) / Rj**2 + Gj * FOUR_PI * Rj**2 * ihj * dp / self.dMint)
        geo['q'] = q
        geo['a'] = a
        return F, a

    def force_R(self, Rin):
        """acceleration for static (U = 0) configuration with radii Rin[1:]; used for equilibrium solves"""
        R = np.concatenate([[0.0], Rin])
        geo = self.geometry(R, np.zeros(self.N + 1), Gguess=np.ones(self.N + 1))
        F, a = self.accel(geo)
        return F, a

    def dt_cfl(self, geo, C):
        R, G, e, p, g1 = geo['R'], geo['G'], geo['e'], geo['p'], geo['g1']
        if self.newt:
            cs = np.sqrt(g1 * p / e)
            a = 1.0; Gz = 1.0
        else:
            cs = np.sqrt(g1 * p / (e + p))
            a = rN_h = geo['rN'] / (e + p)      # 1/h
            Gz = 0.5 * (G[:-1] + G[1:])
        dr = R[1:] - R[:-1]
        return C * np.min(a * dr / (Gz * np.maximum(cs, 1e-30)))

    # ------------------------------------------------------------------------------------------------------
    def start(self):
        self.geo = self.geometry(self.R, self.U, Gguess=self._G0())
        self.Gint = self.geo['G'].copy()
        for _ in range(3):   # converge initial Gamma (depends on U=0, m)
            self.geo = self.geometry(self.R, self.U)
            self.Gint = self.geo['G'].copy()
        self.F, self.a = self.accel(self.geo)

    def _G0(self):
        return np.ones(self.N + 1)

    def step(self, dt):
        N = self.N
        R, U = self.R, self.U
        Uh = U.copy()
        Uh[1:] += 0.5 * dt * self.F
        # drift with a-trapezoid (predictor-corrector in a)
        Rp = R.copy(); Rp[1:] += dt * self.a * Uh[1:]
        geop = self.geometry(Rp, Uh)
        _, ap = self.accel(geop)
        Rn = R.copy(); Rn[1:] += dt * 0.5 * (self.a + ap) * Uh[1:]
        Upred = Uh.copy(); Upred[1:] += 0.5 * dt * self.F
        geon = self.geometry(Rn, Upred, Gguess=geop['G'])
        Fn, an = self.accel(geon)
        Un = Uh.copy(); Un[1:] += 0.5 * dt * Fn
        # refresh Gamma with the final U (cheap: only changes through U^2)
        self.R, self.U, self.t = Rn, Un, self.t + dt
        self.geo, self.F, self.a = geon, Fn, an
        self.Gint = geon['G'].copy()
        self.geo['U'] = Un

    # ------------------------------------------------------------------------------------------------------
    def diagnostics(self, fresh=True):
        if fresh:
            geo = self.geometry(self.R, self.U)
            self.accel(geo)
        else:
            geo = self.geo
        R, U = self.R, self.U
        dRl = np.diff(R) / np.diff(self.r_lab)       # R' = dR/dr with label r = initial radius
        k = int(np.argmin(dRl))
        d = dict(t=self.t, rho_c=float(geo['e'][0]), Rp_min=float(dRl[k]), k_min=k, R_surf=float(R[-1]),
                 M=float(geo['m'][-1]), umax=float(np.max(np.abs(U))), e_max=float(geo['e'].max()),
                 Rp_bulk=float(np.min(dRl[: int(0.9 * self.N)])))
        return d

    def energy_newt(self, geo=None):
        geo = geo or self.geometry(self.R, self.U)
        R, U = self.R, self.U
        KE = 0.5 * np.sum(self.dMint * U[1:]**2)
        # internal energy per unit mass u(rho) = int p/rho^2 drho : tabulate by quadrature on the fly (cache)
        if not hasattr(self, '_u_tab'):
            rN = np.geomspace(self.eos.rN_tab[0], self.eos.rN_tab[-1], 4000)
            _, p, _ = self.eos.thermo(rN)
            integ = p / rN**2
            lr = np.log(rN)
            u = np.concatenate([[0], np.cumsum(0.5 * (integ[1:] * rN[1:] + integ[:-1] * rN[:-1]) * np.diff(lr))])
            u0 = 1.5 * p[0] / rN[0]
            self._u_tab = (lr, u + u0)
        lr, u = self._u_tab
        uk = np.interp(np.log(geo['rN']), lr, u)
        IE = np.sum(self.dM * uk)
        W = -np.sum(geo['m'][1:] * self.dMint / R[1:])
        return KE, IE, W, KE + IE + W


# ----------------------------------------------------------------------------------------------------------
def run(star, t_end, C=0.4, dt_fixed=None, out_every_t=None, stop_on_crossing=True, max_steps=10**8, keep=None,
        rho_stop=None, verbose=False):
    """advance; returns list of diagnostic dicts (at output times and at the end)"""
    star.start()
    hist = [star.diagnostics(fresh=False)]
    if dt_fixed is None:
        dt0 = star.dt_cfl(star.geo, C)
    else:
        dt0 = dt_fixed
    next_out = out_every_t if out_every_t else None
    status = 'tmax'
    n = 0
    while star.t < t_end - 1e-12:
        dt = dt0 if dt_fixed is not None else star.dt_cfl(star.geo, C)
        dt = min(dt, t_end - star.t)
        star.step(dt); n += 1
        if np.any(np.diff(star.R) <= 0):
            status = 'crossing'
            hist.append(star.diagnostics(fresh=False)); break
        if next_out is not None and star.t >= next_out - 1e-12:
            hist.append(star.diagnostics(fresh=False)); next_out += out_every_t
        if rho_stop is not None and star.geo['e'][0] > rho_stop:
            status = 'rho_stop'; hist.append(star.diagnostics(fresh=False)); break
        if n >= max_steps:
            status = 'max_steps'; break
    if status == 'tmax' and (not hist or hist[-1]['t'] < star.t - 1e-12):
        hist.append(star.diagnostics(fresh=False))
    return hist, status, n


# ---------------------------------------------------------------------------------------------------------
# equilibrium solve / linearisation
def relax_equilibrium(star, tol=1e-13, maxit=30, verbose=False):
    """Newton iteration for the discrete static state F(R) = 0 at fixed zone contents dM (U = 0)."""
    N = star.N
    x = star.R[1:].copy()
    for it in range(maxit):
        F, a = star.force_R(x)
        res = np.max(np.abs(F))
        if verbose: print('  newton', it, res)
        if res < tol * 1.0: break
        J = fd_jacobian(star, x)
        dx = solve(J, -F)
        # damping to keep ordering
        lam = 1.0
        while np.any(np.diff(np.concatenate([[0.0], x + lam * dx])) <= 0) and lam > 1e-4:
            lam *= 0.5
        x = x + lam * dx
    star.R[1:] = x
    return res, it


def fd_jacobian(star, x, rel=1e-6):
    N = len(x)
    J = np.zeros((N, N))
    for i in range(N):
        d = rel * max(abs(x[i]), 1e-3 / N)
        xp = x.copy(); xp[i] += d
        xm = x.copy(); xm[i] -= d
        J[:, i] = (star.force_R(xp)[0] - star.force_R(xm)[0]) / (2 * d)
    return J


def linear_modes(star, nmodes=6):
    """eigenvalues omega^2 of the linearised semi-discrete system about the current (static) state:
       d^2 dR/dt^2 = diag(a) J dR  ->  omega^2 = -eig(diag(a) J).  Returns sorted omega^2 (ascending) and the eigenvectors"""
    x = star.R[1:].copy()
    J = fd_jacobian(star, x)
    F, a = star.force_R(x)
    A = a[:, None] * J
    w, v = eig(A)
    om2 = -w.real
    idx = np.argsort(om2)
    return om2[idx][:nmodes], v[:, idx][:, :nmodes].real, np.max(np.abs(w.imag))
