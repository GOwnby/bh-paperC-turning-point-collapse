"""
Spherically symmetric general-relativistic Lagrangian hydrodynamics (Misner & Sharp 1964; May & White 1966),
cold barotropic EOS (Chandrasekhar degenerate electrons, thermodynamically consistent energy density).

Lagrangian (rest-mass) coordinate mu. Interface variables r_j, U_j = dr/dtau, m_j (Misner-Sharp mass);
zone variables rho, P, e.  Exact relations used:
  Gamma^2 = 1 + U^2/c^2 - 2 G m/(r c^2)
  zone:  rho = Gamma_bar * dmu / V,  V=(4pi/3)(r_{j+1}^3-r_j^3)      (from dmu = 4 pi r^2 rho dr / Gamma)
  constraint: dm = e V / c^2                                        (from dm = 4 pi r^2 e dr)
  dU/dtau = -Gamma^2 c^2 (dP/dr)/(e+P) - G (m + 4 pi r^3 P/c^2)/r^2
  dm/dtau = -4 pi r^2 P U / c^2
Lapse: for a barotropic fluid in comoving slicing, alpha*h = const (h=(e+P)/rho), fixed by matching the
exterior Schwarzschild time at the surface: alpha_s = (1-2GM/(R c^2))/Gamma_s.  Time t = exterior (infinity-normalised)
coordinate time, the same normalisation as the Chandrasekhar pulsation eigenvalue sigma.
Well-balancing: the discrete static residual S_j of the initial equilibrium is subtracted (constant in time), the
same technique as the R0 subtraction of the homologous model.
"""
import numpy as np
from scipy.optimize import brentq
import background as bg
from metric import build_metric
G, c = bg.G, bg.c

class Lag:
    def __init__(self, rho_c, N=400, p=2.0, pert=None):
        self.rho_c = rho_c; self.N = N
        d = build_metric(rho_c); self.bgd = d
        r, rho, Mr = d['r'], d['rho'], d['Mr']
        dmu = 4*np.pi*r**2*rho/np.sqrt(1-2*G*Mr/(r*c**2))
        mu = np.concatenate([[0], np.cumsum(0.5*(dmu[1:]+dmu[:-1])*np.diff(r))])
        self.mu_tot = mu[-1]
        x = np.arange(N+1)/N
        self.mu = self.mu_tot*(1-(1-x)**p)           # finer toward the surface
        self.dmu = np.diff(self.mu)
        self.r0 = np.interp(self.mu, mu, r)
        self.R0 = self.r0[-1]
        self.M = d['M']
        self.pert = pert
        self.init_state()

    # ---- thermodynamics -------------------------------------------------
    @staticmethod
    def thermo(rho):
        P = bg.P_of_rho(np.maximum(rho, 1e-3)); P = np.where(rho > 1e-3, P, 0.0)
        e = bg.eps_over_c2_of_rho(np.maximum(rho, 1e-3))*c**2
        e = np.where(rho > 1e-3, e, rho*c**2)
        return P, e

    def solve_m(self, r, U):
        """integrate the Hamiltonian constraint dm = e V/c^2 outward for given r_j, U_j (Gamma implicit)"""
        N = self.N; m = np.zeros(N+1)
        for j in range(N):
            V = 4*np.pi/3*(r[j+1]**3-r[j]**3)
            def Gam(mj1):
                g0 = 1.0 if j == 0 else np.sqrt(1+U[j]**2/c**2-2*G*m[j]/(r[j]*c**2))
                g1 = np.sqrt(max(1+U[j+1]**2/c**2-2*G*mj1/(r[j+1]*c**2), 1e-12))
                return 0.5*(g0+g1)
            def F(mj1):
                rho = Gam(mj1)*self.dmu[j]/V
                P, e = self.thermo(np.array([rho]))
                return mj1-m[j]-e[0]*V/c**2
            guess = m[j]+self.dmu[j]*1.0
            m[j+1] = brentq(F, m[j], 1.5*guess+1e-3*self.mu_tot, xtol=1e-14*self.mu_tot, rtol=1e-15)
        return m

    def init_state(self):
        r = self.r0.copy(); U = np.zeros_like(r)
        if self.pert is not None:
            lam, lamdot = self.pert
            r = r*lam; U = lamdot*self.r0
        m = self.solve_m(r, U)
        self.y0 = np.concatenate([r[1:], U[1:], m[1:]])
        self.S = np.zeros(self.N)
        if self.pert is None:
            self.S = self.f(0.0, self.y0, subtract=False, bracket=True)
        else:
            # residual of the UNPERTURBED discrete equilibrium (same grid) is subtracted
            r0 = self.r0.copy(); U0 = np.zeros_like(r0); m0 = self.solve_m(r0, U0)
            y0 = np.concatenate([r0[1:], U0[1:], m0[1:]])
            self.S = self.f(0.0, y0, subtract=False, bracket=True)
        self.mtot = self.y0[-1]

    # ---- right-hand side ------------------------------------------------
    def unpack(self, y):
        N = self.N
        r = np.concatenate([[0.0], y[:N]]); U = np.concatenate([[0.0], y[N:2*N]]); m = np.concatenate([[0.0], y[2*N:]])
        return r, U, m

    def zone(self, r, U, m):
        Gam2 = np.ones_like(r)
        Gam2[1:] = 1+U[1:]**2/c**2-2*G*m[1:]/(r[1:]*c**2)
        Gam = np.sqrt(np.maximum(Gam2, 1e-12))
        Gb = 0.5*(Gam[1:]+Gam[:-1])
        V = 4*np.pi/3*(r[1:]**3-r[:-1]**3)
        rho = Gb*self.dmu/V
        P, e = self.thermo(rho)
        return Gam, rho, P, e

    def f(self, t, y, subtract=True, bracket=False):
        N = self.N
        r, U, m = self.unpack(y)
        Gam, rho, P, e = self.zone(r, U, m)
        enth = e+P
        h = enth/rho
        Ms = m[-1]
        alpha_s = (1-2*G*Ms/(r[-1]*c**2))/Gam[-1]
        if getattr(self,'freeze_alpha',None) is not None: alpha_s=self.freeze_alpha
        a_zone = alpha_s*c**2/h
        alpha = np.empty(N+1); alpha[1:N] = 0.5*(a_zone[1:]+a_zone[:-1]); alpha[N] = alpha_s; alpha[0] = a_zone[0]
        Pj = np.zeros(N+1); Pj[1:N] = 0.5*(P[1:]+P[:-1]); Pj[N] = 0.0; Pj[0] = P[0]
        enj = np.empty(N+1); enj[1:N] = 0.5*(enth[1:]+enth[:-1]); enj[N] = enth[-1]; enj[0] = enth[0]
        # mass-based (compatible) pressure force: a_P = -Gamma c^2 4 pi r^2 dP /(h_node M_node)
        dmu = self.dmu
        Mnode = np.empty(N); Mnode[:N-1] = 0.5*(dmu[:-1]+dmu[1:]); Mnode[N-1] = 0.5*dmu[-1]
        hz = h*dmu
        hnode = np.empty(N); hnode[:N-1] = 0.5*(hz[:-1]+hz[1:])/Mnode[:N-1]; hnode[N-1] = h[-1]
        dP = np.empty(N); dP[:N-1] = P[1:]-P[:-1]; dP[N-1] = 0.0-P[-1]
        acc = -Gam[1:]*c**2*4*np.pi*r[1:]**2*dP/(hnode*Mnode) - G*(m[1:]+4*np.pi*r[1:]**3*Pj[1:]/c**2)/r[1:]**2
        if subtract: acc = acc-self.S      # residual subtracted INSIDE the bracket so alpha-coupling stays consistent
        dU = alpha[1:]*acc
        if bracket: return acc
        dr_dt = alpha[1:]*U[1:]
        dm_dt = -alpha[1:]*4*np.pi*r[1:]**2*Pj[1:]*U[1:]/c**2
        return np.concatenate([dr_dt, dU, dm_dt])

    def dt_cfl(self, y, cfl=0.4):
        r, U, m = self.unpack(y)
        Gam, rho, P, e = self.zone(r, U, m)
        rho_s = np.maximum(rho, 1e-3)
        cs = c*np.sqrt(np.clip(bg.Gamma1_of_rho(rho_s)*P/(e+P), 1e-12, 1.0))
        dr = np.diff(r)
        h = (e+P)/rho
        alpha_s = (1-2*G*m[-1]/(r[-1]*c**2))/Gam[-1]
        al = alpha_s*c**2/h
        return cfl*np.min(dr/(al*cs))

    def run(self, T, cfl=0.4, nout=400, events=None, y0=None, t0=0.0):
        y = self.y0.copy() if y0 is None else y0.copy(); t = t0
        out = []; nextout = t0; step = 0
        while t < T:
            if t >= nextout:
                out.append(self.diag(t, y)); nextout += (T-t0)/nout
                if events is not None and events(t, y, out[-1]): break
            dt = self.dt_cfl(y, cfl)
            k1 = self.f(t, y); k2 = self.f(t+dt/2, y+dt/2*k1); k3 = self.f(t+dt/2, y+dt/2*k2); k4 = self.f(t+dt, y+dt*k3)
            y = y+dt/6*(k1+2*k2+2*k3+k4); t += dt; step += 1
        out.append(self.diag(t, y))
        return np.array(out), y, t, step

    def diag(self, t, y):
        r, U, m = self.unpack(y)
        Gam, rho, P, e = self.zone(r, U, m)
        # half-rest-mass radius ratio, central density ratio, max compactness, constraint violation
        jh = np.searchsorted(self.mu, 0.5*self.mu_tot)
        lam_h = r[jh]/self.r0[jh]
        comp = np.max(2*G*m[1:]/(r[1:]*c**2))
        V = 4*np.pi/3*(r[1:]**3-r[:-1]**3)
        viol = np.max(np.abs(np.diff(m)-e*V/c**2))/m[-1]
        return [t, lam_h, rho[0]/self.rho_c, comp, viol, r[-1]/self.R0]
