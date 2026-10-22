"""Barotropic equation-of-state model built from the report's EOS routines (copies in ../ref/, unmodified).

The report's closure: p = P(rho), where 'rho' is used in the Misner-Sharp equations as the total
mass-energy density e (m' = 4 pi e R^2 R', report Eq. (mr); the report never distinguishes rest-mass and
energy density).  We keep exactly that p(e) relation (report functions IP.P_of_rho / EN.P_of_rho) and add the
thermodynamics that a barotropic (isentropic, cold) fluid must obey and that a Lagrangian scheme needs:

   first law at fixed entropy:  d e = (e + p) d rho_N / rho_N      (rho_N = conserved number-density function)
   => ln rho_N(e) = ln e - int_0^e  p /(e'+p) d ln e'       (normalised so rho_N -> e as e -> 0)
   => enthalpy per particle h = (e+p)/rho_N,  H = ln h = int dp/(e+p),   adiabatic index Gamma1 = d ln p / d ln rho_N,
   => sound speed squared (units of c^2)  c_s^2 = dp/de = Gamma1 p /(e+p).

rho_N is the quantity whose zone-integral (baryon number in zone) is a Lagrangian invariant of the evolution.
For the report's EOS p/e ~ 1e-3 so rho_N/e - 1 ~ -1e-3: the distinction is a 0.1% effect but is needed for an
exactly conservative Lagrangian scheme.  Newtonian mode (newt=True): rho_N := e := rho, h := 1 (c -> infinity).

Units in the dynamical code: G = c = 1, length unit L, so  e_d = G L^2 e/c^2 ,  p_d = G L^2 p/c^4  (as in the report).
"""
import os, sys, io, contextlib
import numpy as np
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, '..', 'ref')
CACHE = os.path.join(HERE, '..', 'results')

G_SI = 6.674e-11; C_SI = 2.998e8; MSUN = 1.989e30   # as eos.py of the report


def _report_modules(kind):
    sys.path.insert(0, REF)
    with contextlib.redirect_stdout(io.StringIO()):
        import eos as E
        if kind == 'pure':
            import initial_profile as IP
            P = IP.P_of_rho
        else:
            import eos_neutronized as EN
            import initial_profile_neu as IP
            P = EN.P_of_rho
    return E, IP, P


def build_table(kind, n=6001, lo=1e4, hi=1e16):
    """Tabulate the report's P(e) for e >= lo = 1e4 kg/m^3 and derive rho_N, h (SI).  Below lo a pure power law p = K e^(5/3) (Gamma = 5/3)
    anchored at lo is used: the report's own tabulated P(rho) is polluted by cancellation noise (slope errors of several %) in f_cs(x) for x < 1e-3
    (rho < ~1e3 kg/m^3), see logs/eos_checks.log."""
    E, IP, P = _report_modules(kind)
    e = np.geomspace(lo, hi, n)
    p = np.atleast_1d(P(e)).astype(float)
    lne = np.log(e)
    f = p / (e * C_SI**2 + p)                         # integrand of c(e) wrt ln e  (e = mass density, energy density = e c^2)
    # c(e) = int_0^e f dln e' ; for e<lo: p ~ K e^(5/3) (non-relativistic power law of both modules) -> f = p/e ~ e^(2/3)
    # int_0^lo f dln e' = (3/2) f(lo)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(lne))]) + 1.5 * f[0]
    rhoN = e * np.exp(-cum)
    return dict(e=e, p=p, rhoN=rhoN, cum=cum)


def get_table(kind):
    fn = os.path.join(CACHE, f'eos_table_{kind}.npz')
    if os.path.exists(fn):
        d = np.load(fn)
        return {k: d[k] for k in d.files}
    d = build_table(kind)
    os.makedirs(CACHE, exist_ok=True)
    np.savez(fn, **d)
    return d


class EOS:
    """Barotropic EOS in dimensionless (G=c=1, length L) units.  scale: e_d = Ge*e_SI, p_d = Gp*p_SI."""

    def __init__(self, kind, L, newt=False, p_floor=0.0, table=None):
        t = table if table is not None else get_table(kind)
        self.kind, self.L, self.newt = kind, L, newt
        self.Ge = G_SI * L**2 / C_SI**2
        self.Gp = G_SI * L**2 / C_SI**4
        e, p, rN = t['e'] * self.Ge, t['p'] * self.Gp, t['rhoN'] * self.Ge
        if newt:
            rN = e.copy()
        self.e_tab, self.p_tab, self.rN_tab = e, p, rN
        x = np.log(rN)
        self._le = CubicSpline(x, np.log(e))     # ln e   as function of ln rho_N
        self._lp = CubicSpline(x, np.log(p))     # ln p
        self._dlp = self._lp.derivative()         # Gamma1 = dln p/dln rho_N
        self.x_lo, self.x_hi = x[0], x[-1]
        # NR power-law tail below the table: p = K rho_N^(5/3) (e = rho_N)
        self.K_nr = p[0] / rN[0]**(5.0 / 3.0)
        self.p_floor = p_floor * 1.0
        # enthalpy H = ln((e+p)/rho_N) as function of ln H (for the TOV solver)
        H = (t['cum'] + np.log1p(p / e)) if not newt else np.log1p(p / e)   # = ln((e+p)/rho_N), accurate at low density
        self.H_tab = H
        if not newt:
            ok = H > 0
            lH = np.log(H[ok])
            self._lH_x = CubicSpline(lH, x[ok])       # ln rho_N as function of ln H

    # -- thermodynamics at given number density rho_N (dimensionless, array) -------------------------------
    def thermo(self, rN):
        """returns e, p, Gamma1 as functions of rho_N"""
        x = np.log(np.maximum(rN, 1e-300))
        low = x < self.x_lo
        xc = np.where(low, self.x_lo, x)
        e = np.exp(self._le(xc)); p = np.exp(self._lp(xc)); g1 = self._dlp(xc)
        if np.any(low):
            e = np.where(low, rN, e)
            p = np.where(low, self.K_nr * np.maximum(rN, 0.0)**(5.0 / 3.0), p)
            g1 = np.where(low, 5.0 / 3.0, g1)
        if self.newt:
            e = rN
        return e, p, g1

    def e_of_rN(self, rN): return self.thermo(rN)[0]

    def rN_of_e(self, e):
        """inverse of e(rho_N): table interpolation (monotone), polished by Newton steps on thermo()"""
        e = np.asarray(e, float)
        if self.newt:
            return e.copy()
        ee = np.maximum(e, 1e-300)
        rN = np.exp(np.interp(np.log(ee), np.log(self.e_tab), np.log(self.rN_tab)))
        low = ee < self.e_tab[0]
        rN = np.where(low, ee, rN)
        for _ in range(3):
            et, pt, _g = self.thermo(rN)
            rN = rN * (ee / et)**(1.0 / (1.0 + pt / np.maximum(et, 1e-300)))
        return rN

    def state_from_H(self, H):
        """(rho_N, e, p) as function of enthalpy H>0 (for the stellar-structure integrations)"""
        H = np.maximum(np.asarray(H, float), 1e-300)
        lH = np.log(H)
        lH0 = np.log(self.H_tab[self.H_tab > 0][0])
        x = np.where(lH >= lH0, self._lH_x(np.maximum(lH, lH0)), 0.0)
        # below the table (NR tail): H = (5/2) K rho_N^(2/3)  (p = K rho^(5/3), e = rho, H ~ int dp/(rho) )
        rN_low = (H / (2.5 * self.K_nr))**1.5
        rN = np.where(lH >= lH0, np.exp(x), rN_low)
        e, p, _ = self.thermo(rN)
        return rN, e, p


def kind_L(kind, rho_c=5e12):
    """The report's stellar radius (its Newtonian build_star) = length unit L of the report for this EOS.
    rmax = 2e6 m, n = 8000 (collapse2.py) for rho_c = 5e12; for rho_c < 1e12 the star is larger than 2e6 m (R ~ 4.8e6 m at 1e11) and
    the report's linear/Eulerian control scripts, which use rmax = 3e6 m, then silently TRUNCATE the star (see logs/control_star_truncation.log);
    here rmax is enlarged with the same step (250 m)."""
    E, IP, P = _report_modules(kind)
    rmax = 2e6 if rho_c >= 1e12 else 2.5e7
    sol = IP.build_star(rho_c, rmax=rmax, n=int(rmax / 250))
    assert sol.status == 1, 'surface event not reached: increase rmax'
    return sol.t[-1], sol


if __name__ == '__main__':
    for kind in ('pure', 'neu'):
        t = get_table(kind)
        rN = t['rhoN']; e = t['e']
        print(kind, 'table ok; rho_N/e at e=5e12:', np.interp(5e12, e, rN / e), 'at 1e10:', np.interp(1e10, e, rN / e))
