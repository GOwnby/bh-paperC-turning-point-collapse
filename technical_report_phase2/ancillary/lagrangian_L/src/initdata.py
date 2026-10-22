"""Initial-data builders (dimensionless units, length unit L = the report's R_star for the EOS)."""
import numpy as np
import eosmodel as M, structure as S
from lag import Lag


class Setup:
    """EOS + exact TOV solution + the report's Newtonian profile for rho_c = 5e12 kg/m^3 (cached per kind)."""
    _cache = {}

    def __init__(self, kind, rho_c=5e12, newt=False):
        key = (kind, rho_c, newt)
        if key in Setup._cache:
            self.__dict__.update(Setup._cache[key].__dict__); return
        self.kind, self.rho_c, self.newt = kind, rho_c, newt
        self.prof = S.newtonian_profile(kind, rho_c)
        self.L = self.prof['L']
        self.eos = M.EOS(kind, self.L, newt=newt)
        self.tov = S.tov(self.eos, rho_c) if not newt else None
        Setup._cache[key] = self

    def tov_data(self, N):
        st = self.tov
        r = np.linspace(0.0, st['R'], N + 1)
        m, H, Mb, rN, e, p = S.tov_at(st, self.eos, r)
        return dict(R0=r, dM=np.diff(Mb), m=m)

    def newton_data(self, N, newt=False, linear=False):
        """The report's initial state (R_j = r_j, m_j = Newtonian enclosed mass, U = 0), mapped to zone contents."""
        r = np.linspace(0.0, 1.0, N + 1)
        m = self.prof['m_of_r'](r, linear=linear); m[0] = 0.0
        V = 4 * np.pi / 3 * (r[1:]**3 - r[:-1]**3)
        e = np.diff(m) / V
        if newt:
            return dict(R0=r, dM=np.diff(m), m=m)
        rN = self.eos.rN_of_e(e)
        G = np.ones(N + 1); G[1:] = np.sqrt(1 - 2 * m[1:] / r[1:])
        Gz = 0.5 * (G[:-1] + G[1:])
        return dict(R0=r, dM=rN * V / Gz, m=m)

    def make(self, data, **kw):
        return Lag(self.eos, data['dM'], data['R0'], newt=self.newt, **kw)
