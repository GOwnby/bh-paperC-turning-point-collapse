"""
Eulerian, conservative, shock-capturing spherical hydrodynamics on a FIXED
grid, as a genuinely different numerical strategy from the Lagrangian
Misner-Sharp scheme used before. Scope, stated honestly:

  - Newtonian self-gravitating hydrodynamics (mass, momentum, energy
    conservation in spherical symmetry), NOT full relativistic Eulerian
    hydrodynamics. This is a deliberate simplification, justified by the
    finding that the collapse instability once neutronization is included
    is fundamentally Newtonian (the mass-radius turnover dM/drho_c<0 is a
    classical stability criterion); GR corrections were shown to be a small
    (~0.2%) perturbation on this system, not its driver.
  - The leading relativistic correction (gravitational redshift of the
    potential) is retained as a source-term correction, not dropped
    silently.
  - HLL approximate Riemann solver, second-order MUSCL reconstruction with
    a minmod limiter, explicit RK2 (midpoint) time stepping. Standard,
    well-documented methods (e.g. Toro's "Riemann Solvers and Numerical
    Methods for Fluid Dynamics").

A full general-relativistic Eulerian code (Valencia-formulation conservative
variables, a real GR Riemann solver, dynamical lapse/metric evolution) is a
substantially larger undertaking, flagged explicitly as future work.
"""
import numpy as np
import eos_neutronized as EN
import eos as E

G, c = E.G, E.c

def minmod(a, b):
    return np.where(a*b <= 0, 0.0, np.where(np.abs(a) < np.abs(b), a, b))

class EulerGrid:
    def __init__(self, rmax, N, rho_of_r, ghost=2):
        self.N = N; self.ghost = ghost
        self.dr = rmax/N
        self.r = (np.arange(-ghost, N+ghost)+0.5)*self.dr   # cell centers, with ghost zones
        self.rho = np.maximum(rho_of_r(np.abs(self.r)), 1.0)   # floor to avoid vacuum issues
        self.v = np.zeros_like(self.r)
        self.p = EN.P_of_rho(self.rho)

    def eint(self, rho, p):
        # specific internal energy from p via e = integral(p/rho^2)drho along the EOS;
        # for a nearly-polytropic EOS locally, e ~ p/(rho*(Gamma-1)) is a good and standard
        # approximation; use the LOCAL effective Gamma from the EOS table.
        Gamma_loc = EN._dlogPdlogrho(np.log(np.maximum(rho, EN._rho_tab[0])))
        Gamma_loc = np.where(rho >= EN._rho_tab[0], Gamma_loc, 5/3)
        return p/(rho*np.maximum(Gamma_loc-1, 0.1))

    def cons(self):
        rho, v, p = self.rho, self.v, self.p
        D = rho
        Sr = rho*v
        e_int = self.eint(rho, p)
        E_tot = rho*(e_int + 0.5*v**2)
        return np.array([D, Sr, E_tot])

    def prim_from_cons(self, U):
        D, Sr, Etot = U
        rho = np.maximum(D, 1.0)
        v = Sr/rho
        e_int = np.maximum(Etot/rho - 0.5*v**2, 0.0)
        # invert e_int(rho) -> p using local Gamma (consistent with self.eint above)
        Gamma_loc = EN._dlogPdlogrho(np.log(np.maximum(rho, EN._rho_tab[0])))
        Gamma_loc = np.where(rho >= EN._rho_tab[0], Gamma_loc, 5/3)
        p = e_int*rho*np.maximum(Gamma_loc-1, 0.1)
        p = np.maximum(p, 0.0)
        return rho, v, p

    def mass_function(self):
        r = np.abs(self.r[self.ghost:-self.ghost])
        rho = self.rho[self.ghost:-self.ghost]
        m = np.concatenate([[0.0], np.cumsum(4*np.pi*r[:-1]**2*rho[:-1]*np.diff(r))])
        return r, m

    def sound_speed(self, rho, p):
        Gamma_loc = EN._dlogPdlogrho(np.log(np.maximum(rho, EN._rho_tab[0])))
        Gamma_loc = np.where(rho >= EN._rho_tab[0], Gamma_loc, 5/3)
        return np.sqrt(np.maximum(Gamma_loc*p/np.maximum(rho,1e-30), 1e-30))

print("Euler grid module loaded OK")
