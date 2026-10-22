"""General-relativistic linear radial pulsation spectrum of a TOV star, independent of the Lagrangian scheme.

Equations: Chandrasekhar (1964, ApJ 140, 417), in the Sturm-Liouville form of Gondek, Haensel & Zdunik (1997, A&A 325, 217):
   ( Pi z' )' + ( Q + omega^2 Wt ) z = 0,           z = r^3 xi e^{-nu/2},  xi = Delta r / r
   Pi = e^{(3 nu + lambda)/2} Gamma1 p / r^2
   Q  = e^{(3 nu + lambda)/2} r^-2 [ p'^2/(e+p) - 4 p'/r - 8 pi e^lambda (e+p) p ]
   Wt = e^{(nu + 3 lambda)/2} (e+p) / r^2
with e^lambda = 1/(1-2m/r), e^nu = (1-2M/R)/h^2 (a h = const, exterior-matched), omega the frequency wrt the
Schwarzschild time t_inf.  Boundary conditions: z ~ r^3 at the centre, Pi z' -> 0 at the surface (Delta p = 0).
The equation is re-derived in the report-independent note derivations.md from the first-order form (xi, Delta p).
The weak form  int [Pi z'^2 - Q z^2] = omega^2 int Wt z^2  is solved with linear finite elements on a graded mesh,
3-point Gauss quadrature with exact (not interpolated) background values.
omega_mine^2 = omega_inf^2 / (1 - 2M/R)  converts to the proper time of the surface shell used in lag.py.
"""
import numpy as np
from scipy.linalg import eigh
import structure as S

GL3 = (np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)]), np.array([5 / 9, 8 / 9, 5 / 9]))


def background(eos, st, r):
    m, H, Mb, rN, e, p = S.tov_at(st, eos, r)
    _, _, g1 = eos.thermo(rN)
    f = 1 - 2 * m / r
    E = e + p
    dp = -E * (m + 4 * np.pi * r**3 * p) / (r**2 * f)
    h = np.exp(H)
    enu = (1 - 2 * st['M'] / st['R']) / h**2
    elam = 1 / f
    Pi = enu**1.5 * np.sqrt(elam) * g1 * p / r**2
    Q = enu**1.5 * np.sqrt(elam) / r**2 * (dp**2 / E - 4 * dp / r - 8 * np.pi * elam * E * p)
    Wt = np.sqrt(enu) * elam**1.5 * E / r**2
    return Pi, Q, Wt


def modes(eos, st, n=1000, beta=2.0, nmodes=6):
    x = np.linspace(0, 1, n + 1)
    r = st['R'] * (1 - (1 - x)**beta)
    r[0] = 0.0
    K = np.zeros((n + 1, n + 1)); Mm = np.zeros((n + 1, n + 1))
    xg, wg = GL3
    a, b = r[:-1], r[1:]
    h = b - a
    rg = 0.5 * (a[:, None] + b[:, None]) + 0.5 * h[:, None] * xg[None, :]
    Pi, Q, Wt = background(eos, st, rg.ravel())
    Pi, Q, Wt = [v.reshape(rg.shape) for v in (Pi, Q, Wt)]
    wq = 0.5 * h[:, None] * wg[None, :]
    t = (rg - a[:, None]) / h[:, None]
    phi = [1 - t, t]
    dphi = [-1 / h[:, None] * np.ones_like(t), 1 / h[:, None] * np.ones_like(t)]
    for i in range(2):
        for j in range(2):
            kij = np.sum(wq * (Pi * dphi[i] * dphi[j] - Q * phi[i] * phi[j]), axis=1)
            mij = np.sum(wq * Wt * phi[i] * phi[j], axis=1)
            np.add.at(K, (np.arange(n) + i, np.arange(n) + j), kij)
            np.add.at(Mm, (np.arange(n) + i, np.arange(n) + j), mij)
    K = K[1:, 1:]; Mm = Mm[1:, 1:]
    w, v = eigh(K, Mm, subset_by_index=[0, nmodes - 1])
    fac = 1.0 / (1 - 2 * st['M'] / st['R'])
    return w[:nmodes] * fac, v[:, :nmodes], r[1:]


if __name__ == '__main__':
    from initdata import Setup
    for kind in ('neu', 'pure'):
        Sx = Setup(kind)
        s = (2.998e8 / Sx.L)**2
        for n in (250, 500, 1000, 2000):
            w, v, r = modes(Sx.eos, Sx.tov, n=n)
            print(kind, n, np.array2string(w * s, precision=4))
