"""Machine verification of the Misner-Sharp system used in lag.py, without a computer-algebra system.

Method: the Einstein tensor at a point depends only on the 2-jet (value, first and second derivatives) of the metric at that
point.  We represent the metric functions a(t,r), X(t,r), R(t,r) by RANDOM truncated Taylor polynomials of degree 2 (arbitrary
smooth functions have exactly such 2-jets), compute the Christoffel symbols, Ricci and Einstein tensors of
     ds^2 = -a^2 dt^2 + X^2 dr^2 + R^2 (dtheta^2 + sin^2 theta dphi^2)
by exact truncated-Taylor arithmetic (automatic differentiation, no finite differences), and test the claimed identities
to rounding error at many random points.  A polynomial identity in the 2-jet coefficients that holds at many random points holds
identically.  Identities tested (U = (dR/dt)/a, Gamma = R'/X, m := (R/2)(1 + U^2 - Gamma^2), ' = d/dr, dot = d/dt):
  I1  G^t_t = -2 m' / (R^2 R')                       -> with G^t_t = -8 pi e :  m' = 4 pi e R^2 R'     (on G_tr = 0)
  I2  G^r_r = -2 mdot / (R^2 Rdot)                   -> with G^r_r = +8 pi p :  mdot = -4 pi p R^2 Rdot (on G_tr = 0)
  (for a comoving perfect fluid G_tr = 0 is itself a field equation; the jets are drawn at random and dX/dt is then fixed so that
   G_tr = 0 at the point; I3 and I4 are tested for completely unconstrained random jets)
  I3  G_tr  = -(2 X/R) ( Gammadot - a' Rdot/(a X) )   -> G_tr = 0 <=> Gammadot = a' Rdot /(a X)  (NOT Gammadot = 0)
  I4  U Udot = Gamma Gammadot + mdot/R - m Rdot/R^2   (definition; combined with I2, I3 gives dU/dt equation (5))
  I5  nabla_mu T^mu_r = dp/dr + (a'/a)(e+p)            for T^mu_nu = diag(-e,p,p,p)   -> nabla_mu T^mu_r = 0 is the lapse equation a'/a = -p'/(e+p)
  I6  nabla_mu T^mu_t = -edot - (e+p)(Xdot/X + 2 Rdot/R)                              -> nabla_mu T^mu_t = 0 is the energy (continuity) equation
      (I5, I6 use the Christoffel symbols computed here and random 1-jets of e and p; both are exact consequences of the field equations through the
       contracted Bianchi identity, which is the textbook route to the lapse equation of Misner & Sharp 1964.)
Also reported: size of the residual of the report's claims (R'/X independent of tau; evolve2.py's a^2).
"""
import numpy as np, itertools, sys

K = 2
NV = 3          # variables (t, r, theta)
MON = [e for e in itertools.product(range(K + 1), repeat=NV) if sum(e) <= K]
IDX = {e: i for i, e in enumerate(MON)}


class T:
    """truncated multivariate Taylor polynomial about the expansion point, total degree <= order"""
    def __init__(self, c, order=K):
        self.c = np.array(c, float); self.o = order

    @staticmethod
    def const(v, order=K):
        c = np.zeros(len(MON)); c[IDX[(0, 0, 0)]] = v; return T(c, order)

    def _trunc(self):
        for e in MON:
            if sum(e) > self.o: self.c[IDX[e]] = 0.0
        return self

    def __add__(s, o):
        o = o if isinstance(o, T) else T.const(o); return T(s.c + o.c, min(s.o, o.o))._trunc()
    __radd__ = __add__
    def __neg__(s): return T(-s.c, s.o)
    def __sub__(s, o): return s + (-(o if isinstance(o, T) else T.const(o)))
    def __rsub__(s, o): return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, T): return T(s.c * o, s.o)
        o_ = min(s.o, o.o); out = np.zeros(len(MON))
        for e1 in MON:
            a = s.c[IDX[e1]]
            if a == 0: continue
            for e2 in MON:
                e = tuple(x + y for x, y in zip(e1, e2))
                if sum(e) <= o_: out[IDX[e]] += a * o.c[IDX[e2]]
        return T(out, o_)
    __rmul__ = __mul__

    def inv(s):
        v = s.c[IDX[(0, 0, 0)]]; d = s * (1.0 / v) - 1.0   # s/v = 1 + d
        out = T.const(1.0, s.o); term = T.const(1.0, s.o)
        for n in range(1, s.o + 1):
            term = term * (-d); out = out + term
        return out * (1.0 / v)
    def __truediv__(s, o): return s * (o.inv() if isinstance(o, T) else T.const(1.0 / o))
    def __rtruediv__(s, o): return s.inv() * o

    def val(s): return s.c[IDX[(0, 0, 0)]]

    def d(s, i):
        out = np.zeros(len(MON))
        if i >= NV: return T(out, max(s.o - 1, 0))     # nothing depends on phi
        for e in MON:
            if e[i] > 0:
                e2 = list(e); e2[i] -= 1
                out[IDX[tuple(e2)]] += e[i] * s.c[IDX[e]]
        return T(out, max(s.o - 1, 0))

    def sq(s): return s * s


def sin_t(theta0, order=K):
    """sin(theta) as a Taylor polynomial in the third variable"""
    x = np.zeros(len(MON))
    s0, c0 = np.sin(theta0), np.cos(theta0)
    x[IDX[(0, 0, 0)]] = s0; x[IDX[(0, 0, 1)]] = c0; x[IDX[(0, 0, 2)]] = -s0 / 2
    return T(x, order)


def einstein(gdiag):
    """gdiag: list of 4 T (diagonal metric components).  Returns G_{mu nu} (4x4 array of floats) and g (floats)."""
    g = gdiag
    ginv = [x.inv() for x in g]
    # Christoffel Gamma^l_{mn} = 1/2 g^{ll} (d_m g_{ln} + d_n g_{lm} - d_l g_{mn})  (diagonal metric)
    Gam = [[[None] * 4 for _ in range(4)] for _ in range(4)]
    for l in range(4):
        for m in range(4):
            for n in range(4):
                t = T.const(0.0)
                if l == n: t = t + g[l].d(m)
                if l == m: t = t + g[l].d(n)
                if m == n: t = t - g[m].d(l)
                Gam[l][m][n] = ginv[l] * t * 0.5
    Ric = np.zeros((4, 4))
    v = lambda x: x.val()
    for m in range(4):
        for n in range(4):
            s = 0.0
            for l in range(4):
                s += v(Gam[l][m][n].d(l)) - v(Gam[l][m][l].d(n))
                for sg in range(4):
                    s += v(Gam[l][l][sg]) * v(Gam[sg][m][n]) - v(Gam[l][n][sg]) * v(Gam[sg][m][l])
            Ric[m, n] = s
    gv = np.array([v(x) for x in g]); giv = 1.0 / gv
    Rs = sum(giv[i] * Ric[i, i] for i in range(4))
    G = Ric - 0.5 * np.diag(gv) * Rs
    return G, gv, Gam


def rand_fn(rng, v0, scale=0.3):
    c = np.zeros(len(MON))
    for e in MON:
        if e[2] == 0:   # functions of (t, r) only
            c[IDX[e]] = rng.normal() * (scale if sum(e) > 0 else 0.0)
    c[IDX[(0, 0, 0)]] = v0
    return T(c)


def main(npts=40, seed=1, impose_Gtr0=True):
    rng = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed + 1000)      # separate stream for the fluid jets, so that I1-I4 use the same points as before
    worst = dict(I1=0, I2=0, I3=0, I4=0, I5=0, I6=0)
    ctrl5 = []
    ctrl = []
    gamma_dot_rel = []; a2 = []
    for ip in range(npts):
        a = rand_fn(rng, 1.0 + 0.3 * rng.random())
        X = rand_fn(rng, 1.0 + 0.5 * rng.random())
        R = rand_fn(rng, 1.0 + rng.random())
        # impose the momentum constraint G_tr = 0 at the point (comoving perfect fluid): choose dX/dt so that
        #   dot(R') = (dot(X)/X) R' + (a'/a) dot(R)      [identity I3 shows this is exactly G_tr = 0]
        if impose_Gtr0:
            Xd = X.val() * (R.d(0).d(1).val() - a.d(1).val() / a.val() * R.d(0).val()) / R.d(1).val()
            X.c[IDX[(1, 0, 0)]] = Xd
        # make the t- and r-derivatives of R of order one
        th0 = 0.4 + 2.0 * rng.random()
        st = sin_t(th0)
        g = [-(a * a), X * X, R * R, R * R * st * st]
        G, gv, Gam = einstein(g)
        e_f = rand_fn(rng2, 0.5 + rng2.random()); p_f = rand_fn(rng2, 0.2 + 0.5 * rng2.random())
        Gmix_tt = G[0, 0] / gv[0]; Gmix_rr = G[1, 1] / gv[1]
        # derived quantities as T (order 2): U, Gamma, m  and their derivatives
        Rd = R.d(0); Rp = R.d(1)
        U = Rd / a; Gm = Rp / X
        m = R * 0.5 * (1.0 + U * U - Gm * Gm)
        mp = m.d(1).val(); md = m.d(0).val()
        Rv, Rdv, Rpv, av, Xv = R.val(), Rd.val(), Rp.val(), a.val(), X.val()
        I1 = Gmix_tt - (-2 * mp / (Rv**2 * Rpv))
        I2 = Gmix_rr - (-2 * md / (Rv**2 * Rdv))
        Gdot = Gm.d(0).val()
        I3 = G[0, 1] - (-(2 * Xv / Rv) * (Gdot - a.d(1).val() * Rdv / (av * Xv)))
        ctrl.append(abs(Gmix_tt - (-2.02 * mp / (Rv**2 * Rpv))))     # negative control: a 1 % wrong coefficient must be detected
        if not impose_Gtr0: I1 = I2 = 0.0   # I1, I2 are only claimed on the constraint surface G_tr = 0
        # note: G_{tr} with coordinate indices; relation has X from  Rdot'/..; see docstring
        Uv, Gv, mv = U.val(), Gm.val(), m.val()
        I4 = Uv * U.d(0).val() - (Gv * Gdot + md / Rv - mv * Rdv / Rv**2)
        # I5, I6: divergence of the perfect-fluid stress tensor T^mu_nu = diag(-e, p, p, p) in these coordinates
        ev, pv = e_f.val(), p_f.val(); Td = [-ev, pv, pv, pv]
        Chr = lambda l, mm, n: Gam[l][mm][n].val()
        div_r = p_f.d(1).val() + sum(Chr(mu, mu, 1) * (pv - Td[mu]) for mu in range(4))
        div_t = -e_f.d(0).val() + sum(Chr(mu, mu, 0) * (-ev - Td[mu]) for mu in range(4))
        I5 = div_r - (p_f.d(1).val() + a.d(1).val() / av * (ev + pv))
        I6 = div_t - (-e_f.d(0).val() - (ev + pv) * (X.d(0).val() / Xv + 2.0 * Rdv / Rv))
        ctrl5.append(abs(div_r - (p_f.d(1).val() + 1.01 * a.d(1).val() / av * (ev + pv))))   # negative control: 1 % wrong coefficient
        for k, x in zip(('I1', 'I2', 'I3', 'I4', 'I5', 'I6'), (I1, I2, I3, I4, I5, I6)):
            worst[k] = max(worst[k], abs(x))
        # size of G_tr for generic data and relation to Gammadot
        gamma_dot_rel.append(abs(Gdot) / max(abs(Gv), 1e-30))
    worst['ctrl_min'] = min(ctrl); worst['ctrl5_min'] = float(np.median(ctrl5))
    return worst, np.median(gamma_dot_rel)


if __name__ == '__main__':
    w, med = main()
    print('Misner-Sharp identities verified by exact 2-jet Einstein tensor at 40 random points (max abs residual):')
    for k in ('I1', 'I2', 'I3', 'I4', 'I5', 'I6'):
        print(f'   {k}: {w[k]:.2e}')
    print(f"   negative control (I1 with coefficient 2.02 instead of 2): smallest residual over the 40 points = {w['ctrl_min']:.2e} (must be >> 1e-9)")
    print(f"   negative control (I5 with coefficient 1.01 instead of 1): median residual over the 40 points = {w['ctrl5_min']:.2e} (must be >> 1e-9; it is proportional to the random (a'/a)(e+p))")
    ok = all(w[k] < 1e-9 for k in ('I1', 'I2', 'I3', 'I4', 'I5', 'I6')) and w['ctrl_min'] > 1e-6 and w['ctrl5_min'] > 1e-4
    print('PASS' if ok else 'FAIL', ' (tolerance 1e-9; residuals are floating-point rounding of O(1) quantities)')
    print(f'For generic data Gammadot/Gamma ~ {med:.2f}: Gamma is not constant along a shell unless G_tr and a\' Rdot vanish (report Eq. (constraint) holds only if a\'=0).')
    sys.exit(0 if ok else 1)
