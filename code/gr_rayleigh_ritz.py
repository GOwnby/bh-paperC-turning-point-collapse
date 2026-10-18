"""
Phase A2 part 2: Chandrasekhar's (1964) GR variational radial-pulsation
equation, Rayleigh-Ritz solve (reproducing Paper C's validated sigma^2=-41.05
s^-2), this time also extracting the eigenVECTOR (not just the eigenvalue)
so we have the true fundamental eigenfunction to seed the nonlinear code with.

Variational (Rayleigh-quotient) form, Eq. (chandra-variational) in Paper C:
  sigma^2 * Integral e^{-a+3b}(1+P/(rho c^2)) rho r^2 xi^2 dr
    = Integral Gamma1 P e^{3a+b}/r^2 [d/dr(r^2 e^{-a} xi)]^2 dr
      + 4 Integral P' r e^{a+b} xi^2 dr
      - Integral P'^2 r^2 e^{a+b} / (P+rho c^2) xi^2 dr

Define zeta_i = e^{-a} r^2 xi_i.  Then, with
  Pi(r) = Gamma1 * P * e^{3a+b} / r^2          (coefficient of zeta_i' zeta_j')
  Q(r)  = 4 P' r e^{a+b} - P'^2 r^2 e^{a+b}/(P+rho c^2)   (coefficient of xi_i xi_j)
  S(r)  = e^{-a+3b} (1+P/(rho c^2)) rho r^2      (coefficient of xi_i xi_j, LHS)

  A_ij = Integral [ Pi * zeta_i' zeta_j' + Q * xi_i * xi_j ] dr
  B_ij = Integral S * xi_i * xi_j dr
  Generalized eigenproblem: A v = sigma^2 B v.

Basis: xi_1 = r,  xi_k = r*(1-(r/R)^2)^(k-1) for k=2..5 -- same
regularity-respecting basis as Paper C Section 5.8.
"""
import numpy as np
from scipy.linalg import eigh
import background as bg
from metric import build_metric

G, c = bg.G, bg.c
INCLUDE_8PI = True

def make_basis(R, n_basis=5):
    basis = []
    basis.append(lambda r: r)
    basis.append(lambda r: np.ones_like(r))  # derivative helper not needed; use numeric diff
    # Actually build as explicit list of callables xi_k(r):
    fns = [lambda r: r]
    for k in range(1, n_basis):
        fns.append(lambda r, k=k: r*(1-(r/R)**2)**k)
    return fns[:n_basis]

def build_RR_matrices(bgd, n_basis=5):
    r = bgd['r']; rho = bgd['rho']; P = bgd['P']; a = bgd['a']; b = bgd['b']
    R = bgd['R']; dPdr = bgd['dPdr']

    ea, eb = np.exp(a), np.exp(b)
    Gamma1 = np.array([bg.Gamma1_of_rho(rr) if rr>1 else 5/3 for rr in rho])

    Pi = Gamma1*P*np.exp(3*a+b)/r**2
    eps = bgd['eps']
    Q  = 4*dPdr*r*np.exp(a+b) - dPdr**2*r**2*np.exp(a+b)/(P+eps*c**2)
    if INCLUDE_8PI: Q = Q + 8*np.pi*G/c**4*(P+eps*c**2)*P*r**2*np.exp(a+3*b)
    S  = np.exp(-a+3*b)*(1+P/(eps*c**2))*eps*r**2

    fns = make_basis(R, n_basis)
    xis = [f(r) for f in fns]
    zetas = [np.exp(-a)*r**2*xi for xi in xis]
    dzetas = [np.gradient(z, r) for z in zetas]

    n = n_basis
    A = np.zeros((n,n)); B = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = np.trapezoid(Pi*dzetas[i]*dzetas[j] + Q*xis[i]*xis[j], r)
            B[i,j] = np.trapezoid(S*xis[i]*xis[j], r)
    return A, B, xis, r

if __name__ == "__main__":
    bgd = build_metric(5e13)
    for n_basis in [1,2,3,4,5]:
        A, B, xis, r = build_RR_matrices(bgd, n_basis)
        asym = np.max(np.abs(A-A.T))
        bsym = np.max(np.abs(B-B.T))
        evals, evecs = eigh(A, B)
        sigma2 = evals[0]  # most negative / smallest eigenvalue
        print(f"n_basis={n_basis}  sigma^2={sigma2: .4f}  (A sym err={asym:.2e}, B sym err={bsym:.2e})")
    print()
    print("Corrected Paper C values (12000-pt grid): -27.85 -> -28.52 -> -28.55 -> -28.55 -> -28.55 across n=1..5")
