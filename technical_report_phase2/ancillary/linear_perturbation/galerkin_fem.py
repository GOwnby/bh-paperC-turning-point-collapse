"""
Linear finite-element (Galerkin) assembly for the Cowling-approximation radial
pulsation problem, built from the symbolically-derived and verified symmetric
weak form:
  integrand = A(r)*xi*eta + C(r)*(xi*eta' + xi'*eta) + D(r)*xi'*eta'
  A = 4*Gamma1*P0 + P0'*r^2*rho0'/rho0 + 4*P0'*r
  C = r*(2*Gamma1*P0 + P0'*r)
  D = Gamma1*P0*r^2
Mass (kinetic) term: rho0*r^2*xi*eta
Symmetric by construction -- this is the whole point of building from the weak form.
"""
import numpy as np

def assemble(r, rho0, P0, Gamma1, G=None, include_self_gravity=False):
    N = len(r)
    rho0p = np.gradient(rho0, r)
    P0p = np.gradient(P0, r)
    A = 4*Gamma1*P0 + P0p*r**2*rho0p/np.maximum(rho0,1e-300) + 4*P0p*r
    if include_self_gravity:
        # Exact closed form for l=0 (radial) perturbations: phi' = -4 pi G rho0 xi,
        # verified symbolically (verify_exact_poisson.py) to solve the perturbed
        # Poisson equation exactly -- no coupled field or separate solve needed.
        # This adds -4 pi G rho0^2 r^2 to the xi*eta coefficient (no derivatives,
        # so manifestly symmetric).
        assert G is not None, "G required when include_self_gravity=True"
        A = A - 4*np.pi*G*rho0**2*r**2
    C = r*(2*Gamma1*P0 + P0p*r)
    D = Gamma1*P0*r**2

    M = np.zeros((N,N))
    Dmass = np.zeros((N,N))
    # 2-point Gauss quadrature per element for accuracy
    gauss_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]
    for k in range(N-1):
        rk, rk1 = r[k], r[k+1]
        h = rk1 - rk
        # local matrices for element k (nodes k, k+1), linear basis
        Mloc = np.zeros((2,2)); Dloc = np.zeros((2,2))
        for xi_g in gauss_pts:
            rr = 0.5*(rk1-rk)*xi_g + 0.5*(rk1+rk)
            w = h/2  # weight for 2-pt Gauss on [-1,1] mapped, each point weight 1 in reference -> w=h/2 total /2pts... use standard: weight=1 each, jacobian=h/2
            # basis values and derivatives at rr
            phi = [(rk1-rr)/h, (rr-rk)/h]
            dphi = [-1/h, 1/h]
            # interpolate A, C, D at rr (linear interpolation from node values)
            tA = np.interp(rr, [rk,rk1], [A[k],A[k+1]])
            tC = np.interp(rr, [rk,rk1], [C[k],C[k+1]])
            tD = np.interp(rr, [rk,rk1], [D[k],D[k+1]])
            trho = np.interp(rr, [rk,rk1], [rho0[k],rho0[k+1]])
            for i in range(2):
                for j in range(2):
                    Mloc[i,j] += (h/2)*1.0*( tA*phi[i]*phi[j] + tC*(phi[i]*dphi[j]+dphi[i]*phi[j]) + tD*dphi[i]*dphi[j] )
                    Dloc[i,j] += (h/2)*1.0*( trho*rr**2*phi[i]*phi[j] )
        M[k:k+2, k:k+2] += Mloc
        Dmass[k:k+2, k:k+2] += Dloc

    # CORRECTED surface boundary term: the integration by parts of the d(P')/dr
    # piece assumed [eta*r^2*P']_surface = 0, which is WRONG -- the physical free-
    # surface condition is Delta_P(R)=0, not P'(R)=0, and these differ by xi*P0':
    # P'(R) = -xi(R)*P0'(R) (verified symbolically). This contributes a genuine
    # boundary term -R^2*P0'(R)*xi(R)*eta(R) to the weak form, which was missing.
    M[-1,-1] += -r[-1]**2 * P0p[-1]
    return M, Dmass

print("galerkin_fem module loaded OK")
