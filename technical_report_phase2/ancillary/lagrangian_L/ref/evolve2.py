import numpy as np, pickle
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
import initial_profile as IP, eos as E

d = pickle.load(open('initdata2.pkl','rb'))
r_lab, m0, Gamma_shell, L = d['r_lab'], d['m_units'], d['Gamma_shell'], d['L']
N = len(r_lab)
Gscale = E.G*L**2/E.c**2

def p_dimless_of_rho_dimless(rho_d):
    rho_SI = np.maximum(rho_d,0)/Gscale
    P_SI = IP.P_of_rho(rho_SI)
    return Gscale*(P_SI/E.c**2)

def derivs(y_of_r, r):
    cs = CubicSpline(r, y_of_r)
    return cs(r, 1)   # first derivative, evaluated smoothly via the spline

def rho_from_state(R, m):
    Rp = derivs(R, r_lab)
    mp = derivs(m, r_lab)
    rho = mp/(4*np.pi*np.maximum(R,1e-12)**2*np.where(np.abs(Rp)>1e-12,Rp,1e-12))
    return np.maximum(rho,0), Rp

def lapse_from_state(rho_d, p_d):
    dp_dr = derivs(p_d, r_lab)
    denom = rho_d + p_d
    threshold = 1e-6*np.max(rho_d)
    denom_safe = np.maximum(denom, 1e-300)
    weight = denom**2 / (denom**2 + threshold**2)
    integrand = weight * (-dp_dr/denom_safe)
    ln_a_rel = cumulative_trapezoid(integrand[::-1], r_lab[::-1], initial=0.0)[::-1]
    return np.clip(np.exp(ln_a_rel), 1e-6, 1e2)

def rhs(t, y):
    R, m, U = y[0:N], y[N:2*N], y[2*N:3*N]
    rho, Rp = rho_from_state(R, m)
    p_d = p_dimless_of_rho_dimless(rho)
    a = lapse_from_state(rho, p_d)
    X = Rp/np.maximum(Gamma_shell,1e-12)
    dp_dr = derivs(p_d, r_lab)
    # Near the free surface, rho+p -> 0 and the pressure-gradient force is physically
    # finite but numerically ill-conditioned (0/0). A HARD threshold mask here is worse
    # than the disease: it flips discontinuously between RK4 sub-stages (which sample
    # slightly different intermediate states), injecting a spurious one-step velocity
    # kick exactly where a shell's density crosses the cutoff. Use a smooth weight
    # instead so the RHS stays continuous.
    denom = rho + p_d
    threshold = 1e-6 * np.max(rho)
    denom_safe = np.maximum(denom, 1e-300)
    weight = denom**2 / (denom**2 + threshold**2)   # -> 1 for denom >> threshold, -> 0 smoothly for denom << threshold
    pgrad = weight * (Gamma_shell/denom_safe) * (dp_dr/np.where(np.abs(X) > 1e-12, X, 1e-12))
    dRdt = a*U
    dmdt = -4*np.pi*a*p_d*R**2*U
    dUdt = -a**2*((m+4*np.pi*R**3*p_d)/np.maximum(R,1e-12)**2 + pgrad)
    return np.concatenate([dRdt, dmdt, dUdt]), rho, p_d, a

y0 = np.concatenate([r_lab.copy(), m0.copy(), np.zeros(N)])
dy, rho0, p0, a0 = rhs(0.0, y0)
print("max|dR/dt| =", np.max(np.abs(dy[0:N])))
print("max|dU/dt| =", np.max(np.abs(dy[2*N:3*N])), "  (expect O(m/R^2) ~ 1e-3, Newtonian-scale)")
print("rho_d[0] (center) =", rho0[0], " vs expected ~9.5e-3")
print("a range:", a0.min(), a0.max())
print("Any NaN?", np.any(np.isnan(dy)))
i = np.argmax(np.abs(dy[2*N:3*N]))
print(f"worst shell i={i}, r_lab={r_lab[i]:.4f}, rho_d={rho0[i]:.3e}, p_d={p0[i]:.3e}, a={a0[i]:.4f}")
