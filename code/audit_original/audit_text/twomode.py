"""Two-mode nonlinear extension: r_phys = r*[lam + mu*(1-(r/R)^2)]. Reuses single-mode f_cal."""
import numpy as np
from homologous_calibrated import CalibratedHomologousModel
import background as bg
G, c = bg.G, bg.c

class TwoMode:
    def __init__(self, rho_c=5e13, f_cal=None):
        self.m = CalibratedHomologousModel(rho_c)
        self.f_cal = f_cal if f_cal is not None else self.m.calibrate_to_RR()[0]
        m = self.m
        self.r, self.Mr, self.R = m.r, m.Mr_bg, m.R
        self.x2 = (self.r/self.R)**2
        self.w = 4*np.pi*self.r**2*m.rho_bg        # dM per dr
        self.phi = [np.ones_like(self.r), 1-self.x2]
        self.dphi = [np.zeros_like(self.r), -2*self.r/self.R**2]   # r*phi' handled below
        r2 = self.r**2
        self.Mmat = np.array([[np.trapezoid(r2*self.phi[i]*self.phi[j]*self.w, self.r) for j in range(2)] for i in range(2)])
        self.R0 = np.zeros(2)
        self.R0 = self.Q_raw(np.array([1.0, 0.0]))

    def fields(self, q):
        lam, mu = q
        s = lam + mu*(1-self.x2)
        J = lam + mu*(1-3*self.x2)       # s + r s'
        return s, J

    def valid(self, q):
        lam, mu = q
        return (lam+mu > 0) and (lam-2*mu > 0)

    def Q_raw(self, q):
        r, Mr, w = self.r, self.Mr, self.w
        s, J = self.fields(q)
        rho = self.m.rho_bg/(s**2*J)
        P = np.where(rho > 1e-3, bg.P_of_rho(np.maximum(rho, 1e-3)), 0.0)
        rf = r*s
        mask = rf > 0
        gN = np.zeros_like(rf); gF = np.zeros_like(rf)
        gN[mask] = -G*Mr[mask]/rf[mask]**2
        num = G*(bg.eps_over_c2_of_rho(rho[mask])+P[mask]/c**2)/rho[mask]*(Mr[mask]+4*np.pi*rf[mask]**3*P[mask]/c**2)
        den = rf[mask]**2*(1-2*G*Mr[mask]/(rf[mask]*c**2))
        gF[mask] = -num/den
        geff = gN + self.f_cal*(gF-gN)
        Q = np.zeros(2)
        for k in range(2):
            phik = self.phi[k]; rdphik = r*self.dphi[k]/ (r+1e-300) * 0 + (0 if k == 0 else -2*self.x2)
            bracket = 2*phik/s + (phik + rdphik)/J
            Qint = np.trapezoid(np.where(rho > 0, P/rho, 0)*bracket*w, r)
            Qgrav = np.trapezoid(geff*r*phik*w, r)
            Q[k] = Qint + Qgrav
        return Q - self.R0

    def acc(self, q):
        return np.linalg.solve(self.Mmat, self.Q_raw(q))

    def jac(self, eps=1e-5):
        q0 = np.array([1.0, 0.0]); K = np.zeros((2, 2))
        for j in range(2):
            d = np.zeros(2); d[j] = eps
            K[:, j] = (self.acc(q0+d)-self.acc(q0-d))/(2*eps)
        return K

if __name__ == "__main__":
    t = TwoMode()
    print("f_cal", t.f_cal, "Mmat", t.Mmat)
    print("equilibrium acc", t.acc(np.array([1.0, 0.0])))
    for eps in (1e-4, 1e-5, 1e-6):
        K = t.jac(eps); w, v = np.linalg.eig(K)
        print(eps, "K=", K.tolist(), "eig(sigma^2)=", w, "\nvec", v)
    # single-mode check: mu=0 slice
    print("single-mode k:", t.m.linearized_k(t.f_cal), " two-mode K00:", t.jac()[0, 0])
