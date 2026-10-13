"""Derivation numerics: lone-electron branch satisfying (a) H norm, (b) H IP theorem / zero charge-slope
mismatch at q=1 (dE^SL/dq|_1 = -2J), (c) F'(s0)=0 (finite potential at the nucleus).
Form F = c0/(1 + c1 (s^2 - s0^2)^2). Then variational-stability test on trial orbitals."""
import numpy as np
from scipy.optimize import fsolve, minimize_scalar, minimize
r = np.linspace(1e-7, 40, 200001); dr = r[1]-r[0]
S0SQ = (6*np.pi)**(-2/3); a_hf = 0.25
CX = -(3/4)*(6/np.pi)**(1/3)
def Fform(s2, p):
    c0, c1 = p; return c0/(1 + c1*(s2-S0SQ)**2)
def EX(n, gn, p, F=Fform):
    s2 = gn**2/(4*(6*np.pi**2)**(2/3)*n**(8/3)+1e-300)
    return np.sum(4*np.pi*r**2*CX*n**(4/3)*F(s2, p))*dr
nH = np.exp(-2*r)/np.pi; gH = -2*nH
def conds(p):
    E1 = EX(nH, gH, p); h = 1e-4
    dEdq = (EX((1+h)*nH,(1+h)*gH,p) - EX((1-h)*nH,(1-h)*gH,p))/(2*h)
    return [E1 + 5/16, dEdq + 5/8]
p = fsolve(conds, [2.0, 5.0]); print("c0, c1 =", p, conds(p))
s = np.linspace(0, 3, 13); print("F(s):", np.round(Fform(s**2, p),3))
# charge line
for q in [0.25, 0.5, 0.75, 1.0]:
    print(q, "E_Hx(q) =", q*q*5/16 + EX(q*nH, q*gH, p))
def energy(phi, dphi, p=p, F=Fform):
    norm = np.sum(4*np.pi*r**2*phi**2)*dr; phi = phi/np.sqrt(norm); dphi = dphi/np.sqrt(norm)
    n = phi**2; gn = 2*phi*dphi
    T = 0.5*np.sum(4*np.pi*r**2*dphi**2)*dr; V = -np.sum(4*np.pi*r*n)*dr
    q = 4*np.pi*r**2*n*dr; vH = np.cumsum(q)/r + (np.sum(4*np.pi*r*n*dr) - np.cumsum(4*np.pi*r*n*dr))
    J = 0.5*np.sum(q*vH)
    return T+V+(1-a_hf)*(J+EX(n, gn, p, F))
# flexible radial variational minimization: phi = sum_k c_k exp(-z_k r)
zs = np.array([0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.0, 8.0])
E = np.exp(-np.outer(zs, r)); dE = -zs[:,None]*E
def fvar(c): return energy(c@E, c@dE)
c0 = np.zeros(len(zs)); c0[3] = 1.0
res = minimize(fvar, c0, method='BFGS', options={'maxiter':400})
print("H atom: E(exact density) =", fvar(c0), " variational min =", res.fun, " lowering (mEh) =", 1000*(res.fun-fvar(c0)))
def Fpbe(s2, p): return 1.804-0.804/(1+0.2195*s2/0.804)
res2 = minimize(lambda c: energy(c@E, c@dE, None, Fpbe), c0, method='BFGS', options={'maxiter':400})
print("PBE-x reference: E(exact)=", energy(c0@E, c0@dE, None, Fpbe), "min =", res2.fun, "lowering(mEh)=", 1000*(res2.fun-energy(c0@E, c0@dE, None, Fpbe)))
