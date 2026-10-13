"""Derivation check: is the H-atom exactness of F_1 variationally stable? Radial trial orbitals
phi = e^{-r} + c g(r); energy E = <T+V> + (1-a)(J - Ex^{F1}...) i.e. E = <h> + J + a(-J) + (1-a)Ex^SL[n]."""
import numpy as np
from scipy.optimize import minimize_scalar
r = np.linspace(1e-7, 40, 400001); dr = r[1]-r[0]
A = 5/9; S0SQ = (6*np.pi)**(-2/3); a_hf = 0.25; m = 1.75
def energy(phi, dphi):
    norm = np.sum(4*np.pi*r**2*phi**2)*dr; phi = phi/np.sqrt(norm); dphi = dphi/np.sqrt(norm)
    n = phi**2; gn = 2*phi*dphi
    T = 0.5*np.sum(4*np.pi*r**2*dphi**2)*dr; V = -np.sum(4*np.pi*r*n)*dr
    q = 4*np.pi*r**2*n*dr; Qin = np.cumsum(q); vH = Qin/r + (np.sum(4*np.pi*r*n*dr) - np.cumsum(4*np.pi*r*n*dr))
    J = 0.5*np.sum(q*vH)
    s2 = gn**2/(4*(6*np.pi**2)**(2/3)*n**(8/3)+1e-300); hi = s2 >= S0SQ
    F = np.where(hi, A/np.where(hi, s2, 1), A/S0SQ*(s2/S0SQ)**(m/2))
    Ex = np.sum(4*np.pi*r**2*(-(3/4)*(6/np.pi)**(1/3)*n**(4/3))*F)*dr
    return T+V+J-a_hf*J+(1-a_hf)*Ex, T+V, J, Ex
for gname, g, dg in [("e^-2r", np.exp(-2*r), -2*np.exp(-2*r)), ("r e^-r", r*np.exp(-r), (1-r)*np.exp(-r)),
                     ("e^-r/2", np.exp(-r/2), -0.5*np.exp(-r/2)), ("r^2 e^-1.5r", r**2*np.exp(-1.5*r), (2*r-1.5*r**2)*np.exp(-1.5*r))]:
    f = lambda c: energy(np.exp(-r)+c*g, -np.exp(-r)+c*dg)[0]
    res = minimize_scalar(f, bounds=(-0.6, 0.6), method='bounded')
    print(f"{gname:12s} E(c=0)={f(0):.6f}  min E={res.fun:.6f} at c={res.x:+.4f}   lowering={1000*(res.fun-f(0)):.2f} mEh")
# scaled exponentials stay exact:
for z in [0.8, 1.2]:
    print("zeta",z, energy(np.exp(-z*r), -z*np.exp(-z*r)))
