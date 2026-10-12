import numpy as np, sys
from scipy.optimize import fsolve
sys.path.insert(0, '.')
import submission as S
d = np.load('he_hf.npz')
rH = np.linspace(1e-5, 40, 200001); rhoH = np.exp(-2*rH)/np.pi; sigH = (2*rhoH)**2
def exH(mu0, a1):
    S.MU0 = mu0
    e = 0.5*S._ex_unpol(2*rhoH, 4*sigH, 2*sigH/(8*rhoH), a1)
    return np.trapezoid(e*4*np.pi*rH**2, rH)
def exHe(mu0, a1):
    S.MU0 = mu0
    r, rho, dr = d['r'], d['rho'], d['drho']; sig = dr**2
    e = S._ex_unpol(rho, sig, sig/(8*rho), a1)
    return np.trapezoid(e*4*np.pi*r**2, r)
def f(x):
    return [exH(*x)+0.3125, exHe(*x)-float(d['ex'])]
for g in [(0.3,5.),(1.,5.),(0.5,8.),(2.,3.)]:
    x, info, ier, msg = fsolve(f, g, full_output=True)
    print(g, x, f(x), ier)
