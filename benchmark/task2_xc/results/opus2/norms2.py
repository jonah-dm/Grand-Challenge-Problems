import numpy as np, sys
from scipy.optimize import fsolve
sys.path.insert(0, '.')
import submission as S
d = np.load('he_hf.npz')
rH = np.linspace(1e-5, 40, 200001); rhoH = np.exp(-2*rH)/np.pi; sigH = (2*rhoH)**2
def exH(k0, a1):
    S.KAPPA0 = k0
    e = 0.5*S._ex_unpol(2*rhoH, 4*sigH, 2*sigH/(8*rhoH), a1)
    return np.trapezoid(e*4*np.pi*rH**2, rH)
def exHe(k0, a1):
    S.KAPPA0 = k0
    r, rho, dr = d['r'], d['rho'], d['drho']; sig = dr**2
    e = S._ex_unpol(rho, sig, sig/(8*rho), a1)
    return np.trapezoid(e*4*np.pi*r**2, r)
f = lambda x: [exH(*x)+0.3125, exHe(*x)-float(d['ex'])]
x = fsolve(f, (0.174, 4.9)); print(x, f(x))
S.KAPPA0=0.174
from scipy.optimize import brentq
a=brentq(lambda a: exH(0.174,a)+0.3125, 1, 20); print('H-only a1', a, 'He err', exHe(0.174,a)-float(d['ex']))
for k0 in [0.12,0.15,0.174,0.2,0.25,0.3]:
    try:
        a=brentq(lambda a: exH(k0,a)+0.3125, 0.5, 60); print(k0, a, exHe(k0,a)-float(d['ex']))
    except Exception as e: print(k0, e)
