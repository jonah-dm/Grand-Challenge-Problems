"""Bathtub (rearrangement) bound on the SR exchange energy density of a one-orbital channel, and omega* (L4 with it)."""
import numpy as np, importlib.util, os
from scipy.optimize import brentq
from scipy.special import erfc, erf
spec = importlib.util.spec_from_file_location('o6', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../opus6c/submission.py'))
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
def G(y):   # int_0^y t erfc(t) dt  (closed form)
    return 0.5*y**2*erfc(y) - y*np.exp(-y**2)/(2*np.sqrt(np.pi)) + 0.25*erf(y)
def Fcap(n, w):
    """sup |e_x^SR|/|e_x^LDA,pol| over one-orbital channels with max density n (attained at the maximum point)"""
    R = (3/(4*np.pi*n))**(1/3)
    e = 0.5*n*4*np.pi*n*G(w*R)/w**2
    return e/(0.75*(6/np.pi)**(1/3)*n**(4/3))
def Fbr_H(w):
    s0 = (6*np.pi)**(-1/3)
    return 5/9*s0**-2*float(M.Qfun(np.array([1/w**2]))[0])/w**2
if __name__ == '__main__':
    from scipy.integrate import quad
    print('check G', G(1.3), quad(lambda t: t*erfc(t), 0, 1.3)[0])
    print('high-density limit', Fcap(1e6, 0.5), (4*np.pi**2/18)**(1/3))
    n = 1e-6; w = 0.5; print('low-density / on-top ratio', 0.5*n*4*np.pi*n*G(50)/w**2/(np.pi/(2*w**2)*n*n))
    ws = brentq(lambda w: Fbr_H(w) - Fcap(1/np.pi, w), 0.2, 3)
    print('omega* =', repr(ws), 'Fbr=Fcap=', Fbr_H(ws), Fcap(1/np.pi, ws))
    # check the bound holds for H at all r (b constant): branch F(r) vs cap F(n(r))
    for r in [0, 0.5, 1, 2, 4]:
        nr = np.exp(-2*r)/np.pi; k = (6*np.pi**2*nr)**(1/3)
        Fb = 5/9*(k/ws)**2*float(M.Qfun(np.array([1/ws**2]))[0])
        print(r, round(Fb, 4), round(Fcap(nr, ws), 4))
