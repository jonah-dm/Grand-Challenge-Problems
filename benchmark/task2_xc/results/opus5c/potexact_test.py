import numpy as np
from scipy.optimize import brentq, minimize
from scipy.interpolate import interp1d
from potexact import shoot, growth, S0
import ipth
from ipth import energy, E, dE, zs
def make_F(c, rmax=8.0):
    F0 = brentq(lambda f: growth(f, c), -5, 10)
    sol = shoot(F0, c, rmax)
    rr = np.linspace(1e-4, rmax, 4000); FF = sol.sol(rr)[0]
    ss2 = (S0*np.exp(2*rr/3))**2
    fint = interp1d(np.log(ss2), FF, bounds_error=False, fill_value=(FF[0], FF[-1]))
    return F0, (lambda s2, p=None: fint(np.log(np.maximum(s2, 1e-30))))
# find c* where F''(s0)=0 (monotone)
for c in [0.25, 0.2795, 0.32]:
    F0, Ff = make_F(c)
    s = np.array([0.2, 0.376, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0])
    c0 = np.zeros(len(zs)); c0[3] = 1.0
    e_exact = energy(c0@E, c0@dE, None, Ff)
    res = minimize(lambda cc: energy(cc@E, cc@dE, None, Ff), c0, method='BFGS', options={'maxiter':300})
    print(f"c={c}: F(s0)={F0:.3f} F(s)={np.round(Ff(s**2),3)}  E(exact dens)={e_exact:.6f}  var.min={res.fun:.6f}  lowering={1000*(res.fun-e_exact):.3f} mEh")
