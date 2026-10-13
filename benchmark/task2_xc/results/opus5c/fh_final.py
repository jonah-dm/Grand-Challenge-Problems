import numpy as np
from scipy.optimize import brentq, minimize
from scipy.interpolate import CubicSpline
from potexact import shoot, growth, S0
from fh_scan import Dh2p, Dh2, ex, S0SQ
from ipth import energy, E, dE, zs
def build(c, rmax=6.0):
    F0 = brentq(lambda f: growth(f, c, rmax), -5, 10)
    sol = shoot(F0, c, rmax)
    rr = np.linspace(1e-4, rmax, 3000); FF = sol.sol(rr)[0]
    # truncate at first interior minimum (F' = 0) beyond the core region; constant continuation
    dF = np.gradient(FF, rr); idx = np.where((dF[1:] > 0) & (dF[:-1] <= 0) & (rr[1:] > 0.3))[0]
    k = idx[0]+1 if len(idx) else len(rr)-1
    t = 2*rr[:k+1]/3*2          # t = ln(s^2/s0^2) = 4r/3
    return F0, t, FF[:k+1], rr[k]
def Ffun(F0, t, FF):
    cs = CubicSpline(t, FF, bc_type=((1, 0.0), (1, 0.0)))
    def F(s2):
        tt = np.log(np.maximum(s2, 1e-300)/S0SQ)
        return np.where(tt <= 0, FF[0], np.where(tt >= t[-1], FF[-1], cs(np.clip(tt, 0, t[-1]))))
    return F
for rmax in [6.0, 7.0]:
  for c in [0.16, 0.17, 0.18, 0.19, 0.20, 0.21]:
    try:
        F0, t, FF, rk = build(c, rmax)
    except Exception as e:
        print(c, "fail", e); continue
    F = Ffun(F0, t, FF)
    a,b = ex(Dh2p, F); a2,b2 = ex(Dh2, F)
    c0 = np.zeros(len(zs)); c0[3] = 1.0
    e0 = energy(c0@E, c0@dE, None, lambda s2,p=None: F(s2))
    res = minimize(lambda cc: energy(cc@E, cc@dE, None, lambda s2,p=None: F(s2)), c0, method='BFGS', options={'maxiter':300})
    print(f"rmax={rmax} c={c}: F(s0)={F0:.3f} Fmin={FF[-1]:.3f} at s={S0*np.exp(2*rk/3):.2f}; H2+ err={(b-a)*1e3:+.2f} H2ch err={(b2-a2)*1e3:+.2f} mEh; H: E(exact)={e0:.6f} SCF-like min={res.fun:.6f}")
