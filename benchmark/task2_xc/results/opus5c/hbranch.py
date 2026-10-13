"""One-orbital exchange branch F_H(s): GGA enhancement whose exchange potential on the H-atom density
equals the exact exchange potential v_x = -v_H up to a constant c (derivation in NOTES.md, Sec. 6).
ODE (s = s0 e^{2r/3}):  (3/4) F'' + F' (3/(2r) - 5/2) + (4/3) F = e^{2r/3} (v_x(r) + c) / (C pi^{-1/3}).
Regular at r = 0 (F'(0) = 0); F(0) chosen to suppress the e^{8r/3} (s^4) mode at r = rmax.
F is held constant below s0 and beyond its first minimum (both points have F' = 0, so C1)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
CX = -(3/4)*(6/np.pi)**(1/3); K = CX*np.pi**(-1/3)
S0 = (6*np.pi)**(-1/3); S0SQ = S0**2
def vx(r): return -(1 - np.exp(-2*r)*(1+r))/r
def _rhs(r, y, c):
    F, Fp = y
    return [Fp, (np.exp(2*r/3)*(vx(r)+c)/K - Fp*(1.5/r - 2.5) - (4/3)*F)/0.75]
def shoot(F0, c, rmax):
    r0 = 1e-5; a = ((vx(r0)+c)/K - (4/3)*F0)/4.5
    return solve_ivp(_rhs, [r0, rmax], [F0 + a*r0**2, 2*a*r0], args=(c,), rtol=1e-11, atol=1e-13, dense_output=True)
def build(c, rmax=6.0, npts=2000):
    g = lambda f: shoot(f, c, rmax).y[0, -1]
    F0 = brentq(g, -5, 10, xtol=1e-13)
    sol = shoot(F0, c, rmax)
    rr = np.linspace(1e-5, rmax, npts); FF = sol.sol(rr)[0]; dF = sol.sol(rr)[1]
    idx = np.where((dF[1:] > 0) & (dF[:-1] <= 0) & (rr[1:] > 0.3))[0]
    k = idx[0]+1 if len(idx) else npts-1
    return F0, 4*rr[:k+1]/3, FF[:k+1]          # t = ln(s^2/s0^2) = 4r/3
