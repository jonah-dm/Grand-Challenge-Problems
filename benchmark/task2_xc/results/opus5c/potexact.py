"""Derivation numerics: GGA iso-orbital enhancement F(s) whose exchange potential on the H density
equals v_x^exact(r) + c.  ODE in r (s = s0 e^{2r/3}):
 (3/4) F'' + F' (3/(2r) - 5/2) + (4/3) F = e^{2r/3} (v_x^ex(r) + c) / (C pi^{-1/3})
Regular solution at r=0 (F'(0)=0); shoot on F(0) to kill the s^4 (e^{8r/3}) mode."""
import numpy as np
from scipy.integrate import solve_ivp
C = -(3/4)*(6/np.pi)**(1/3); K = C*np.pi**(-1/3)
S0 = (6*np.pi)**(-1/3)
def vx(r): return -(1 - np.exp(-2*r)*(1+r))/r
def rhs(r, y, c):
    F, Fp = y
    src = np.exp(2*r/3)*(vx(r)+c)/K
    Fpp = (src - Fp*(1.5/r - 2.5) - (4/3)*F)/0.75
    return [Fp, Fpp]
def shoot(F0, c, rmax=9.0):
    r0 = 1e-4
    # series: near 0, F ~ F0 + a r^2 ; plug: (3/4)2a + 2a*(3/2) ... -> (1.5a + 3a) + (4/3)F0 = src(0)
    src0 = (vx(r0)+c)/K
    a = (src0 - (4/3)*F0)/4.5
    sol = solve_ivp(rhs, [r0, rmax], [F0 + a*r0**2, 2*a*r0], args=(c,), rtol=1e-10, atol=1e-12, dense_output=True)
    return sol
def growth(F0, c, rmax=9.0):
    sol = shoot(F0, c, rmax); return sol.y[0,-1]*np.exp(-8*rmax/3)
from scipy.optimize import brentq
for c in [0.0, 0.1, 0.2, 0.3]:
    Fs = np.linspace(-5, 10, 61); g = [growth(f, c) for f in Fs]
    roots = [brentq(lambda f: growth(f,c), Fs[i], Fs[i+1]) for i in range(len(Fs)-1) if np.sign(g[i])!=np.sign(g[i+1])]
    for F0 in roots:
        sol = shoot(F0, c, 8.0)
        rr = np.array([0.0001, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8])
        F = sol.sol(rr)[0]
        print(f"c={c}: F(s0)={F0:.4f}; s=", np.round(S0*np.exp(2*rr/3),3), "F=", np.round(F,3))
