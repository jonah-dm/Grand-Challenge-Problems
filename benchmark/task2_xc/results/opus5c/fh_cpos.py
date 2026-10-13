import numpy as np
from scipy.optimize import brentq
from potexact_test import make_F
from fh_scan import Dh2p, Dh2, ex, S0SQ
s = np.exp(np.linspace(np.log(0.376), np.log(200), 400))
def minF(c, rmax=12.0):
    F0, FH = make_F(c, rmax); return FH(s**2).min()
for c in [0.17, 0.175, 0.18, 0.185, 0.19, 0.2]:
    F0, FH = make_F(c, 12.0); Fx = lambda s2: np.where(s2 < S0SQ, F0, FH(s2))
    a,b = ex(Dh2p, Fx); a2,b2 = ex(Dh2, Fx)
    print(f"c={c}: minF(s<200)={FH(s**2).min():.3f} at s={s[np.argmin(FH(s**2))]:.2f}; F(s0)={F0:.3f}; H2+ err={(b-a)*1000:.2f} mEh; H2chan err={(b2-a2)*1000:.2f} mEh")
