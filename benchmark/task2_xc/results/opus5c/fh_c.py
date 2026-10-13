import numpy as np
from scipy.optimize import brentq
from potexact import shoot, growth, S0
from fh_scan import Dh2p, Dh2, ex, S0SQ
from potexact_test import make_F
def h2p_err(c):
    F0, FH = make_F(c); Fx = lambda s2: np.where(s2 < S0SQ, F0, FH(s2))
    a, b = ex(Dh2p, Fx); return b - a
cs = brentq(h2p_err, 0.12, 0.2, xtol=1e-5); print("c from H2+ norm:", cs)
F0, FH = make_F(cs, rmax=9.0)
s = np.array([0.376, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 20])
print("s:", s); print("F:", np.round(FH(s**2), 4))
Fx = lambda s2: np.where(s2 < S0SQ, F0, FH(s2)); print("H2 chan:", ex(Dh2, Fx))
