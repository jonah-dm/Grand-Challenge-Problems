"""Linear blend of Minnesota-family hybrid meta-GGAs (libxc): MN15, M06-2X, PW6B95.
Each component's semilocal part is weighted; exact-exchange fraction is the weighted sum."""
from pyscf.dft import libxc
import numpy as np

COMP = {"mn15": 0.5, "m062x": 0.5, "pw6b95": 0.0}
HYB = float(sum(w * libxc.hybrid_coeff(n) for n, w in COMP.items()))
META = True


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    out = None
    for n, w in COMP.items():
        if w == 0:
            continue
        r = libxc.eval_xc(n, rho, spin, relativity, deriv, omega, verbose)
        out = [None if x is None else w * np.asarray(x) for x in r] if out is None else \
            [None if a is None else a + w * np.asarray(b) for a, b in zip(out, r)]
    return tuple(out)
