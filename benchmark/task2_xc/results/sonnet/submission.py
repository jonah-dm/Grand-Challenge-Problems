"""Linear blend of Minnesota-family hybrid meta-GGAs (libxc): MN15, M06-2X, PW6B95.
Each component's semilocal part is weighted; exact-exchange fraction is the weighted sum."""
from pyscf.dft import libxc
import numpy as np

COMP = {"mn15": 0.4, "m062x": 0.4, "pw6b95": 0.2}
HYB = float(sum(w * libxc.hybrid_coeff(n) for n, w in COMP.items()))
META = True


def _comb(a, b, w):
    if b is None:
        return a
    if isinstance(b, (tuple, list)):
        return tuple(_comb(None if a is None else a[k], b[k], w) for k in range(len(b)))
    b = w * np.asarray(b)
    return b if a is None else a + b


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    out = None
    rho = np.asarray(rho)
    if rho.shape[-2] == 5:  # (rho, grad, tau) -> insert zero laplacian row for libxc
        z = np.zeros_like(rho[..., :1, :])
        rho = np.concatenate([rho[..., :4, :], z, rho[..., 4:, :]], axis=-2)
    for n, w in COMP.items():
        r = libxc.eval_xc(n, rho, spin, relativity, deriv, omega, verbose)
        out = tuple(_comb(None if out is None else out[k], r[k], w) for k in range(len(r)))
    v = list(out[1])
    if np.ndim(v[2]) == 0:
        v[2] = np.zeros_like(v[0])
    return out[0], tuple(v), out[2], out[3]
