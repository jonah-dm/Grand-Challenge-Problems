"""Pure GGA candidate: Becke-88 exchange + LYP correlation (BLYP-type). No exact exchange."""
from pyscf.dft import libxc

HYB = 0.0
META = False
XC_SEMI = "GGA_X_B88 + GGA_C_LYP"


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    return libxc.eval_xc(XC_SEMI, rho, spin, relativity, deriv, omega, verbose)
