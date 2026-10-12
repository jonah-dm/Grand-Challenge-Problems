"""Submission wrapper: global meta-GGA hybrid built from libxc pieces (see NOTES.md)."""
from pyscf.dft import libxc

XC = "PW6B95"
HYB = float(libxc.hybrid_coeff(XC))
META = libxc.is_meta_gga(XC)


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    # libxc.eval_xc returns only the semilocal (DFT) part; the exact-exchange
    # fraction HYB is added by PySCF through define_xc_(..., hyb=HYB).
    return libxc.eval_xc(XC, rho, spin, relativity, deriv, omega, verbose)
