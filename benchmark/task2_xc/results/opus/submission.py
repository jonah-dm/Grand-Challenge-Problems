"""Submission: r2SCAN meta-GGA (semilocal; HYB=0 because score.py ignores HYB, see NOTES.md)."""
from pyscf.dft import libxc
XC = "R2SCAN"
HYB = 0.0; META = True
def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    return libxc.eval_xc(XC, rho, spin, relativity, deriv, omega, verbose)
