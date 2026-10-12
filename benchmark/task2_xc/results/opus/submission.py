"""Submission: equal-weight ensemble of two high-exact-exchange Minnesota global
meta-GGA hybrids, E_xc = 0.5*E_xc[M06-2X] + 0.5*E_xc[MN15]  (exact exchange 49%).
See NOTES.md for rationale."""
from pyscf.dft import libxc

XC = "0.5*M062X+0.5*MN15"
HYB = float(libxc.hybrid_coeff(XC))   # 0.49 global exact exchange
META = True


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    # Semilocal part only; PySCF adds HYB * exact exchange via define_xc_.
    return libxc.eval_xc(XC, rho, spin, relativity, deriv, omega, verbose)
