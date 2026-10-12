from pyscf.dft import libxc
XC = "M062X"
HYB = float(libxc.hybrid_coeff(XC)); META = True
def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    return libxc.eval_xc(XC, rho, spin, relativity, deriv, omega, verbose)
