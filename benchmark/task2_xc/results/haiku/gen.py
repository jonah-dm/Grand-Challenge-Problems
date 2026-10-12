import sys
a, name = float(sys.argv[1]), sys.argv[2]
src = f'''"""Hybrid GGA: (1-a)*PBE_x + PBE_c + a*HF exchange (PBE0-type family)."""
from pyscf.dft import libxc
HYB = {a}
XC_SEMI = "{1-a}*GGA_X_PBE + GGA_C_PBE"
META = False

def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    return libxc.eval_xc(XC_SEMI, rho, spin, relativity, deriv, omega, verbose)
'''
open(f"scan/{name}.py","w").write(src)
