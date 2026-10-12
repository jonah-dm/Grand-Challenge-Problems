from pyscf import gto, dft
from pyscf.dft import libxc
mol = gto.M(atom="H 0 0 0; H 0 0 0.741", basis="def2-svp", verbose=0)
def mk(semi):
    def ev(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
        return libxc.eval_xc(semi, rho, spin, relativity, deriv, omega, verbose)
    return ev
for semi, hyb in [("GGA_X_PBE + GGA_C_PBE", 0.0), ("0.8*GGA_X_PBE + GGA_C_PBE", 0.0), ("0.8*GGA_X_PBE + GGA_C_PBE", 0.2)]:
    mf = dft.RKS(mol); mf.grids.level = 3
    mf = mf.define_xc_(mk(semi), "GGA", hyb=hyb)
    print(semi, hyb, mf.kernel())
