import sys, numpy as np, importlib.util
from pyscf import gto, scf, dft
from diag2 import Fc
spec = importlib.util.spec_from_file_location("sub", "submission.py"); sub = importlib.util.module_from_spec(spec); spec.loader.exec_module(sub)
data = {}
for at, bas in [("He","aug-cc-pvqz"),("Ne","aug-cc-pvqz"),("Ar","aug-cc-pvqz")]:
    mol = gto.M(atom=f"{at} 0 0 0", basis=bas, verbose=0); mf = scf.RHF(mol).run()
    dm = mf.make_rdm1(); vk = mf.get_k(); exHF = -0.25*np.einsum('ij,ji', dm, vk)
    g = dft.gen_grid.Grids(mol); g.level = 5; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); r = dft.numint.eval_rho(mol, ao, dm/2, xctype='MGGA', with_lapl=False)
    data[at] = (exHF, r, g.weights)
def ex(at, F):
    exHF, r, w = data[at]; n = r[0]; sig = (r[1:4]**2).sum(0); tau = r[-1]; m = n > 1e-12
    n, sig, tau, w = n[m], sig[m], tau[m], w[m]
    z = np.minimum(sig/(8*n*tau), 1); h = z**3*(4-3*z); s2 = sig/(sub._CS2*n**(8/3))
    return 2*np.sum(w*sub._CX*n**(4/3)*((1-h)*sub._F0(s2) + h*F(s2)))
for at in data: print(at, "exact Ex =", round(data[at][0],5), " W0 (PBE-form) =", round(ex(at, sub._F0),5))
for c in [0.1895, 0.20, 0.21, 0.22, 0.23, 0.25, 0.2795]:
    F = Fc(c); print(c, {at: round(ex(at,F)-data[at][0],4) for at in data})
from scipy.optimize import brentq
cNe = brentq(lambda c: ex("Ne", Fc(c)) - data["Ne"][0], 0.22, 0.23, xtol=1e-6)
cAr = brentq(lambda c: ex("Ar", Fc(c)) - data["Ar"][0], 0.22, 0.23, xtol=1e-6)
F = Fc(cNe)
print("cNe =", cNe, " cAr =", cAr, " errors at cNe:", {at: round(ex(at,F)-data[at][0],5) for at in data})
