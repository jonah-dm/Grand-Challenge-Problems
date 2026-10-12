# Solve for (kappa1, mu1) of the iso-orbital (z=1) exchange limit from two exact norms:
#  E_x[H atom] = -5/16 (spin-scaled), E_x[He HF density] = exact HF exchange.
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from pyscf import gto, scf, dft
CX = -0.75 * (3 / np.pi) ** (1 / 3)
F = lambda s2, k, m: 1 + k - k / (1 + m * s2 / k)
def exH(k, m):
    def f(r):
        n2 = 2 * np.exp(-2 * r) / np.pi
        kf = (3 * np.pi**2 * n2) ** (1 / 3); s2 = (2 * n2) ** 2 / (4 * kf**2 * n2**2)
        return 0.5 * CX * n2 ** (4 / 3) * F(s2, k, m) * 4 * np.pi * r * r
    return quad(f, 1e-12, 60, limit=400)[0]
mol = gto.M(atom='He', basis='aug-cc-pv5z', verbose=0)
mf = scf.RHF(mol).run(); dm = mf.make_rdm1()
exx = -0.25 * np.einsum('ij,ji', dm, mf.get_k(dm))
g = dft.gen_grid.Grids(mol); g.level = 6; g.build()
ao = dft.numint.eval_ao(mol, g.coords, deriv=1)
r = dft.numint.eval_rho(mol, ao, dm, xctype='GGA')
m_ = r[0] > 1e-14; n = r[0][m_]; sig = (r[1:4]**2).sum(0)[m_]; w = g.weights[m_]
kf = (3*np.pi**2*n)**(1/3); s2 = sig/(4*kf**2*n**2)
exHe = lambda k, m: (CX * n**(4/3) * F(s2, k, m) * w).sum()
sol = fsolve(lambda p: [exH(*p) + 5/16, exHe(*p) - exx], [0.804, 0.27])
print('exact He Ex', exx, ' kappa1, mu1 =', sol, ' check', exH(*sol), exHe(*sol))
from scipy.optimize import brentq
muHe = brentq(lambda m: exHe(0.804, m) - exx, 0.1, 1.0)
print('He-normed mu1 (kappa=0.804):', muHe, ' resulting E_x[H] =', exH(0.804, muHe))
