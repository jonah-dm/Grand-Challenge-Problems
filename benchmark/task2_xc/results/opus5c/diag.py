import os, sys, numpy as np, importlib.util
from pyscf import gto, dft
sys.path.insert(0,'../../kit'); from geom import to_str; import dev_set as S
spec = importlib.util.spec_from_file_location("sub", "submission.py"); sub = importlib.util.module_from_spec(spec); spec.loader.exec_module(sub)
def dens(name):
    atoms, ch, sp = S.SPECIES[name]
    mol = gto.M(atom=to_str(atoms), basis='def2-svp', charge=ch, spin=sp, verbose=0)
    mf = dft.UKS(mol); mf.xc='pbe0'; mf.grids.level=3; mf.kernel()
    ao = dft.numint.eval_ao(mol, mf.grids.coords, deriv=1)
    dm = mf.make_rdm1()
    ra = dft.numint.eval_rho(mol, ao, dm[0], xctype='MGGA', with_lapl=False)
    rb = dft.numint.eval_rho(mol, ao, dm[1], xctype='MGGA', with_lapl=False)
    return ra, rb, mf.grids.weights
def dE(name):
    ra, rb, w = dens(name)
    na, nb = ra[0], rb[0]; out = {}
    for sch in [(ra,), (rb,)]:
        pass
    # per-channel exchange difference F_H - F0 weighted by L
    tot = np.zeros(4); bins = [(0,0.4),(0.4,1),(1,3),(3,1e9)]
    for r in (ra, rb):
        n = r[0]; sig = (r[1:4]**2).sum(0); tau = r[-1]; m = n > 1e-10
        n, sig, tau, ww = n[m], sig[m], tau[m], w[m]
        z = np.minimum(sig/(8*n*tau), 1); h = z**3*(4-3*z)
        s2 = sig/(sub._CS2*n**(8/3))
        d = 0.75*sub._CX*n**(4/3)*h*(sub._FH(s2+0j).real - sub._F0(s2))
        s = np.sqrt(s2)
        for i,(a,b) in enumerate(bins): tot[i] += np.sum(ww*d*((s>=a)&(s<b)))
    return tot
for mol, atoms in [("N2", {"N_atom":2}), ("CH4", {"C_atom":1,"H_atom":4})]:
    dm = dE(mol); da = sum(c*dE(a) for a,c in atoms.items())
    print(mol, "Delta(W2-W0) xc energy by s-bin [0,.4),[.4,1),[1,3),[3,inf): molecule", np.round(dm*627.5,1), "atoms", np.round(da*627.5,1), "AE change", np.round((da-dm)*627.5,1))
