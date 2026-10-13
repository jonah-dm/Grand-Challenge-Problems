import sys, numpy as np, importlib.util
from scipy.interpolate import CubicSpline
from pyscf import gto, dft
sys.path.insert(0,'../../kit'); from geom import to_str; import dev_set as S
from hbranch import shoot, S0SQ
from scipy.optimize import brentq
spec = importlib.util.spec_from_file_location("sub", "submission.py"); sub = importlib.util.module_from_spec(spec); spec.loader.exec_module(sub)
def Fc(c, rmax=6.0):
    F0 = brentq(lambda f: shoot(f, c, rmax).y[0,-1], -5, 10, xtol=1e-13); sol = shoot(F0, c, rmax)
    rr = np.linspace(1e-5, rmax, 2000); FF = sol.sol(rr)[0]; t = 4*rr/3
    cs = CubicSpline(t, FF)
    return lambda s2: np.where(np.log(s2/S0SQ)<=0, FF[0], np.where(np.log(s2/S0SQ)>=t[-1], FF[-1], cs(np.clip(np.log(s2/S0SQ),0,t[-1]))))
cache = {}
def dens(name):
    if name in cache: return cache[name]
    atoms, ch, sp = S.SPECIES[name]
    mol = gto.M(atom=to_str(atoms), basis='def2-svp', charge=ch, spin=sp, verbose=0)
    mf = dft.UKS(mol); mf.xc='pbe0'; mf.grids.level=3; mf.kernel()
    ao = dft.numint.eval_ao(mol, mf.grids.coords, deriv=1); dm = mf.make_rdm1()
    out = [dft.numint.eval_rho(mol, ao, dm[i], xctype='MGGA', with_lapl=False) for i in (0,1)], mf.grids.weights
    cache[name] = out; return out
def dE(name, F):
    (ra, rb), w = dens(name); tot = 0
    for r in (ra, rb):
        n = r[0]; sig = (r[1:4]**2).sum(0); tau = r[-1]; m = n > 1e-10
        n, sig, tau, ww = n[m], sig[m], tau[m], w[m]
        z = np.minimum(sig/(8*n*tau), 1); h = z**3*(4-3*z); s2 = sig/(sub._CS2*n**(8/3))
        tot += np.sum(ww*0.75*sub._CX*n**(4/3)*h*(F(s2) - sub._F0(s2)))
    return tot
tests = [("AE:N2", {"N_atom":2,"N2":-1}), ("AE:CH4", {"C_atom":1,"H_atom":4,"CH4":-1}), ("AE:H2O", {"O_atom":1,"H_atom":2,"H2O":-1}),
         ("AE:HF", {"F_atom":1,"H_atom":1,"HF":-1}), ("AE:C2H2", {"C_atom":2,"H_atom":2,"C2H2":-1}), ("IP:B", {"B_cation":1,"B_atom":-1}), ("IP:F", {"F_cation":1,"F_atom":-1})]
for c in [float(x) for x in sys.argv[1:]]:
    F = Fc(c)
    print(f"c={c}: F(s)= ", np.round(F(np.array([0.376,0.7,1,1.5,2,3,5,10])**2),2))
    print("   non-SCF change vs W0 (kcal/mol):", {k: round(627.5*sum(v*dE(s,F) for s,v in e.items()),1) for k,e in tests})
