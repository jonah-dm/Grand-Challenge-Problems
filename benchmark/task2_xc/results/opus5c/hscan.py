import numpy as np, sys
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from pyscf import gto, dft
from hbranch import build, S0SQ, CX
E6 = [35.52322122,6.513143725,1.822142904,0.6259552659,0.2430767471,0.1001124280]
C6 = [0.00916359628,0.04936149294,0.16853830490,0.37056279970,0.41649152980,0.13033408410]
def basis(zeta):
    sc = (zeta/1.24)**2
    return gto.basis.parse("H S\n" + "\n".join(f"  {e*sc} {c}" for e,c in zip(E6,C6)))
def dens(atom, zeta, spin, charge):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=spin, charge=charge, verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, mol.intor('int2e'))
    g = dft.gen_grid.Grids(mol); g.level = 7; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); m = n > 1e-14
    return J, n[m], sig[m], g.weights[m]
def Ffun(t, FF):
    cs = CubicSpline(t, FF, bc_type=((1, 0.0), (1, 0.0)))
    return lambda s2: np.where(np.log(s2/S0SQ) <= 0, FF[0], np.where(np.log(s2/S0SQ) >= t[-1], FF[-1], cs(np.clip(np.log(s2/S0SQ), 0, t[-1]))))
def ex(D, F):
    J, n, sig, w = D; s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3))
    return np.sum(w*CX*n**(4/3)*F(s2)) + J
r = np.linspace(1e-6, 30, 60001); dr = r[1]-r[0]
def energyH(phi, dphi, F, a_hf=0.25):
    norm = np.sum(4*np.pi*r**2*phi**2)*dr; phi = phi/np.sqrt(norm); dphi = dphi/np.sqrt(norm)
    n = phi**2; gn = 2*phi*dphi
    T = 0.5*np.sum(4*np.pi*r**2*dphi**2)*dr; V = -np.sum(4*np.pi*r*n)*dr
    q = 4*np.pi*r**2*n*dr; vH = np.cumsum(q)/r + (np.sum(4*np.pi*r*n*dr) - np.cumsum(4*np.pi*r*n*dr)); J = 0.5*np.sum(q*vH)
    s2 = gn**2/(4*(6*np.pi**2)**(2/3)*n**(8/3)+1e-300)
    return T+V+(1-a_hf)*(J+np.sum(4*np.pi*r**2*CX*n**(4/3)*F(s2))*dr)
if __name__ == "__main__":
    Dh2p = dens("H 0 0 -1; H 0 0 1", 1.238, 1, 1); Dh2 = dens("H 0 0 -0.7; H 0 0 0.7", 1.166, 0, 0)
    zs = np.array([0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.0, 8.0]); EE = np.exp(-np.outer(zs, r)); dEE = -zs[:,None]*EE
    rmax = float(sys.argv[1])
    for c in [float(x) for x in sys.argv[2:]]:
        F0, t, FF = build(c, rmax); F = Ffun(t, FF)
        c0 = np.zeros(len(zs)); c0[3] = 1.0; f = lambda cc: energyH(cc@EE, cc@dEE, F)
        res = minimize(f, c0, method='BFGS', options={'maxiter':200})
        print(f"rmax={rmax} c={c}: F(s0)={F0:.4f} Fend={FF[-1]:.3f} at s={np.sqrt(S0SQ*np.exp(t[-1])):.2f}; H2+ err={ex(Dh2p,F)*1e3:+.2f} H2ch err={ex(Dh2,F)*1e3:+.2f} mEh; H E={f(c0):.6f} relaxed={res.fun:.6f}", flush=True)
