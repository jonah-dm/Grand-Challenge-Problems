import numpy as np, sys
sys.path.insert(0,'.')
from pyscf import gto, dft
import cand_A as SA, cand_C as SC
mol = gto.M(atom='N 0 0 0', basis='def2-svp', spin=3, verbose=0)
mf = dft.UKS(mol); mf.grids.level=4; mf = mf.define_xc_(SA.eval_xc, 'MGGA', hyb=0.25); mf.xc='HF'
mf = dft.UKS(mol); mf.grids.level=4; mf.xc='HF'; mf = mf.define_xc_(lambda *a, **k: SA.eval_xc(*a, **k), 'MGGA', hyb=0.25); mf.kernel()
g = mf.grids; ao = dft.numint.eval_ao(mol, g.coords, deriv=1)
dm = mf.make_rdm1()
ra = dft.numint.eval_rho(mol, ao, dm[0], xctype='MGGA'); rb = dft.numint.eval_rho(mol, ao, dm[1], xctype='MGGA')
na, nb = ra[0], rb[0]; gaa=(ra[1:4]**2).sum(0); gbb=(rb[1:4]**2).sum(0); gab=(ra[1:4]*rb[1:4]).sum(0)
ta, tb = ra[-1], rb[-1]
eA = SA.edens(na,nb,gaa,gab,gbb,ta,tb,mode='A'); eC = SC.edens(na,nb,gaa,gab,gbb,ta,tb,mode='C')
r = np.linalg.norm(g.coords,axis=1); n=na+nb; zeta=(na-nb)/np.maximum(n,1e-30)
tw=(gaa+2*gab+gbb)/(8*np.maximum(n,1e-30)); zt=np.minimum(tw/np.maximum(ta+tb,1e-30),1)
d=(eC-eA)*g.weights
for lo,hi in [(0,0.1),(0.1,0.3),(0.3,0.6),(0.6,1),(1,2),(2,4),(4,100)]:
    m=(r>=lo)&(r<hi)
    print(f"r {lo}-{hi}: dE={d[m].sum()*627.5:8.2f} kcal  <zeta>={np.average(zeta[m],weights=(n*g.weights)[m]):.3f} <zt>={np.average(zt[m],weights=(n*g.weights)[m]):.3f}")
