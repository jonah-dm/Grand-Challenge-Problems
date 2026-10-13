exec(open('hdiag.py').read().split('mf2=')[0])
import scipy.optimize as so
# direct minimization over the occupied alpha orbital coefficients (normalized), to find the true minimum
mf2=dft.UKS(mol); mf2.grids.level=4; mf2.xc='HF'; mf2=mf2.define_xc_(mz.eval_xc,'MGGA',hyb=mz.HYB)
S=mol.intor('int1e_ovlp'); c0=mf.mo_coeff[0][:,0]
def en(c):
    c=c/np.sqrt(c@S@c); D=np.array([np.outer(c,c),np.zeros_like(S)]); return mf2.energy_tot(D)
print('start',en(c0))
res=so.minimize(en,c0,method='BFGS',options={'maxiter':60,'gtol':1e-6}); print('min',res.fun,res.nit)
c=res.x/np.sqrt(res.x@S@res.x)
mf3=dft.UKS(mol); mf3.grids.level=4; mf3.xc='HF'; mf3=mf3.define_xc_(mz.eval_xc,'MGGA',hyb=mz.HYB)
f=mf3.get_fock(dm=np.array([np.outer(c,c),np.zeros_like(S)]))
print('orbital gradient norm', np.linalg.norm((f[0]@c) - (c@f[0]@c)*(S@c)))
