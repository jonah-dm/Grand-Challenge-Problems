import numpy as np
from pyscf import gto, scf, dft
from scipy.optimize import brentq, fsolve
# He HF density, exact exchange
m=gto.M(atom='He 0 0 0',basis='aug-cc-pv5z',verbose=0); mf=scf.RHF(m).run()
g=dft.gen_grid.Grids(m); g.level=9; g.build()
dm=mf.make_rdm1()
ao=dft.numint.eval_ao(m,g.coords,deriv=1)
rho=dft.numint.eval_rho(m,ao,dm,xctype='GGA')
n=rho[0]; sig=np.einsum('ik,ik->k',rho[1:4],rho[1:4]); w=g.weights
Ex_He=-0.25*2*np.einsum('ij,ji',dm,np.einsum('ijkl,kl->ij',m.intor('int2e'),dm)*0+0) if False else None
vj,vk=scf.hf.get_jk(m,dm)
Ex_He=-0.25*np.einsum('ij,ji',dm,vk); print('Ex He exact',Ex_He, 'N',np.sum(n*w))
# H atom exact density
r=np.linspace(1e-6,40,400001); nH=np.exp(-2*r)/np.pi; sH=(2*nH)**2*1.0
wH=4*np.pi*r**2*(r[1]-r[0])
def exu(n): return -0.75*(3/np.pi)**(1/3)*n**(4/3)
def s_of(n,sg): return np.sqrt(sg)/(2*(3*np.pi**2)**(1/3)*n**(4/3))
def F1(s,A,B): return (1+A*s**2)/(1+B*s**2.5)
def hnorm(A,B):  # one spin channel, n_s=nH: Ex=1/2 Ex_unpol[2 n_s]
    nn=2*nH; ss=4*sH
    return np.sum(0.5*exu(nn)*F1(s_of(nn,ss),A,B)*wH)
def henorm(A,B):
    msk=n>1e-12
    return np.sum(exu(n[msk])*F1(s_of(n[msk],sig[msk]),A,B)*w[msk])
print('H unif, GGA0', hnorm(0,0)); print('He',henorm(0,0))
f=lambda p:[hnorm(*p)+0.3125, henorm(*p)-Ex_He]
for g0 in [(0.2,0.1),(0.5,0.5),(1,1),(0.3,0.3)]:
    sol,info,ier,msg=fsolve(f,g0,full_output=True); print(g0,sol,ier,f(sol))
