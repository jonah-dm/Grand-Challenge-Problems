"""C1 violation table on the exact fractional H-atom density lambda*e^{-2r}/pi (fully polarized).
f(lam) = E_xc[lam n_H]/E_x[n_H]; exact f = lam^2 (E_xc[n_H] = -5/16 at lam=1)."""
import numpy as np, importlib.util, os
from pyscf.dft import libxc
r=np.linspace(1e-6,40,400001); dr=r[1]-r[0]; wts=4*np.pi*r*r*dr
def rho_pol(lam, meta):
    n=lam*np.exp(-2*r)/np.pi; g=-2*n; tau=g*g/(8*n)
    z=np.zeros_like(n); nb=np.full_like(n,1e-30)
    a=np.array([n,g,z,z,tau]) ; b=np.array([nb,z,z,z,z])
    return np.array([a,b])
def Exc_libxc(code, lam):
    hyb=libxc.hybrid_coeff(code); meta=libxc.xc_type(code)=='MGGA'
    rho=rho_pol(lam,meta)
    t=libxc.xc_type(code)
    if t=='LDA': rho=rho[:,:1]
    elif t!='MGGA': rho=rho[:,:4]
    exc=libxc.eval_xc(code,rho,spin=1,deriv=0)[0]
    return np.sum(exc*rho[0,0]*wts)+hyb*(-lam**2*5/16)
def load(p):
    s=importlib.util.spec_from_file_location('m',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def Exc_sub(m, lam):
    rho=rho_pol(lam,True); rho[1,0]=1e-30
    exc=m.eval_xc('',rho,1,0,1,None,0)[0]
    return np.sum(exc*rho[0,0]*wts)+m.HYB*(-lam**2*5/16)
rows=[]
for code in ['lda,','pbe,','b88,','tpss,','scan,','r2scan,','pbe0','b3lyp','pbe,pbe','r2scan']:
    f=[Exc_libxc(code,l) for l in (1,0.5,0.25)]
    rows.append((code,f))
B='/home/user/Grand-Challenge-Problems/benchmark/task2_xc/results/'
for nm,p,env in [('sonnet3',B+'sonnet3/submission.py',None),('opus4',B+'opus4/submission.py',None),('opus5t',B+'opus5t/submission.py',None)]:
    m=load(p); rows.append((nm,[Exc_sub(m,l) for l in (1,0.5,0.25)]))
print('%-10s %10s %10s %10s %10s %10s'%('functional','Exc[H]','f(1/2)','viol(1/2)','f(1/4)','dE H2+ inf kcal'))
for nm,f in rows:
    f1,fh,fq=f; F1=fh/(-0.3125); F2=fq/(-0.3125)
    # stretched H2+ limit: 2*E_xc[n_H/2] vs exact -5/32 (E_xc only; T,V,J exact for exact density)
    print('%-10s %10.4f %10.4f %10.4f %10.4f %10.1f'%(nm,f1,F1,F1-0.25,F2,(2*fh+5/32)*627.5095))
print('exact      %10.4f %10.4f %10.4f %10.4f %10.1f'%(-0.3125,0.25,0,0.0625,0))
