exec(open('hdiag.py').read().split('mf2=')[0])
S=mol.intor("int1e_ovlp"); c0=mf.mo_coeff[0][:,0]+0.05*mf.mo_coeff[0][:,1]; c0=c0/np.sqrt(c0@S@c0)
rng=np.random.default_rng(0); d=rng.normal(size=c0.size)*0.01
for nm,m in (('off',mo),('zeta2',mz)):
    m2=dft.UKS(mol); m2.grids.level=4; m2.xc='HF'; m2=m2.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB)
    def en(c):
        c=c/np.sqrt(c@S@c); return m2.energy_tot(np.array([np.outer(c,c),np.zeros_like(S)]))
    h=1e-4; num=(en(c0+h*d)-en(c0-h*d))/(2*h)
    f=m2.get_fock(dm=np.array([np.outer(c0,c0),np.zeros_like(S)]))[0]
    # dE/dc = 2 (F c - eps S c) projected, for normalized c0
    g=2*(f@c0-(c0@f@c0)*(S@c0)); print(nm,'numeric',num,'fock',g@d)
