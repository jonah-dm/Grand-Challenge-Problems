exec(open('hdiag.py').read().split('mf2=')[0])
S=mol.intor('int1e_ovlp')
c0=mf.mo_coeff[0][:,0]+0.05*mf.mo_coeff[0][:,1]; c0/=np.sqrt(c0@S@c0)
d=np.random.default_rng(0).normal(size=c0.size)*0.01
for nm,m in (('off',mo),('zeta2',mz)):
    m2=dft.UKS(mol); m2.grids.level=4; m2.xc='HF'; m2=m2.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB); m2.grids.build()
    def exc_of(c):
        c=c/np.sqrt(c@S@c); D=np.array([np.outer(c,c),0*S]); return m2._numint.nr_uks(mol,m2.grids,m2.xc,D)[1]
    h=1e-4; num=(exc_of(c0+h*d)-exc_of(c0-h*d))/(2*h)
    D=np.array([np.outer(c0,c0),0*S]); V=m2._numint.nr_uks(mol,m2.grids,m2.xc,D)[2][0]
    # derivative of normalized c: dc_n = d - c0 (c0 S d)
    dn=d-c0*(c0@S@d)
    print(nm,'num',num,'matrix',2*dn@V@c0)
