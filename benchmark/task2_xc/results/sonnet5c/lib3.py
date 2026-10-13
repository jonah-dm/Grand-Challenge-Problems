"""HF one-orbital densities: He 1s (spin orbital of closed-shell He), Ne 1s, Ar 1s (core orbitals as isolated
one-orbital spin densities).  J exact in the same basis.  Saves lib3.pkl with n, sig, lap, weights, J."""
import numpy as np, pickle
from pyscf import gto, scf, dft
from pyscf.dft import numint
from pyscf.scf import hf as _hf
out=[]
for el,basis in (('He','aug-cc-pv5z'),('Ne','def2-qzvp'),('Ar','def2-qzvp'),('Be','def2-qzvp'),('Li','def2-qzvp')):
    mol=gto.M(atom=el+' 0 0 0',basis=basis,spin=1 if el=='Li' else 0,verbose=0)
    mf=(scf.UHF if el=='Li' else scf.RHF)(mol); mf.conv_tol=1e-11; mf.kernel()
    C=mf.mo_coeff[0] if el=='Li' else mf.mo_coeff
    c=C[:,0]
    D=np.outer(c,c); J=0.5*np.einsum('ij,ij',D,_hf.get_jk(mol,D,hermi=1,with_k=False)[0] if False else mf.get_j(mol,D))
    g=dft.gen_grid.Grids(mol); g.level=7; g.build()
    ao=numint.eval_ao(mol,g.coords,deriv=2)
    phi=ao[0]@c; dphi=np.stack([ao[i]@c for i in (1,2,3)]); lphi=(ao[4]+ao[7]+ao[9])@c
    n=phi**2; sig=np.sum((2*phi*dphi)**2,0); lap=2*np.sum(dphi**2,0)+2*phi*lphi
    out.append(dict(frame=el+'_1s',J=J,w=g.weights,n=n,sig=sig,lap=lap,k=0,e=mf.mo_energy[0] if el!='Li' else mf.mo_energy[0][0]))
    print(el,J,np.sum(g.weights*n))
pickle.dump(out,open('lib3.pkl','wb'))
