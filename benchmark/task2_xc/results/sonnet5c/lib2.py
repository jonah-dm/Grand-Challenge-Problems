"""Generate library of exact one-orbital densities (eigenstates of one-electron Hamiltonians
with bare nuclei) and tabulate data needed for exchange-enhancement design. No SCF involved."""
import numpy as np, pickle
from pyscf import gto, dft
from pyscf.dft import numint
def basis(l_max=2):
    ex_s=[0.04*2.5**k for k in range(12)]
    ex_p=[0.08*2.5**k for k in range(7)]
    ex_d=[0.15*3.0**k for k in range(4)]
    b=[[0,[e,1.0]] for e in ex_s]+[[1,[e,1.0]] for e in ex_p]+[[2,[e,1.0]] for e in ex_d]
    return b
sym={1:'H',2:'He',3:'Li'}
frames=[]
frames.append([(1,(0,0,0))])
frames.append([(3,(0,0,0))])
for R in (1.0,2.0,3.0,4.5,7.0): frames.append([(1,(0,0,0)),(1,(0,0,R))])
for R in (2.5,4.0): frames.append([(3,(0,0,0)),(1,(0,0,R))])
for R in (2.0,3.5): frames.append([(2,(0,0,0)),(1,(0,0,R))])
rows=[]
for fr in frames:
    atoms=[(sym[z],c) for z,c in fr]
    mol=gto.M(atom=atoms,basis={s:basis() for s in set(sym[z] for z,_ in fr)},unit='Bohr',spin=1 if sum(z for z,_ in fr)%2 else 0,verbose=0)
    mol.cart=False
    T=mol.intor('int1e_kin'); S=mol.intor('int1e_ovlp')
    V=0
    for z,c in fr:
        mol.set_rinv_origin(c)
        V=V-z*mol.intor('int1e_rinv')
    import scipy.linalg as sl
    e,C=sl.eigh(T+V,S)
    from pyscf.scf import hf as _hf
    g=dft.gen_grid.Grids(mol); g.level=6; g.build()
    ao=numint.eval_ao(mol,g.coords,deriv=2)
    nst=1
    for k in range(nst):
        c=C[:,k]
        D=np.outer(c,c); vj=_hf.get_jk(mol,D,hermi=1,with_k=False)[0] if False else mol.intor('int2e') if False else None; J=0.5*np.einsum('ij,ij',D,_hf.get_jk(mol,D,hermi=1)[0])
        phi=ao[0]@c; dphi=np.stack([ao[i]@c for i in (1,2,3)]); lphi=(ao[4]+ao[7]+ao[9])@c
        n=phi**2; gn=2*phi*dphi
        sig=np.sum(gn**2,0)
        rows.append(dict(frame=str(fr),k=k,e=e[k],J=J,w=g.weights,n=n,sig=sig,lap=2*np.sum(dphi**2,0)+2*phi*lphi))
        print(fr[0][0],len(fr),k,round(e[k],4),round(J,5),round(np.sum(g.weights*n),5),flush=True)
pickle.dump(rows,open('lib2.pkl','wb'))
