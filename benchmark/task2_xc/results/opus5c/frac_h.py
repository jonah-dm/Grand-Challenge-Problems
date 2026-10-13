import os, importlib.util, numpy as np
from pyscf import gto, dft
spec = importlib.util.spec_from_file_location("sub", "submission.py"); sub = importlib.util.module_from_spec(spec); spec.loader.exec_module(sub)
mol = gto.M(atom="H 0 0 0", basis="aug-cc-pvtz", spin=1, verbose=0)
for xc in ['sub','pbe0']:
    out=[]
    for q in [0.25,0.5,0.75,1.0]:
        mf = dft.UKS(mol); mf.grids.level=5
        if xc=='sub': mf.xc='HF'; mf=mf.define_xc_(sub.eval_xc,'MGGA',hyb=sub.HYB)
        else: mf.xc=xc
        def get_occ(mo_energy=None, mo_coeff=None, q=q):
            o=np.zeros_like(mo_energy); o[0][np.argmin(mo_energy[0])]=q; return o
        mf.get_occ=get_occ; mf.conv_tol=1e-9; mf.kernel()
        out.append(mf.e_tot)
    print(xc, "E(q)+q/2:", [round(e+q/2,5) for e,q in zip(out,[0.25,0.5,0.75,1.0])], " E(1/2)-E(1)/2 =", round(out[1]-out[3]/2,5))
