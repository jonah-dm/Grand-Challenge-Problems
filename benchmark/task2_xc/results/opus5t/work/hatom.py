import sys, importlib.util
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location('s',sys.argv[1]); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for bas in ('def2-svp','aug-cc-pvtz'):
    mol=gto.M(atom='H 0 0 0',basis=bas,spin=1,verbose=0)
    mf=dft.UKS(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB); mf.conv_tol=1e-9
    mf.run(); print(bas, mf.e_tot, mf.converged, mf.mo_energy[0][:2])
