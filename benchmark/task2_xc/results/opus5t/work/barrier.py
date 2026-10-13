"""Self-generated validation (not dev set): H+H2 barrier (exact 9.6), H2+ at R=4 bohr fractional test."""
import os, sys, importlib.util
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location('s',sys.argv[1]); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def E(atom, spin, charge=0, unit='A'):
    mol=gto.M(atom=atom,basis='aug-cc-pvtz',spin=spin,charge=charge,verbose=0,unit=unit)
    mf=(dft.RKS if spin==0 else dft.UKS)(mol); mf.grids.level=4
    mf.xc='HF'
    mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB); mf.conv_tol=1e-9; mf.max_cycle=100; mf.run()
    if not mf.converged:
        mf.damp, mf.level_shift, mf.diis_space, mf.max_cycle = 0.5, 0.5, 12, 300; mf.run(dm0=mf.make_rdm1())
    if not mf.converged:
        mf = mf.newton().run()
    if not mf.converged: print('NOT CONVERGED', atom)
    return mf.e_tot
eh=E('H 0 0 0',1); eh2=E('H 0 0 0; H 0 0 0.741',0); ets=E('H 0 0 -0.93; H 0 0 0; H 0 0 0.93',1)
print(os.environ.get('O5_W'),'E(H) %.5f'%eh,'barrier %.2f kcal/mol (exact 9.6)'%((ets-eh-eh2)*627.5095), 'H2 De %.2f (exact 109.5)'%((2*eh-eh2)*627.5095),flush=True)
for R in (2.0,4.0,6.0):
    e=E('H 0 0 0; H 0 0 %f'%R,1,charge=1,unit='B')
    print('H2+ R=%.1f bohr: E-E(H) = %.2f kcal/mol'%(R,(e-eh)*627.5095),flush=True)
