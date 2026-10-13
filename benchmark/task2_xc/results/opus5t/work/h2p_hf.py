from pyscf import gto, scf
def E(atom,charge,spin):
    return scf.UHF(gto.M(atom=atom,basis='aug-cc-pvtz',charge=charge,spin=spin,unit='B',verbose=0)).run().e_tot
eh=E('H 0 0 0',0,1)
for R in (2.,4.,6.): print('exact(HF) H2+ R=%.1f: %.2f'%(R,(E('H 0 0 0; H 0 0 %f'%R,1,1)-eh)*627.5095))
