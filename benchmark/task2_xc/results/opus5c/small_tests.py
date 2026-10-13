"""Phase-2 checks on one/two-electron systems (not dev set): SCF energies and HOMO eigenvalues."""
import sys, os, importlib.util, numpy as np
from pyscf import gto, dft
spec = importlib.util.spec_from_file_location("sub", os.path.join(os.path.dirname(os.path.abspath(__file__)), "submission.py"))
sub = importlib.util.module_from_spec(spec); spec.loader.exec_module(sub)
def run(atom, spin, charge, basis, xc, unit='Angstrom'):
    mol = gto.M(atom=atom, basis=basis, spin=spin, charge=charge, unit=unit, verbose=0)
    mf = (dft.RKS if spin == 0 else dft.UKS)(mol); mf.grids.level = 5
    if xc == 'sub':
        mf.xc = 'HF'; mf = mf.define_xc_(sub.eval_xc, 'MGGA', hyb=sub.HYB)
    else:
        mf.xc = xc
    mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.kernel()
    eh = mf.mo_energy[0][mf.mo_occ[0] > 0].max() if spin else mf.mo_energy[mf.mo_occ > 0].max()
    return mf.e_tot, eh, mf.converged
basis = sys.argv[1] if len(sys.argv) > 1 else 'aug-cc-pvqz'
for label, atom, spin, charge, unit, exact in [("H atom", "H 0 0 0", 1, 0, 'Angstrom', -0.5),
                                ("H2+ R=2.0 bohr", "H 0 0 0; H 0 0 2.0", 1, 1, 'Bohr', -0.6026342),
                                ("He atom", "He 0 0 0", 0, 0, 'Angstrom', -2.903724),
                                ("H2 R=1.4 bohr", "H 0 0 0; H 0 0 1.4", 0, 0, 'Bohr', -1.174476)]:
    for xc in ['sub', 'pbe0']:
        e, eh, conv = run(atom, spin, charge, basis, xc, unit)
        print(f"{label:16s} {xc:5s} E={e:.6f} (exact {exact:.6f}, err {627.5*(e-exact):+.2f} kcal/mol)  eps_HOMO={eh:.4f}  conv={conv}", flush=True)
