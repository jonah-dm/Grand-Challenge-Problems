"""H2+ at R=6 bohr (aug-cc-pVTZ, exact = UHF same basis) and H+H2 barrier (def2-QZVP, harness settings).
usage: O6_K=K1 python work/barrier.py"""
import sys, os
here = os.path.dirname(os.path.abspath(__file__)); root = os.path.abspath(os.path.join(here, '..', '..', '..'))
sys.path.insert(0, root); sys.path.insert(0, os.path.join(root, 'kit'))
import score
from pyscf import gto, scf
K = 627.5094740631
sub = score.load_submission(os.path.join(here, '..', 'submission.py')); print('cand', os.environ.get('O7_K', 'S5'), sub.OMEGA, sub.F_LO)
B = 0.529177210903
out = {}
for R in [2.0, 6.0]:
    h2p = [('H', (0, 0, 0)), ('H', (0, 0, R * B))]
    e, c, _ = score.run_species(h2p, 1, 1, 'aug-cc-pvtz', None, sub)
    eh, c2, _ = score.run_species([('H', (0, 0, 0))], 0, 1, 'aug-cc-pvtz', None, sub)
    m = gto.M(atom='H 0 0 0; H 0 0 %f' % (R * B), basis='aug-cc-pvtz', charge=1, spin=1, verbose=0)
    ex = scf.UHF(m).run().e_tot
    mh = gto.M(atom='H 0 0 0', basis='aug-cc-pvtz', spin=1, verbose=0)
    exh = scf.UHF(mh).run().e_tot
    print('H2+ R=%.1f bohr: DFA %.2f  exact %.2f  error %.2f kcal/mol  conv %s %s' % (R, (e - eh) * K, (ex - exh) * K, (e - eh - ex + exh) * K, c, c2), flush=True)
r = 0.930
ts = [('H', (0, 0, -r)), ('H', (0, 0, 0)), ('H', (0, 0, r))]
from geom import diatomic
eTS, cT, _ = score.run_species(ts, 0, 1, 'def2-qzvp', None, sub)
eH2, c1, _ = score.run_species(diatomic('H', 'H', 0.741), 0, 0, 'def2-qzvp', None, sub)
eH, c2, _ = score.run_species([('H', (0, 0, 0))], 0, 1, 'def2-qzvp', None, sub)
print('H+H2 barrier (def2-QZVP): %.2f kcal/mol (ref 9.6)  conv %s %s %s  E(H)=%.7f' % ((eTS - eH2 - eH) * K, cT, c1, c2, eH), flush=True)
