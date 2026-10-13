"""H+H2 barrier (linear TS r=0.930 A, ref 9.6) and H2+ at R=6 bohr (exact = UHF same basis), harness settings."""
import sys, os, importlib.util
base = '/home/user/Grand-Challenge-Problems/benchmark/task2_xc'
sys.path.insert(0, base); sys.path.insert(0, base + '/kit')
import score
from pyscf import gto, scf
K = 627.5094740631
path = sys.argv[1]; basis = sys.argv[2] if len(sys.argv) > 2 else 'def2-qzvp'
sub = score.load_submission(path) if path != 'ref' else None
xc = None if sub else 'b3lyp'
R = 6 * 0.529177210903
geo = {'H': ([('H', (0, 0, 0))], 0, 1), 'H2': ([('H', (0, 0, 0)), ('H', (0, 0, 0.741))], 0, 0),
       'TS': ([('H', (0, 0, -0.930)), ('H', (0, 0, 0)), ('H', (0, 0, 0.930))], 0, 1),
       'H2p': ([('H', (0, 0, 0)), ('H', (0, 0, R))], 1, 1)}
E = {}
for k, (at, c, s) in geo.items():
    E[k], conv, dt = score.run_species(at, c, s, basis, xc, sub)
    print(k, E[k], conv, f"{dt:.0f}s", flush=True)
ex = {}
for k in ['H', 'H2p']:
    at, c, s = geo[k]
    mol = gto.M(atom=score.to_str(at), basis=basis, charge=c, spin=s, verbose=0)
    ex[k] = scf.UHF(mol).run().e_tot
print(f"{path}: barrier = {(E['TS']-E['H']-E['H2'])*K:.2f} (ref 9.6); H2+(6 bohr) - H = {(E['H2p']-E['H'])*K:.2f} "
      f"(exact {(ex['H2p']-ex['H'])*K:.2f}); E(H) = {E['H']:.6f}")
