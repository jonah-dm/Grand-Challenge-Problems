import sys, json
sys.path.insert(0, "/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit")
sys.path.insert(0, "/home/user/Grand-Challenge-Problems/benchmark/task2_xc")
import importlib.util
from pyscf import gto, dft
from geom import to_str
spec = importlib.util.spec_from_file_location("s", "scan/a20.py"); s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
mol = gto.M(atom="H 0 0 0; H 0 0 0.741", basis="def2-svp", verbose=0)
for name, f in [("pbe", None), ("sub", s)]:
    mf = dft.RKS(mol); mf.grids.level = 3
    if f is None: mf.xc = "pbe"
    else: mf = mf.define_xc_(f.eval_xc, "GGA", hyb=f.HYB)
    print(name, mf.kernel())
