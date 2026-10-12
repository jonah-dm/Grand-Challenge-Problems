"""Compute CCSD(T)/CBS reference energies for a species set.

E_ref = E_HF/cc-pVQZ + E_corr extrapolated from cc-pVTZ/cc-pVQZ (Helgaker X^-3),
frozen-core, nonrelativistic, at fixed geometry (no ZPE).

usage: python reference.py <set_module> <out.json>
"""
import importlib, json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "kit"))
from pyscf import gto, scf, cc, lib
from geom import to_str

HARTREE2KCAL = 627.5094740631


def ccsd_t(atoms, charge, spin, basis):
    mol = gto.M(atom=to_str(atoms), basis=basis, charge=charge, spin=spin, verbose=0)
    mf = (scf.RHF if spin == 0 else scf.UHF)(mol)
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.run()
    if not mf.converged:
        mf = scf.newton(mf).run()
    ncore = sum(1 for a, _ in atoms if a not in ("H", "He"))  # 1s core per heavy atom
    if ncore == 0 or mol.nelectron - 2 * ncore <= 0:
        frozen = 0 if ncore == 0 else None
    else:
        frozen = ncore
    if mol.nelectron - 2 * ncore <= 0 or mol.nelectron <= 1:
        return mf.e_tot, 0.0
    mycc = cc.CCSD(mf, frozen=frozen)
    mycc.conv_tol = 1e-8
    mycc.max_cycle = 200
    mycc.run()
    return mf.e_tot, mycc.e_corr + mycc.ccsd_t()


def main(setname, out):
    lib.num_threads(int(os.environ.get("NTHREADS", "4")))
    S = importlib.import_module(setname)
    res = json.load(open(out)) if os.path.exists(out) else {}
    for name, (atoms, charge, spin) in S.SPECIES.items():
        if name in res:
            continue
        t = time.time()
        hf3, c3 = ccsd_t(atoms, charge, spin, "cc-pvtz")
        hf4, c4 = ccsd_t(atoms, charge, spin, "cc-pvqz")
        ccbs = (4**3 * c4 - 3**3 * c3) / (4**3 - 3**3)
        res[name] = {"E_tot": hf4 + ccbs, "hf_qz": hf4, "corr_tz": c3, "corr_qz": c4}
        print(f"{name:14s} {hf4 + ccbs:16.8f}  {time.time() - t:6.1f}s", flush=True)
        json.dump(res, open(out, "w"), indent=1)
    items = {}
    for k, (cls, expr) in S.ITEMS.items():
        items[k] = {"class": cls, "kcal": sum(c * res[s]["E_tot"] for s, c in expr.items()) * HARTREE2KCAL}
    json.dump({"species": res, "items": items}, open(out.replace(".json", "_items.json"), "w"), indent=1)
    for k, v in items.items():
        print(f"{k:16s} {v['kcal']:9.2f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
