"""Score an XC functional on a reference set.

usage: python score.py <set_module> <reference.json> <functional> [--basis def2-qzvp] [--out results.json]

<functional> is either a PySCF/libxc functional string (e.g. 'pbe', 'b3lyp', 'r2scan',
'wb97m-v') or 'file:path/to/submission.py'. A submission module must define
    eval_xc(xc_code, rho, spin, relativity, deriv, omega, verbose) -> (exc, vxc, fxc, kxc)
following pyscf.dft.libxc.define_xc_ conventions, and optionally
    HYB  = global exact-exchange fraction (default 0.0)
    RSH  = (omega, alpha, beta) range separation: alpha*HF(full) + beta*HF(short-range, erfc(omega r)/r)
           (default (0,0,0)); eval_xc then receives omega and must return the matching semilocal part
    NLC  = (b, C) parameters of a VV10 nonlocal correlation term added by the harness (default None)
    META = True if eval_xc uses tau
Metrics: WTMAD (classes scaled by mean |reference|) and MCM (plain mean of class MAEs, kcal/mol).
MCM is the ranking metric on hidden set B, where small-magnitude classes dominate WTMAD.

Metric: MAE per class and weighted total MAD (WTMAD-2 style: each class scaled by its
mean |reference| so energy classes with large magnitudes do not dominate).
"""
import argparse, importlib, importlib.util, json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "kit"))
from pyscf import gto, dft, lib
from geom import to_str

HARTREE2KCAL = 627.5094740631


def load_submission(path):
    spec = importlib.util.spec_from_file_location("submission", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_species(atoms, charge, spin, basis, xc, sub):
    mol = gto.M(atom=to_str(atoms), basis=basis, charge=charge, spin=spin, verbose=0)
    mf = (dft.RKS if spin == 0 else dft.UKS)(mol)
    mf.grids.level = 4
    if sub is not None:
        hyb = float(getattr(sub, "HYB", 0.0))
        rsh = tuple(getattr(sub, "RSH", (0.0, 0.0, 0.0)))   # (omega, alpha, beta): alpha*HF + beta*SR-HF
        nlc = getattr(sub, "NLC", None)   # (b, C) VV10 parameters, or None
        # PySCF decides whether to build exact exchange / VV10 from mf.xc's *name*, not from the
        # coefficients passed to define_xc_; a dummy name of the right kind makes it honor them.
        if nlc:
            mf.xc = "wB97X-V"
        elif hyb > 0 or rsh[0] != 0:
            mf.xc = "HF" if rsh[0] == 0 else "wB97X"
        meta = bool(getattr(sub, "META", False))
        fn = sub.eval_xc
        if meta:
            # PySCF drops None entries of vxc and then truncates to (rho, sigma, tau): a zeros array in the
            # Laplacian slot would displace vtau. Normalize the Laplacian slot to None.
            def fn(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None, _f=sub.eval_xc):
                out = list(_f(xc_code, rho, spin, relativity, deriv, omega, verbose))
                if deriv > 0 and out[1] is not None and len(out[1]) == 4:
                    v = list(out[1]); v[2] = None; out[1] = tuple(v)
                return tuple(out)
        mf = mf.define_xc_(fn, "MGGA" if meta else "GGA", hyb=hyb, rsh=rsh)
        if nlc:
            mf.nlc = "VV10"; mf.nlcgrids.level = 1
            mf._numint.nlc_coeff = lambda *a, **k: (((nlc[0], nlc[1]), 1.0),)
    else:
        mf.xc = xc
        if "-v" in xc.lower() or "vv10" in xc.lower():
            mf.nlc = "vv10"
    mf.conv_tol = 1e-8
    mf.max_cycle = 100
    t = time.time()
    mf.run()
    if not mf.converged:  # derivative-free fallback first (custom functionals rarely provide fxc)
        mf.damp, mf.level_shift, mf.diis_space, mf.max_cycle = 0.5, 0.5, 12, 300
        mf.run(dm0=mf.make_rdm1())
    if not mf.converged:
        try:
            mf = mf.newton().run()
        except Exception:
            pass
    return mf.e_tot, mf.converged, time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set"); ap.add_argument("reference"); ap.add_argument("functional")
    ap.add_argument("--basis", default="def2-qzvp"); ap.add_argument("--out")
    a = ap.parse_args()
    lib.num_threads(int(os.environ.get("NTHREADS", "4")))
    S = importlib.import_module(a.set)
    ref = json.load(open(a.reference))["items"]
    sub = load_submission(a.functional[5:]) if a.functional.startswith("file:") else None
    E, meta = {}, {}
    for name, (atoms, charge, spin) in S.SPECIES.items():
        E[name], conv, dt = run_species(atoms, charge, spin, a.basis, a.functional, sub)
        meta[name] = {"converged": conv, "seconds": dt}
        print(f"  {name:14s} {E[name]:16.8f} {'' if conv else 'NOT CONVERGED'} {dt:5.1f}s", flush=True)
    errs, by_class = {}, {}
    for k, (cls, expr) in S.ITEMS.items():
        val = sum(c * E[s] for s, c in expr.items()) * HARTREE2KCAL
        errs[k] = {"class": cls, "value": val, "ref": ref[k]["kcal"], "error": val - ref[k]["kcal"]}
        by_class.setdefault(cls, []).append(errs[k])
    summary = {}
    for cls, rows in by_class.items():
        mae = sum(abs(r["error"]) for r in rows) / len(rows)
        mean_abs_ref = sum(abs(r["ref"]) for r in rows) / len(rows)
        summary[cls] = {"MAE": mae, "n": len(rows), "scaled": mae / mean_abs_ref * 100}
    wtmad = sum(v["scaled"] * v["n"] for v in summary.values()) / sum(v["n"] for v in summary.values())
    mean_class_mae = sum(v["MAE"] for v in summary.values()) / len(summary)
    out = {"functional": a.functional, "basis": a.basis, "energies": E, "meta": meta,
           "items": errs, "summary": summary, "WTMAD": wtmad, "MCM": mean_class_mae,
           "all_converged": all(m["converged"] for m in meta.values()),
           "total_seconds": sum(m["seconds"] for m in meta.values())}
    print(f"\n{a.functional}: " + "  ".join(f"{c} MAE={v['MAE']:.2f}" for c, v in summary.items())
          + f"  WTMAD={wtmad:.2f}  MCM={mean_class_mae:.2f}  time={out['total_seconds']:.0f}s")
    if a.out:
        json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
