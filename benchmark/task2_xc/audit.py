"""Automated audit of a submission: wrapper/mix detection and exact-constraint checks.

usage: python audit.py file:results/<model>/submission.py [--out audit.json]

Wrapper check: evaluate the submission's exc on random (rho, sigma, tau) samples and fit it as a
linear combination of common libxc functionals (non-negative least squares). R^2 ~ 1 with a
sparse fit means the submission is (close to) an existing functional or a mix of them.

Constraint checks (heuristic, pointwise, spin-unpolarized):
  UEG   : at sigma=0, tau=tau_UEG the enhancement factor over LDA xc should be 1.
  LO    : Lieb-Oxford: e_xc >= -2.27 * e_x^LDA (pointwise form, lambda=2.27).
  NEG   : e_xc <= 0 everywhere.
  SMOOTH: finite, no NaN, over the sampled range.
"""
import argparse, importlib.util, json, sys, os
import numpy as np
from pyscf.dft import libxc
from scipy.optimize import nnls

LIB = ["LDA,PW", "PBE", "PBE0", "BLYP", "B3LYP", "TPSS", "SCAN", "R2SCAN", "M06-L", "M06", "M06-2X",
       "M05", "M05-2X", "M08-HX", "M11", "MN12-L", "MN15-L", "MN15", "PW6B95", "B97-1", "B97-2", "B97-D",
       "revPBE", "RPBE", "PW91", "BP86", "B88,P86", "PBEsol", "wB97X-V", "wB97M-V", "SCAN0", "HSE06",
       "B97M-V", "TPSSh", "revTPSS", "HCTH", "OLYP", "X3LYP", "BHANDHLYP", "CAM-B3LYP", "MPW1K", "BMK"]


def load(path):
    spec = importlib.util.spec_from_file_location("submission", path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def sample(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    rho = 10 ** rng.uniform(-4, 1.0, n)
    s = 10 ** rng.uniform(-2, 1.0, n)                   # reduced gradient
    kf = (3 * np.pi**2 * rho) ** (1 / 3)
    sigma = (2 * kf * rho * s) ** 2
    tau_w = sigma / (8 * rho)
    tau_ueg = 0.3 * (3 * np.pi**2) ** (2 / 3) * rho ** (5 / 3)
    alpha = 10 ** rng.uniform(-2, 1.0, n)
    tau = tau_w + alpha * tau_ueg
    lapl = np.zeros(n)
    return rho, sigma, tau, lapl, tau_ueg


def exc_of(fn_or_code, rho, sigma, tau, lapl, kind):
    """kind: 'LDA' | 'GGA' | 'MGGA' — the rho format the callee expects."""
    z = np.zeros_like(rho)
    if kind == "LDA":
        r = rho
    elif kind == "GGA":
        r = np.array([rho, np.sqrt(sigma), z, z])
    else:
        r = np.array([rho, np.sqrt(sigma), z, z, lapl, tau])
    if callable(fn_or_code):
        return fn_or_code("", r, 0, 0, 1, None, 0)[0]
    return libxc.eval_xc(fn_or_code, r, 0, 0, 1)[0]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("sub"); ap.add_argument("--out"); a = ap.parse_args()
    sub = load(a.sub[5:] if a.sub.startswith("file:") else a.sub)
    meta = bool(getattr(sub, "META", False))
    rho, sigma, tau, lapl, tau_ueg = sample()
    kind = "MGGA" if meta else "GGA"
    e_sub = np.asarray(exc_of(sub.eval_xc, rho, sigma, tau, lapl, kind))
    res = {"META": meta, "HYB": float(getattr(sub, "HYB", 0.0))}
    res["SMOOTH"] = bool(np.all(np.isfinite(e_sub)))
    res["NEG_frac_positive"] = float(np.mean(e_sub > 1e-10))
    ex_lda = -(3 / 4) * (3 / np.pi) ** (1 / 3) * rho ** (1 / 3)
    res["LO_frac_violating"] = float(np.mean(e_sub < 2.27 * ex_lda - 1e-10))
    # UEG limit: sigma -> 0, tau = tau_ueg
    n = 200
    r0 = 10 ** np.linspace(-4, 1, n)
    t0 = 0.3 * (3 * np.pi**2) ** (2 / 3) * r0 ** (5 / 3)
    e0 = np.asarray(exc_of(sub.eval_xc, r0, np.zeros(n) + 1e-30, t0, np.zeros(n), kind))
    e_lda = np.asarray(exc_of("LDA,PW", r0, np.zeros(n) + 1e-30, t0, np.zeros(n), "LDA"))
    # add exact exchange fraction back for hybrids: libxc energy excludes HF part; compare xc minus hyb*ex
    ex0 = np.asarray(exc_of("LDA_X", r0, np.zeros(n) + 1e-30, t0, np.zeros(n), "LDA"))
    e_lda_eff = e_lda - res["HYB"] * ex0
    res["UEG_max_rel_dev"] = float(np.max(np.abs(e0 - e_lda_eff) / np.abs(e_lda_eff)))
    # wrapper / mix detection
    cols, names = [], []
    for f in LIB:
        try:
            e = np.asarray(exc_of(f, rho, sigma, tau, lapl, libxc.xc_type(f)))
            if np.all(np.isfinite(e)):
                cols.append(e); names.append(f)
        except Exception:
            pass
    A = np.array(cols).T
    w, rnorm = nnls(A, e_sub)
    ss = np.sum((e_sub - e_sub.mean()) ** 2)
    res["mix_R2"] = float(1 - rnorm**2 / ss)
    res["mix_components"] = {n_: round(float(c), 4) for n_, c in zip(names, w) if c > 1e-3}
    best = None
    for n_, c in zip(names, A.T):
        r2 = 1 - np.sum((e_sub - c) ** 2) / ss
        if best is None or r2 > best[1]:
            best = (n_, float(r2))
    res["closest_single"] = {"name": best[0], "R2": best[1]}
    verdict = ("WRAPPER" if best[1] > 0.99999 else
               "MIX" if res["mix_R2"] > 0.9999 and len(res["mix_components"]) <= 3 else
               "NEAR-EXISTING" if res["mix_R2"] > 0.999 else "DISTINCT")
    res["verdict"] = verdict
    print(json.dumps(res, indent=1))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
