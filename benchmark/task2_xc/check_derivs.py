"""Check a submission's vxc against finite differences of its exc (spin-unpolarized, random points).

usage: python check_derivs.py file:results/<model>/submission.py

Wrong derivatives are the most common reason a new functional fails to converge or gives
nonsense energies. vxc must be d(rho*exc)/d(variable) in PySCF's convention:
  GGA : vxc = (vrho, vsigma)                 variables rho, sigma=|grad rho|^2
  MGGA: vxc = (vrho, vsigma, vlapl, vtau)    plus lapl, tau
"""
import sys, importlib.util
import numpy as np


def load(path):
    spec = importlib.util.spec_from_file_location("submission", path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def pack(rho, sigma, tau, meta):
    z = np.zeros_like(rho)
    if meta:
        return np.array([rho, np.sqrt(sigma), z, z, z, tau])
    return np.array([rho, np.sqrt(sigma), z, z])


def main():
    sub = load(sys.argv[1][5:] if sys.argv[1].startswith("file:") else sys.argv[1])
    meta = bool(getattr(sub, "META", False))
    rng = np.random.default_rng(1)
    n = 200
    rho = 10 ** rng.uniform(-3, 0.5, n)
    s = 10 ** rng.uniform(-1.5, 0.7, n)
    kf = (3 * np.pi**2 * rho) ** (1 / 3)
    sigma = (2 * kf * rho * s) ** 2
    tau_w = sigma / (8 * rho); tau_ueg = 0.3 * (3 * np.pi**2) ** (2 / 3) * rho ** (5 / 3)
    tau = tau_w + 10 ** rng.uniform(-1.5, 0.7, n) * tau_ueg

    def E(r, sg, t):
        return r * np.asarray(sub.eval_xc("", pack(r, sg, t, meta), 0, 0, 1, None, 0)[0])

    exc, vxc = sub.eval_xc("", pack(rho, sigma, tau, meta), 0, 0, 1, None, 0)[:2]
    vrho, vsigma = np.asarray(vxc[0]), np.asarray(vxc[1])
    ok = True
    for name, v_an, var in [("vrho", vrho, 0), ("vsigma", vsigma, 1)] + ([("vtau", np.asarray(vxc[3]), 2)] if meta else []):
        h = 1e-5 * np.maximum([rho, sigma, tau][var], 1e-8)
        args = [rho.copy(), sigma.copy(), tau.copy()]
        ap, am = [a.copy() for a in args], [a.copy() for a in args]
        ap[var] += h; am[var] -= h
        v_fd = (E(*ap) - E(*am)) / (2 * h)
        rel = np.abs(v_an - v_fd) / (np.abs(v_fd) + 1e-6)
        bad = np.mean(rel > 1e-3)
        print(f"{name:7s} max rel err {rel.max():.2e}   fraction >1e-3: {bad:.3f}")
        ok &= bad < 0.02
    # spin-polarized: rho has shape (2, nvar, n); vrho (n,2), vsigma (n,3) [aa, ab, bb], vtau (n,2)
    fa = rng.uniform(0.2, 0.8, n)
    ra, rb = rho * fa, rho * (1 - fa)
    sa, sb = sigma * fa**2, sigma * (1 - fa) ** 2
    sab = np.sqrt(sa * sb) * rng.uniform(-0.9, 0.9, n)   # |grad a||grad b| cos
    ta, tb = tau * fa, tau * (1 - fa)

    def packs(ra, rb, sa, sb, sab, ta, tb):
        z = np.zeros_like(ra)
        ga = np.sqrt(sa); gb_x = sab / np.maximum(ga, 1e-30); gb_y = np.sqrt(np.maximum(sb - gb_x**2, 0))
        a = np.array([ra, ga, z, z, z, ta]) if meta else np.array([ra, ga, z, z])
        b = np.array([rb, gb_x, gb_y, z, z, tb]) if meta else np.array([rb, gb_x, gb_y, z])
        return np.array([a, b])

    def Es(*v):
        r = v[0] + v[1]
        return r * np.asarray(sub.eval_xc("", packs(*v), 1, 0, 1, None, 0)[0])

    try:
        exc, vxc = sub.eval_xc("", packs(ra, rb, sa, sb, sab, ta, tb), 1, 0, 1, None, 0)[:2]
    except Exception as e:
        print("spin-polarized call FAILED:", e); print("DERIVATIVES WRONG (spin-polarized)"); return
    vrho = np.asarray(vxc[0]); vsig = np.asarray(vxc[1])
    checks = [("vrho_a", vrho[:, 0], 0), ("vrho_b", vrho[:, 1], 1), ("vsigma_aa", vsig[:, 0], 2),
              ("vsigma_ab", vsig[:, 1], 4), ("vsigma_bb", vsig[:, 2], 3)]
    if meta:
        vt = np.asarray(vxc[3]); checks += [("vtau_a", vt[:, 0], 5), ("vtau_b", vt[:, 1], 6)]
    base = [ra, rb, sa, sb, sab, ta, tb]
    for name, v_an, var in checks:
        h = 1e-5 * np.maximum(np.abs(base[var]), 1e-8)
        ap, am = [a.copy() for a in base], [a.copy() for a in base]
        ap[var] += h; am[var] -= h
        v_fd = (Es(*ap) - Es(*am)) / (2 * h)
        rel = np.abs(v_an - v_fd) / (np.abs(v_fd) + 1e-6)
        bad = np.mean(rel > 1e-3)
        print(f"{name:10s} max rel err {rel.max():.2e}   fraction >1e-3: {bad:.3f}")
        ok &= bad < 0.02
    print("DERIVATIVES OK" if ok else "DERIVATIVES WRONG: fix vxc before running the SCF")


if __name__ == "__main__":
    main()
