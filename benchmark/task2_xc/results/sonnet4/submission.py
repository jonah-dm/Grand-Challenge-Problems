"""sonnet4: 25% hybrid meta-GGA with (i) two-norm iso-orbital exchange, (ii) Stoll-gated correlation.
See NOTES.md.  Env switches XS4_X (x variant) and XS4_G (gate power) exist only for the disclosed ablation."""
import os
import numpy as np
from pyscf.dft import libxc

HYB = 0.25
META = True
_XV = os.environ.get("XS4_X", "twonorm")     # twonorm | pbe | hnorm
_GP = float(os.environ.get("XS4_G", "2"))   # gate exponent on z_s (2 = Stoll z^2); 0 = no gate

KAPPA = 0.804
BETA = 0.06672
MU = BETA * np.pi ** 2 / 3
# iso-orbital branch F1 = (1 + A s^2)/(1 + B s^(5/2)): A,B solved from exact one-electron norms
# (H atom Ex=-5/16 Ha, He HF-density Ex=-1.02573 Ha), see norms.py
A1, B1 = 5.69000301, 5.02070045
if _XV == "hnorm":      # H-norm only: PBE form with mu raised to meet H (ablation)
    pass
CX = -0.75 * (3 / np.pi) ** (1 / 3)
KF3 = (3 * np.pi ** 2) ** (1 / 3)
EPS = 1e-12


Z0 = 0.8


def _cap(z):
    """C^1 cap of the iso-orbital indicator at 1 (cubic Hermite on [Z0,1], identity below Z0)."""
    L = 1 - Z0
    u = (z - Z0) / L
    mid = Z0 + L * (u + u ** 2 - u ** 3)      # value; slope 1 at u=0, slope 0 at u=1
    return np.where(z.real <= Z0, z, np.where(z.real >= 1, 1.0 + 0 * z, mid))


def _dcap(z):
    u = (z - Z0) / (1 - Z0)
    return np.where(z <= Z0, 1.0, np.where(z >= 1, 0.0, 1 + 2 * u - 3 * u ** 2))


def _zs(n, sg, t):
    return _cap(sg / (8 * n * t))


def _ex_spin(n, sg, t):
    """exchange energy density (per volume) of one spin channel (complex-step safe)."""
    n = n + 0j if not np.iscomplexobj(n) else n
    nn = 2 * n
    sg4 = 4 * np.where(sg.real < 1e-20, 1e-20, sg)
    s2 = sg4 / (2 * KF3 * nn ** (4 / 3)) ** 2
    Fg = 1 + KAPPA - KAPPA / (1 + MU * s2 / KAPPA)
    if _XV == "pbe":
        F = Fg
    else:
        z = _zs(n, sg, t)
        F1 = (1 + A1 * s2) / (1 + B1 * s2 ** 1.25)
        F = (1 - z ** 2) * Fg + z ** 2 * F1
    return 0.5 * CX * nn ** (4 / 3) * F


def _ex_all(n, sg, t):
    """returns e, de/dn, de/dsig, de/dtau for one spin channel, zero where n tiny"""
    ok = n > EPS
    n_ = np.where(ok, n, 1.0); s_ = np.where(ok, sg, 0.0); t_ = np.where(ok, np.maximum(t, 1e-20), 1.0)
    h = 1e-30
    e = _ex_spin(n_ + 0j, s_ + 0j, t_ + 0j).real
    dn = _ex_spin(n_ + 1j * h * n_, s_ + 0j, t_ + 0j).imag / (h * n_)
    ds = _ex_spin(n_ + 0j, s_ + 1j * h * np.maximum(s_, 1e-10), t_ + 0j).imag / (h * np.maximum(s_, 1e-10))
    dt = _ex_spin(n_ + 0j, s_ + 0j, t_ + 1j * h * t_).imag / (h * t_)
    z = np.zeros_like(n)
    return [np.where(ok, v, 0.0) for v in (e, dn, ds, dt)]


def _pbec(ra, rb, saa, sab, sbb):
    """PBE correlation (per volume) with derivs wrt na, nb, saa, sab, sbb via libxc."""
    na = np.maximum(ra, 1e-14); nb = np.maximum(rb, 1e-14)
    ga = np.sqrt(saa); gb = np.sqrt(sbb)
    # place b gradient so that sab is reproduced
    gbx = np.where(ga > 1e-14, sab / np.maximum(ga, 1e-14), 0.0)
    gby = np.sqrt(np.maximum(sbb - gbx ** 2, 0))
    z = np.zeros_like(na)
    r = (np.array([na, ga, z, z]), np.array([nb, gbx, gby, z]))
    exc, v = libxc.eval_xc("GGA_C_PBE", r, spin=1, deriv=1)[:2]
    ok = (ra > EPS) & (rb > EPS)
    return (na + nb) * exc, v[0][:, 0], v[0][:, 1], v[1][:, 0], v[1][:, 1], v[1][:, 2]


def _pbec_ss(n, sg):
    """E_c^PBE[n_s, 0] and derivs wrt n_s, sigma_ss"""
    ok = n > EPS
    n_ = np.where(ok, n, 1e-14); g = np.sqrt(np.where(ok, sg, 0.0)); z = np.zeros_like(n_)
    r = (np.array([n_, g, z, z]), np.array([z, z, z, z]))
    exc, v = libxc.eval_xc("GGA_C_PBE", r, spin=1, deriv=1)[:2]
    e = n_ * exc
    return [np.where(ok, q, 0.0) for q in (e, v[0][:, 0], v[1][:, 0])]


def _core(ra, rb, saa, sab, sbb, ta, tb):
    N = ra.shape[0]
    # exchange
    exa = _ex_all(ra, saa, ta); exb = _ex_all(rb, sbb, tb)
    sc = 1 - HYB
    E = sc * (exa[0] + exb[0])
    va = sc * exa[1]; vb = sc * exb[1]
    Vaa = sc * exa[2]; Vbb = sc * exb[2]; Vab = np.zeros(N)
    Vta = sc * exa[3]; Vtb = sc * exb[3]
    # correlation
    ec, c_na, c_nb, c_aa, c_ab, c_bb = _pbec(ra, rb, saa, sab, sbb)
    E = E + ec; va = va + c_na; vb = vb + c_nb; Vaa = Vaa + c_aa; Vab = Vab + c_ab; Vbb = Vbb + c_bb
    if _GP > 0:
        for (n, sg, t, which) in ((ra, saa, ta, 0), (rb, sbb, tb, 1)):
            e, dn, ds = _pbec_ss(n, sg)
            ok = (n > EPS) & (t > 1e-20)
            nn = np.where(ok, n, 1.0); tt = np.where(ok, t, 1.0); ss = np.where(ok, sg, 0.0)
            z0_ = ss / (8 * nn * tt)
            z = _cap(z0_); dc = _dcap(z0_)
            G = np.where(ok, z ** _GP, 0.0)
            dG = np.where(ok, _GP * z ** (_GP - 1), 0.0)   # dG/dz
            dzn = -z0_ * dc / nn; dzt = -z0_ * dc / tt
            dzs = np.where(ss > 0, z0_ * dc / np.maximum(ss, 1e-30), dc / (8 * nn * tt))
            E = E - G * e
            dn_t = G * dn + e * dG * dzn
            ds_t = G * ds + e * dG * dzs
            dt_t = e * dG * dzt
            if which == 0:
                va = va - dn_t; Vaa = Vaa - ds_t; Vta = Vta - dt_t
            else:
                vb = vb - dn_t; Vbb = Vbb - ds_t; Vtb = Vtb - dt_t
    return E, va, vb, Vaa, Vab, Vbb, Vta, Vtb


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    if spin == 0:
        r = np.asarray(rho)
        n = r[0]; g = r[1:4]; sg = np.einsum("ik,ik->k", g, g); t = r[-1]
        E, va, vb, Vaa, Vab, Vbb, Vta, Vtb = _core(n / 2, n / 2, sg / 4, sg / 4, sg / 4, t / 2, t / 2)
        nn = np.maximum(n, 1e-300)
        exc = E / nn
        vrho = va
        vsig = 0.25 * (Vaa + Vab + Vbb)
        vtau = Vta
        z = np.zeros_like(n)
        return exc, (vrho, vsig, z, vtau), None, None
    ra_, rb_ = np.asarray(rho[0]), np.asarray(rho[1])
    ra, rb = ra_[0], rb_[0]
    ga, gb = ra_[1:4], rb_[1:4]
    saa = np.einsum("ik,ik->k", ga, ga); sbb = np.einsum("ik,ik->k", gb, gb); sab = np.einsum("ik,ik->k", ga, gb)
    ta, tb = ra_[-1], rb_[-1]
    E, va, vb, Vaa, Vab, Vbb, Vta, Vtb = _core(ra, rb, saa, sab, sbb, ta, tb)
    exc = E / np.maximum(ra + rb, 1e-300)
    z = np.zeros_like(ra)
    return (exc, (np.stack([va, vb], 1), np.stack([Vaa, Vab, Vbb], 1), np.stack([z, z], 1), np.stack([Vta, Vtb], 1)), None, None)
