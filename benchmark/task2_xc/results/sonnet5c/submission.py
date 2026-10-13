"""sonnet5c: one-orbital-branch hybrid meta-GGA (see NOTES.md for derivation and candidate list).

E_xc = 1/4 E_x^HF + 3/4 E_x^SL + E_c^G
 E_x^SL : spin-resolved (exact spin scaling), channel enhancement
          F = (1-w) F_gga(s) + w F_1(s),  w = z_s^2,  z_s = tau_W,s / tau_s
 F_gga  : PBE-form, kappa = 0.804 (Lieb-Oxford), mu = beta pi^2/3 (GE2-cancellation with PBE-type correlation)
 F_1    : one-orbital branch, constants solved at import from exact one-orbital data (see _one_orbital_constants)
 E_c^G  : PBE correlation with Stoll same-spin removal in one-orbital channels, pair-aware (prior art, NOTES)
All switches via env vars are only for the disclosed candidate comparison; defaults = submitted form.
"""
import os
import numpy as np

X1 = os.environ.get('S5_X1', 'Hs')            # one-orbital branch: 'none' | 'H' | 'HHe' | 'LIB'
CORR = os.environ.get('S5_CORR', 'pair')     # 'pbe' | 'stoll' | 'pair'
ZX = os.environ.get('S5_ZX', 'tot')         # 'chan': per-spin-channel z (exact spin scaling) | 'tot': total z
PW = float(os.environ.get('S5_WP', '2'))     # w = z^PW   (2: smallest integer keeping GE2 of the GGA branch)

HYB = 0.25
META = True
_nlc = os.environ.get('S5_NLC', '')
if _nlc:
    NLC = tuple(float(x) for x in _nlc.split(','))

KAPPA = 0.804
BETA = 0.06672455060314922
MU = BETA * np.pi**2 / 3
GAMMA = (1 - np.log(2)) / np.pi**2
EPS = 1e-14
_CX = -0.75 * (3 / np.pi) ** (1 / 3)


def _h1(s2, p):
    """one-orbital enhancement increment h(s) = a s^2 / (1 + b s^4)   (tail ~ s^-2: fractional-amplitude scaling)"""
    a, b = p
    if X1 == 'Hs':     # s^-1/2 tail, constant chosen so that max_s h = kappa (Lieb-Oxford saturated)
        return a * s2 / (1 + b * s2 * np.sqrt(s2))
    return a * s2 / (1 + b * s2 * s2)


def _H_exact_b():
    """b such that the exchange energy of the H atom 1s density equals -5/16 Ha, with a = mu."""
    from scipy.optimize import brentq
    r = np.linspace(1e-4, 40, 40001)
    w = 4 * np.pi * r**2 * (r[1] - r[0])
    n = np.exp(-2 * r) / np.pi              # spin density of H (fully polarised)
    gn = 2 * n
    n2 = 2 * n
    kf = (3 * np.pi**2 * n2) ** (1 / 3)
    s2 = (4 * gn**2) / (2 * kf * n2) ** 2
    exl = 0.5 * _CX * n2 ** (4 / 3)
    f = lambda b: np.sum(w * exl * (1 + _h1(s2, (MU, b)))) + 5.0 / 16
    return brentq(f, 1e-4, 50.0)


def _one_orbital_constants():
    if X1 == 'H':
        return (MU, _H_exact_b())
    if X1 == 'Hs':
        return (MU, 4.0 / (5.0 * KAPPA / MU) ** 1.25)
    if X1 == 'HHe':       # H, He, Ne1s, Ar1s, Be1s, Li1s HF orbitals all exact (fit7.py)
        return (1.9843, 16.6604)
    if X1 == 'LIB':       # + H2+ (R<=2), HeH2+ (fit7.py)
        return (2.59, 23.772)
    return (MU, 1.0)


_P1 = _one_orbital_constants() if X1 != 'none' else None


def _g(A, a1, b1, b2, b3, b4, rs):
    s = np.sqrt(rs)
    den = 2 * A * (b1 * s + b2 * rs + b3 * rs * s + b4 * rs * rs)
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / den)


def _pw92(rs, zeta):
    ec0 = _g(0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, rs)
    ec1 = _g(0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, rs)
    mac = -_g(0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, rs)
    fz = ((1 + zeta) ** (4 / 3) + (1 - zeta) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    z4 = zeta**4
    return ec0 + mac * fz * (1 - z4) / 1.709921 + (ec1 - ec0) * fz * z4


def _zsat(z, k=40.0):
    return 1 - np.log(1 + np.exp(k * (1 - z))) / k


def _ex_chan(r, sg, t, ztot=None):
    """exchange energy per volume of one spin channel (spin scaling)."""
    n2 = 2 * r
    kf = (3 * np.pi**2 * n2) ** (1 / 3)
    s2 = 4 * sg / (4 * kf**2 * n2**2)
    h_gga = KAPPA * (1 - 1 / (1 + MU * s2 / KAPPA))
    if _P1 is None:
        h = h_gga
    else:
        z = _zsat(sg / (8 * r * t)) if ztot is None else ztot
        w = z**PW
        h = (1 - w) * h_gga + w * _h1(s2, _P1)
    return 0.5 * _CX * n2 ** (4 / 3) * (1 + h)


def _H_pbe(n, sig, zeta, phi, ec):
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi**2 * ks2 * n**2)
    phi3 = phi**3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    return GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))


def _ec_pol(r, sg):
    rs = (3 / (4 * np.pi * r)) ** (1 / 3)
    ec = _pw92(rs, 1.0 + 0 * r)
    phi = 2 ** (-1 / 3)
    return r * (ec + _H_pbe(r, sg, 1.0, phi, ec))


def _E(ra, rb, saa, sab, sbb, ta, tb):
    n = ra + rb
    sig = saa + 2 * sab + sbb
    zt = _zsat(sig / (8 * n * (ta + tb))) if ZX == 'tot' else None
    ex = _ex_chan(ra, saa, ta, zt) + _ex_chan(rb, sbb, tb, zt)
    zeta = (ra - rb) / n
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _pw92(rs, zeta)
    phi = 0.5 * ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3))
    corr = n * (ec + _H_pbe(n, sig, zeta, phi, ec))
    if CORR != 'pbe':
        za = _zsat(saa / (8 * ra * ta)); zb = _zsat(sbb / (8 * rb * tb))
        for r, sg, zs, zo in ((ra, saa, za, zb), (rb, sbb, zb, za)):
            if CORR == 'stoll':
                corr = corr - zs**2 * _ec_pol(r, sg)
            else:
                corr = corr - zs**2 * (1 - zo**2 * (1 - zeta**2)) * _ec_pol(r, sg)
    return ex * (1 - HYB) + corr


def eval_xc(xc_code, rho, spin, relativity=0, deriv=1, omega=None, verbose=None):
    rho = np.asarray(rho)
    if spin == 0:
        n = rho[0]; sg = np.sum(rho[1:4] ** 2, axis=0); tau = rho[-1]
        ra = rb = n / 2; saa = sab = sbb = sg / 4; ta = tb = tau / 2
    else:
        ra, rb = rho[0][0], rho[1][0]
        ga, gb = rho[0][1:4], rho[1][1:4]
        saa, sbb, sab = np.sum(ga**2, 0), np.sum(gb**2, 0), np.sum(ga * gb, 0)
        ta, tb = rho[0][-1], rho[1][-1]
    N = len(ra)
    ra = np.maximum(ra, EPS); rb = np.maximum(rb, EPS)
    ta = np.maximum(ta, 1e-20); tb = np.maximum(tb, 1e-20)
    args = [np.asarray(a, dtype=float) for a in (ra, rb, saa, sab, sbb, ta, tb)]
    mask = (args[0] + args[1]) > 1e-11
    sel = [a[mask] for a in args]
    E0 = np.zeros(N); d = np.zeros((7, N))
    if mask.any():
        E0[mask] = np.real(_E(*sel))
        h = 1e-30
        for i in range(7):
            v = [s.astype(complex) for s in sel]
            v[i] = v[i] + 1j * h
            d[i, mask] = np.imag(_E(*v)) / h
    exc = E0 / (args[0] + args[1])
    if spin == 0:
        vxc = (0.5 * (d[0] + d[1]), 0.25 * (d[2] + d[3] + d[4]), np.zeros(N), 0.5 * (d[5] + d[6]))
    else:
        vxc = (np.stack([d[0], d[1]], 1), np.stack([d[2], d[3], d[4]], 1), np.zeros((N, 2)), np.stack([d[5], d[6]], 1))
    return exc, vxc, None, None
