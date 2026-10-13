"""opus6c: LC hybrid whose short-range semilocal exchange carries a bounded, Laplacian-free branch that
is exactly amplitude-quadratic in one-orbital spin channels (see NOTES.md).

E_xc = E_x^{LR-HF}(w*) + a E_x^{SR-HF}(w*) + (1-a) E_x^{SR,sl}(w*) + E_c(opus4, prior art)
SR semilocal exchange per spin channel: (1-w_s) e^{SR,base} + w_s e^{SR,1}
  e^{SR,base}: sonnet3 z-blended enhancement (prior art), ITYH-attenuated (known)
  e^{SR,1}  = -(pi/(2 w*^2)) n_s^2 * 5 Q(b),  b = |grad n_s|^2/(4 w*^2 n_s^2),  Q = R(b)/b,
              R(b) = 1 - (32/(5 pi)) int_0^inf e^{-b u^2} (1+u^2)^-4 du   (exact SR self-exchange of every
              hydrogenic ion at every fractional occupation; bounded: sup = on-top value), capped smoothly at F<=1.804
  w_s = zeta^2 z_s^2 (C1 premise: single-orbital channel, other spin absent)
w* = 0.368708: smallest omega for which the exact-C1 SR branch respects F <= 1+kappa on the H atom (derived, NOTES L4).
Env O6_K selects the disclosed candidates (default K1 = submitted).
"""
import os
import numpy as np
from scipy.special import erf as _erf, erfcx as _erfcx, gamma as _gamma, comb as _comb

_K = os.environ.get('O6_K', 'K1')
OMEGA = float(os.environ.get('O6_OMEGA', '0.3687078276482022'))
A_SR = 0.0 if _K == 'K2' else 0.25          # SR exact-exchange fraction (PBE0 1/4; K2: pure LC)
HYB = A_SR
RSH = (OMEGA, 1.0, A_SR - 1.0)               # 100% LR-HF, A_SR SR-HF
META = True
F_LO = 1.804
ZREG = 1e-2
PW = 2.0
KAPPA, MU, BETA = 0.804, 0.2195149727645171, 0.06672455060314922
CTAIL = 4.0 / (5.0 * KAPPA / MU) ** 1.25
GAMMA = (1 - np.log(2)) / np.pi**2
EPS = 1e-14

# asymptotic coefficients of J4(b) for large b: 1/2 sum_k (-1)^k C(k+3,3) Gamma(k+1/2) b^(-k-1/2)
_ASY = [0.5 * (-1) ** k * _comb(k + 3, 3) * _gamma(k + 0.5) for k in range(14)]


def _J4(b):
    sb = np.sqrt(b)
    E = _erfcx(sb)
    J1 = np.pi / 2 * E
    J2 = np.pi / 4 * (1 - 2 * b) * E + 0.5 * np.sqrt(np.pi) * sb
    J3 = ((3 - 2 * b) * J2 + 2 * b * J1) / 4
    return ((5 - 2 * b) * J3 + 2 * b * J2) / 6


def _J4asy(b):
    out = 0 * b
    for k, c in enumerate(_ASY):
        out = out + c * b ** (-k - 0.5)
    return out


def Qfun(b):
    """Q(b) = R(b)/b, complex-step safe"""
    br = np.real(b)
    small = br < 1e-5
    large = br > 400.0
    mid = ~(small | large)
    bm = np.where(mid, b, 1.0 + 0 * b)
    bl = np.where(large, b, 1000.0 + 0 * b)
    Qm = (1 - 32 / (5 * np.pi) * _J4(bm)) / bm
    Ql = (1 - 32 / (5 * np.pi) * _J4asy(bl)) / bl
    Qs = 0.2 - b / 10 + b * b / 6
    return np.where(small, Qs, np.where(large, Ql, Qm))


def _att(a):
    """LDA short-range attenuation factor for erfc(omega r)/r (complex-step safe)"""
    big = np.real(a) > 3.0
    ab = np.where(big, 1.0 + 0 * a, a)
    asm = np.where(big, a, 1.0 + 0 * a)
    Fd = 1 - 8 / 3 * ab * (np.sqrt(np.pi) * _erf(1 / (2 * ab)) + (2 * ab - 4 * ab**3) * np.exp(-1 / (4 * ab**2)) - 3 * ab + 4 * ab**3)
    x = 1 / asm**2
    Fs = x / 36 - x**2 / 960 + x**3 / 26880 - x**4 / 829440
    return np.where(big, Fs, Fd)


def _g(A, a1, b1, b2, b3, b4, rs):
    s = np.sqrt(rs)
    den = 2 * A * (b1 * s + b2 * rs + b3 * rs * s + b4 * rs * rs)
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / den)


def _pw92(rs, zeta):
    ec0 = _g(0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, rs)
    ec1 = _g(0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, rs)
    mac = -_g(0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, rs)
    fz = ((1 + zeta) ** (4 / 3) + (1 - zeta) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    fpp0 = 1.709921
    z4 = zeta**4
    return ec0 + mac * fz * (1 - z4) / fpp0 + (ec1 - ec0) * fz * z4


def _zsat(z, k=40.0):
    return 1 - np.log(1 + np.exp(k * (1 - z))) / k


def _ex_base_sr(n2, s2, z):
    """SR exchange energy density of an unpolarized gas of density n2 (sonnet3 enhancement, ITYH)"""
    ex_lda = -0.75 * (3 / np.pi) ** (1 / 3) * n2 ** (4 / 3)
    h_gga = KAPPA * (1 - 1 / (1 + MU * s2 / KAPPA))
    s5 = s2 * np.sqrt(s2)
    h_orb = MU * s2 / (1 + CTAIL * s5)
    w = z**2
    h = (1 - w) * h_gga + w * h_orb
    kf = (3 * np.pi**2 * n2) ** (1 / 3)
    a = OMEGA * np.sqrt(1 + h) / (2 * kf)
    return ex_lda * (1 + h) * _att(a)


def _ex_branch(r, sig_s):
    """one-orbital SR branch for spin density r, in the unpolarized-density convention (n2 = 2r):
    returns e such that the channel energy density is 0.5*e"""
    n2 = 2 * r
    ex_lda = -0.75 * (3 / np.pi) ** (1 / 3) * n2 ** (4 / 3)
    k = (6 * np.pi**2 * r) ** (1 / 3)
    b = sig_s / (4 * OMEGA**2 * r * r)
    F = 5.0 / 9.0 * (k / OMEGA) ** 2 * Qfun(b)
    Fc = F * (1 + (F / F_LO) ** 8) ** (-1 / 8)
    return ex_lda * Fc


def _ec_pol(r, sg):
    rs = (3 / (4 * np.pi * r)) ** (1 / 3)
    ec = _pw92(rs, 1.0 + 0 * r)
    phi = 2 ** (-1 / 3)
    kf = (3 * np.pi**2 * r) ** (1 / 3)
    t2 = sg / (4 * phi**2 * (4 * kf / np.pi) * r**2)
    phi3 = phi**3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    H = GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))
    return r * (ec + H)


def _E(ra, rb, saa, sab, sbb, ta, tb):
    n = ra + rb
    sig = saa + 2 * sab + sbb
    z = _zsat(sig / (8 * n * (ta + tb)))
    zeta = (ra - rb) / n
    ex = 0.0
    for r, s, t in ((ra, saa, ta), (rb, sbb, tb)):
        n2 = 2 * r
        kf = (3 * np.pi**2 * n2) ** (1 / 3)
        s2 = 4 * s / (4 * kf**2 * n2**2)
        e = _ex_base_sr(n2, s2, z)
        if _K != 'K0':
            tu_ = 0.3 * (6 * np.pi**2) ** (2 / 3) * r ** (5 / 3) * ZREG
            zs_ = _zsat((s / (8 * r) + tu_) / (t + tu_))
            w = zs_**2 if _K == 'K3' else zeta**2 * zs_**2
            e = (1 - w) * e + w * _ex_branch(r, s)
        ex = ex + 0.5 * e
    # correlation: opus4 pair-aware Stoll PBE (prior art, unchanged)
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _pw92(rs, zeta)
    phi = 0.5 * ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3))
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi**2 * ks2 * n**2)
    phi3 = phi**3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    H = GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))
    corr = n * (ec + H)
    za_ = _zsat(saa / (8 * ra * ta)); zb_ = _zsat(sbb / (8 * rb * tb))
    for r, sg_s, t_s in ((ra, saa, ta), (rb, sbb, tb)):
        zs = _zsat(sg_s / (8 * r * t_s))
        zo = zb_ if r is ra else za_
        corr = corr - zs**2 * (1 - zo**PW * (1 - zeta**2)) * _ec_pol(r, sg_s)
    return ex * (1 - A_SR) + corr


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
        with np.errstate(all='ignore'):
            E0[mask] = np.real(_E(*sel))
            h = 1e-30
            for i in range(7):
                v = [s.astype(complex) for s in sel]
                v[i] = v[i] + 1j * h
                d[i, mask] = np.imag(_E(*v)) / h
    n = args[0] + args[1]
    exc = E0 / n
    if spin == 0:
        vxc = (0.5 * (d[0] + d[1]), 0.25 * (d[2] + d[3] + d[4]), None, 0.5 * (d[5] + d[6]))
    else:
        vxc = (np.stack([d[0], d[1]], 1), np.stack([d[2], d[3], d[4]], 1), None, np.stack([d[5], d[6]], 1))
    return exc, vxc, None, None
