"""sonnet3: spin-resolved one-orbital-gated correlation + one-orbital-aware exchange, 25% EXX.

NEW CONSTRAINT (exact): in any spin channel that is described by a single orbital
(z_sigma = tau_W,sigma/tau_sigma = 1) there is no same-spin correlation (one electron per spin
orbital cannot correlate with itself or with a same-spin partner). So
   E_c = E_c^PBE[na,nb] - sum_sigma z_sigma^2 * E_c^PBE[n_sigma,0]   (same-spin Stoll part gated)
vanishes for H, and for every fully polarised one-orbital density, AND the same-spin
correlation vanishes for He/H2-like two-electron singlets, which the zeta^2 z^2 gate of
sonnet2 cannot do. Exchange: z-blend of GE2 (mu=10/81) and s^-1/2-tail one-orbital form (see below).
"""
import os
import numpy as np
from scipy.special import erf as _erf
OMEGA = float(os.environ.get('O4_OMEGA', '0'))
_MUMODE=os.environ.get('S3_MU','pbe'); _GATE=os.environ.get('S3_GATE','pair'); PW=float(os.environ.get('O4_PW','2'))

HYB = 0.25
META = True
if OMEGA > 0:
    # long-range corrected: 25% SR exact exchange, 100% LR exact exchange
    RSH = (OMEGA, 1.0, -0.75)


def _att(a):
    """LDA short-range attenuation factor for erfc(omega r)/r (complex-step safe)"""
    big = np.real(a) > 3.0
    ab = np.where(big, 1.0 + 0 * a, a)
    asm = np.where(big, a, 1.0 + 0 * a)
    Fd = 1 - 8 / 3 * ab * (np.sqrt(np.pi) * _erf(1 / (2 * ab)) + (2 * ab - 4 * ab**3) * np.exp(-1 / (4 * ab**2)) - 3 * ab + 4 * ab**3)
    x = 1 / asm**2
    Fs = x / 36 - x**2 / 960 + x**3 / 26880 - x**4 / 829440
    return np.where(big, Fs, Fd)
MU_GE = 10.0 / 81.0
KAPPA, MU, BETA = 0.804, 0.2195149727645171, 0.06672455060314922
CTAIL = 4.0 / (5.0 * KAPPA / MU) ** 1.25   # sets max(F_x-1)=kappa exactly (LO)
GAMMA = (1 - np.log(2)) / np.pi**2
EPS = 1e-14


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


def _ex_unpol_density(n2, s2, z):
    """exchange energy density (per volume) of an unpolarized gas of density n2, s^2=s2"""
    ex_lda = -0.75 * (3 / np.pi) ** (1 / 3) * n2 ** (4 / 3)
    # (i) many-orbital / slowly varying regime (z->0): PBE-form, exact GE2 coefficient mu=10/81
    _m = MU_GE if _MUMODE == 'ge2' else MU
    h_gga = KAPPA * (1 - 1 / (1 + _m * s2 / KAPPA))
    # (ii) one-orbital regime (z->1), e.g. density tails: F_x ~ s^(-1/2) at large s (exact tail
    #      behaviour of the exchange energy density), mu = mu_PBE at small s, max(F_x)-1 = kappa (LO)
    s5 = s2 * np.sqrt(s2)
    h_orb = MU * s2 / (1 + CTAIL * s5)
    w = z**2
    h = (1 - w) * h_gga + w * h_orb
    if OMEGA > 0:
        # ITYH: short-range part with effective Fermi wavevector kF/sqrt(Fx)
        kf = (3 * np.pi**2 * n2) ** (1 / 3)
        a = OMEGA * np.sqrt(1 + h) / (2 * kf)
        return ex_lda * (1 + h) * _att(a)
    return ex_lda * (1 + h)


def _zsat(z, k=40.0):
    """smooth saturation z -> min(z,1) (z>1 is unphysical); analytic, so SCF stays smooth"""
    return 1 - np.log(1 + np.exp(k * (1 - z))) / k


def _E(ra, rb, saa, sab, sbb, ta, tb):
    """semilocal xc energy per volume (complex-step safe); exact-exchange part removed."""
    n = ra + rb
    sig = saa + 2 * sab + sbb
    z = sig / (8 * n * (ta + tb))
    z = _zsat(z)
    # exchange (spin scaling)
    ex = 0.0
    for r, s in ((ra, saa), (rb, sbb)):
        n2 = 2 * r
        kf = (3 * np.pi**2 * n2) ** (1 / 3)
        s2 = 4 * s / (4 * kf**2 * n2**2)       # sigma_unpol = 4 sigma_s
        ex = ex + 0.5 * _ex_unpol_density(n2, s2, z)
    # correlation
    zeta = (ra - rb) / n
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _pw92(rs, zeta)
    phi = 0.5 * ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3))
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi**2 * ks2 * n**2)
    phi3 = phi**3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    H = GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))
    # same-spin correlation of a one-orbital spin channel is zero (Stoll partition)
    tot = n * (ec + H)
    corr = tot
    za_ = _zsat(saa / (8 * ra * ta)); zb_ = _zsat(sbb / (8 * rb * tb))
    for r, sg_s, t_s in ((ra, saa, ta), (rb, sbb, tb)):
        zs = _zsat(sg_s / (8 * r * t_s))
        if _GATE == 'stoll':
            corr = corr - zs**2 * _ec_pol(r, sg_s)
        elif _GATE == 'pair':
            # remove same-spin part only where the OTHER channel is not also one-orbital;
            # in a two-electron singlet region it is reassigned to opposite-spin correlation
            zo = zb_ if r is ra else za_
            corr = corr - zs**2 * (1 - zo**PW) * _ec_pol(r, sg_s)
    return ex * (1 - HYB) + corr


def _ec_pol(r, sg):
    """PBE correlation energy per volume of a fully polarised density (r, sigma=sg)"""
    rs = (3 / (4 * np.pi * r)) ** (1 / 3)
    ec = _pw92(rs, 1.0 + 0 * r)
    phi = 2 ** (-1 / 3)
    kf = (3 * np.pi**2 * r) ** (1 / 3)
    t2 = sg / (4 * phi**2 * (4 * kf / np.pi) * r**2)
    phi3 = phi**3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    H = GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))
    return r * (ec + H)


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
    args = [ra, rb, saa, sab, sbb, ta, tb]
    args = [np.asarray(a, dtype=float) for a in args]
    mask = (ra + rb) > 1e-11
    # evaluate only on mask, fill rest with zero
    sel = [a[mask] for a in args]
    E0 = np.zeros(N); d = np.zeros((7, N))
    if mask.any():
        E0[mask] = np.real(_E(*sel))
        h = 1e-30
        for i in range(7):
            v = [s.astype(complex) for s in sel]
            v[i] = v[i] + 1j * h
            d[i, mask] = np.imag(_E(*v)) / h
    n = args[0] + args[1]
    exc = E0 / n
    if spin == 0:
        vrho = 0.5 * (d[0] + d[1])
        vsigma = 0.25 * (d[2] + d[3] + d[4])
        vtau = 0.5 * (d[5] + d[6])
        vxc = (vrho, vsigma, np.zeros(N), vtau)
    else:
        vrho = np.stack([d[0], d[1]], 1)
        vsigma = np.stack([d[2], d[3], d[4]], 1)
        vtau = np.stack([d[5], d[6]], 1)
        vxc = (vrho, vsigma, np.zeros((N, 2)), vtau)
    return exc, vxc, None, None
