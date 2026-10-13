"""haiku4 candidate B: Stoll same-spin gate + one-orbital exchange normed to the exact H atom.

Built from two routes already in prior art, recombined with one new constraint-derived constant:
  - Correlation (sonnet3 route): E_c = E_c^PBE[n_a,n_b] - sum_s z_s^2 E_c^PBE[n_s,0].
    A spin channel holding one orbital (z_s = 1) has no same-spin correlation, so the gate
    removes the same-spin PBE piece exactly for every one-orbital channel. For H, E_c = 0 exactly.
  - Exchange: z-blend (weight z^2) of a PBE-form branch (z -> 0, mu = PBE value beta*pi^2/3, which
    cancels the PBE correlation gradient term) and a one-orbital branch of the same form,
        F_orb(s) - 1 = kappa (1 - 1/(1 + mu1 s^2 / kappa)),
    which saturates the Lieb-Oxford ceiling F_x <= 1 + kappa at large s.
    mu1 is NOT fitted: it is the unique root of E_x^SL[H atom] = -5/16 Ha (exact one-electron exchange),
    solved once at import. For H, z = 1 exactly, so the one-orbital branch is the one that enters.
  - 25% exact exchange (adiabatic-connection 1/4), semilocal exchange scaled by 0.75.

Derivatives: complex-step differentiation for exchange, correlation and the gate (exact to round-off).
"""
import numpy as np
from scipy.optimize import brentq

HYB = 0.25
META = True
ORB = "hnorm"          # one-orbital branch: 'hnorm' (this candidate) or 'tail' (sonnet3 form, ablation)

KAPPA = 0.804                      # Lieb-Oxford ceiling 1 + kappa <= 1.804
MU_GGA = 0.2195149727645171        # beta*pi^2/3 with beta = 0.06672455 (PBE), z -> 0 branch
BETA = 0.06672455060314922
GAMMA = (1 - np.log(2)) / np.pi ** 2
C_X = 0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
CTAIL = 4.0 / (5.0 * KAPPA / MU_GGA) ** 1.25   # only used by ORB == 'tail'
H_EXCHANGE_EXACT = -5.0 / 16.0
EPS = 1e-14


def _h_orb_hnorm(s2, mu1):
    return KAPPA * (1.0 - 1.0 / (1.0 + mu1 * s2 / KAPPA))


def _h_orb_tail(s2, mu1=None):
    s5 = s2 * np.sqrt(s2)
    return MU_GGA * s2 / (1.0 + CTAIL * s5)


def _h_gga(s2):
    return KAPPA * (1.0 - 1.0 / (1.0 + MU_GGA * s2 / KAPPA))


def _ex_unpol(n2, s2, z, mu1):
    """Exchange energy per volume of an unpolarized gas of density n2 (spin-scaled use below)."""
    ex_lda = -C_X * n2 ** (4.0 / 3.0)
    h_orb = _h_orb_hnorm(s2, mu1) if ORB == "hnorm" else _h_orb_tail(s2)
    w = z * z
    return ex_lda * (1.0 + (1.0 - w) * _h_gga(s2) + w * h_orb)


def _zsat(z, k=40.0):
    return 1 - np.log(1 + np.exp(k * (1 - z))) / k


def _ex_total(ra, rb, saa, sbb, sab, ta, tb, mu1):
    n = ra + rb
    sig = saa + 2 * sab + sbb
    z = _zsat(sig / (8 * n * (ta + tb)))
    ex = 0.0
    for r, sg in ((ra, saa), (rb, sbb)):
        n2 = 2 * r
        kf2 = (3 * np.pi ** 2 * n2) ** (2.0 / 3.0)
        s2 = sg / (kf2 * n2 * n2)          # = |grad n2|^2 / (4 kf^2 n2^2), spin-scaled
        ex = ex + 0.5 * _ex_unpol(n2, s2, z, mu1)
    return ex


def _h_atom_exchange(mu1):
    r = np.linspace(1e-7, 40.0, 400001)
    rho = np.exp(-2.0 * r) / np.pi            # H 1s density, one spin channel occupied
    # grad rho = -2 rho r_hat, tau_W = |grad rho|^2/(8 rho) = tau: z = 1 exactly
    n2 = 2.0 * rho
    s2 = 4.0 * rho ** 2 / ((3 * np.pi ** 2 * n2) ** (2.0 / 3.0) * n2 * n2)
    ex = 0.5 * _ex_unpol(n2, s2, 1.0, mu1)
    f = 4.0 * np.pi * r ** 2 * ex
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(r)))


MU1 = brentq(lambda m: _h_atom_exchange(m) - H_EXCHANGE_EXACT, 0.05, 2.0, xtol=1e-12) \
    if ORB == "hnorm" else MU_GGA


# ---- PBE correlation (PW92 + PBE H), per volume, complex-step safe ----
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
    z4 = zeta ** 4
    return ec0 + mac * fz * (1 - z4) / fpp0 + (ec1 - ec0) * fz * z4


def _pbe_H(rs, zeta, phi, t2, ec):
    phi3 = phi ** 3
    A = BETA / GAMMA / (np.exp(-ec / (GAMMA * phi3)) - 1)
    return GAMMA * phi3 * np.log(1 + BETA / GAMMA * t2 * (1 + A * t2) / (1 + A * t2 + A * A * t2 * t2))


def _ec_pbe_total(ra, rb, saa, sab, sbb):
    n = ra + rb
    sig = saa + 2 * sab + sbb
    zeta = (ra - rb) / n
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _pw92(rs, zeta)
    phi = 0.5 * ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3))
    kf = (3 * np.pi ** 2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi ** 2 * ks2 * n ** 2)
    return n * (ec + _pbe_H(rs, zeta, phi, t2, ec))


def _ec_pol_single(r, sg):
    """PBE correlation per volume of a fully polarized channel (r, |grad r|^2 = sg)."""
    rs = (3 / (4 * np.pi * r)) ** (1 / 3)
    ec = _pw92(rs, 1.0 + 0 * r)
    phi = 2 ** (-1 / 3)
    kf = (3 * np.pi ** 2 * r) ** (1 / 3)
    t2 = sg / (4 * phi ** 2 * (4 * kf / np.pi) * r ** 2)
    return r * (ec + _pbe_H(rs, 1.0, phi, t2, ec))


def _E(ra, rb, saa, sab, sbb, ta, tb):
    """Semilocal xc energy per volume with the exact-exchange share removed (HYB handled by harness)."""
    ex = _ex_total(ra, rb, saa, sbb, sab, ta, tb, MU1)
    corr = _ec_pbe_total(ra, rb, saa, sab, sbb)
    for r, sg, t in ((ra, saa, ta), (rb, sbb, tb)):
        zs = _zsat(sg / (8 * r * t))
        corr = corr - zs * zs * _ec_pol_single(r, sg)
    return (1 - HYB) * ex + corr


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    rho = np.asarray(rho)
    EPS_R = 1e-14
    if spin == 0:
        n = rho[0]; sg = np.sum(rho[1:4] ** 2, axis=0); tau = rho[-1]
        ra = rb = n / 2; saa = sab = sbb = sg / 4; ta = tb = tau / 2
    else:
        ra, rb = rho[0][0], rho[1][0]
        ga, gb = rho[0][1:4], rho[1][1:4]
        saa, sbb, sab = np.sum(ga ** 2, 0), np.sum(gb ** 2, 0), np.sum(ga * gb, 0)
        ta, tb = rho[0][-1], rho[1][-1]
    N = len(ra)
    ra = np.maximum(ra, EPS_R); rb = np.maximum(rb, EPS_R)
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
