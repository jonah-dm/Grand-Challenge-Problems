"""opus6t: spin-blind one-orbital-manifold hybrid meta-GGA (round 6, theory track).
E_xc = 1/4 E_x^HF + int { (1-W)[3/4 e_x^PBE + e_c^pair] + W e_SB }          (MODE 'C', submitted)
  W    = zeta^2 z_tot^2,  z_tot = tau_W[n]/tau  (one-orbital manifold of the TOTAL density, NOTES sec. 4)
  e_SB = [ (1-a)(1+zeta^2)/2 + (1-zeta^2)/2 ] * e_x^q[n,0]   (exact C2 particle-triangle structure, NOTES sec. 3c)
  e_x^q[n,0]: fully polarised exchange of the total density with F_q = F_LO/(1+F_LO s^2/c_q), c_q from the H norm.
  e_c^pair: PBE correlation with opus4's pair-aware Stoll gate (prior art).
MODE 'A' = base only; MODE 'B' = opus5t-type per-channel C1 exchange (W_s = zeta^2 z_s^2), for comparison.
All derivatives by complex step on a complex-safe numpy implementation.
"""
import numpy as np

HYB = 0.25
META = True
MODE = 'B'
A_HF = 0.25
KAPPA = 0.804
MU = 0.2195149727645171
FLO = 1.0 + KAPPA
CQ = 2.208864      # H-atom norm, exact density (work/cq.py)
EPSZ = 1e-3          # z regulator (numerical only)
PI = np.pi
CX = -0.75 * (3.0 / PI) ** (1.0 / 3.0)

def _pw(rs, A, a1, b1, b2, b3, b4):
    sr = np.sqrt(rs)
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / (2 * A * (b1 * sr + b2 * rs + b3 * rs * sr + b4 * rs * rs)))

def _ec_unif(rs, zeta):
    ec0 = _pw(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
    ec1 = _pw(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    mac = _pw(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
    fz = ((1 + zeta) ** (4 / 3) + (1 - zeta) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    z4 = zeta ** 4
    return ec0 - mac * fz * (1 - z4) / 1.709921 + (ec1 - ec0) * fz * z4

def _ec_pbe(na, nb, g):
    """PBE correlation energy per volume; g = |grad n|^2. zeta clipped away from +-1 by floors."""
    n = na + nb
    zeta = (na - nb) / n
    rs = (3 / (4 * PI * n)) ** (1 / 3)
    ec = _ec_unif(rs, zeta)
    phi = ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3)) / 2
    kf = (3 * PI * PI * n) ** (1 / 3)
    ks = np.sqrt(4 * kf / PI)
    t2 = g / (4 * phi * phi * ks * ks * n * n)
    gam = (1 - np.log(2)) / PI ** 2
    beta = 0.06672455060314922
    A = beta / gam / (np.exp(-ec / (gam * phi ** 3)) - 1)
    at2 = A * t2
    H = gam * phi ** 3 * np.log(1 + beta / gam * t2 * (1 + at2) / (1 + at2 + at2 * at2))
    return n * (ec + H)

def _ex_ch(ns, gss, F):
    """spin-scaled exchange of one channel: 1/2 e_x^unif(2 ns) F(s), s from 2 ns."""
    n2 = 2 * ns
    kf = (3 * PI * PI * n2) ** (1 / 3)
    s2 = 4 * gss / (4 * kf * kf * n2 * n2)
    return 0.5 * CX * n2 ** (4 / 3) * F(s2)

F_pbe = lambda s2: 1 + KAPPA - KAPPA / (1 + MU * s2 / KAPPA)
F_q = lambda s2: FLO / (1 + FLO * s2 / CQ)

def _where(c, a, b):
    return np.where(c, a, b)

def _z(tw, t, tu):
    z = (tw + EPSZ * tu) / (t + EPSZ * tu)
    return _where(z.real > 1, 1.0 + 0 * z, z)

def edens(na, nb, gaa, gab, gbb, ta, tb, mode=None):
    mode = MODE if mode is None else mode
    tiny = 1e-14
    na = na + tiny; nb = nb + tiny
    gaa = gaa + 1e-30; gbb = gbb + 1e-30
    n = na + nb
    g = gaa + 2 * gab + gbb
    g = _where(g.real < 1e-30, 1e-30 + 0 * g, g)
    zeta = (na - nb) / n
    c6 = 0.3 * (6 * PI * PI) ** (2 / 3)
    za = _z(gaa / (8 * na), ta + 1e-30, c6 * na ** (5 / 3))
    zb = _z(gbb / (8 * nb), tb + 1e-30, c6 * nb ** (5 / 3))
    # correlation: opus4 pair-aware Stoll gate on PBE (prior art)
    ec = _ec_pbe(na, nb, g)
    eca = _ec_pbe(na, tiny + 0 * na, gaa)
    ecb = _ec_pbe(tiny + 0 * nb, nb, gbb)
    pair = 1 - zeta * zeta
    ec = ec - za * za * (1 - zb * zb * pair) * eca - zb * zb * (1 - za * za * pair) * ecb
    if mode == 'B':
        wa = zeta * zeta * za * za; wb = zeta * zeta * zb * zb
        ex = _ex_ch(na, gaa, lambda s2: (1 - wa) * F_pbe(s2) + wa * F_q(s2)) + \
             _ex_ch(nb, gbb, lambda s2: (1 - wb) * F_pbe(s2) + wb * F_q(s2))
        return (1 - A_HF) * ex + ec
    ex = _ex_ch(na, gaa, F_pbe) + _ex_ch(nb, gbb, F_pbe)
    base = (1 - A_HF) * ex + ec
    if mode == 'A':
        return base
    c3 = 0.3 * (3 * PI * PI) ** (2 / 3)
    zt = _z(g / (8 * n), ta + tb + 1e-30, c3 * n ** (5 / 3))
    W = zeta * zeta * zt * zt
    z2 = zeta * zeta
    eSB = ((1 - A_HF) * (1 + z2) / 2 + (1 - z2) / 2) * _ex_ch(n, g, F_q)
    return (1 - W) * base + W * eSB

H_CS = 1e-20

def _cs(f, args, i):
    a = [x.astype(complex) for x in args]
    h = H_CS * (np.abs(args[i]) + 1e-12)
    a[i] = a[i] + 1j * h
    return f(*a).imag / h

def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    rho = np.asarray(rho)
    if spin == 0:
        n = rho[0]; s = rho[1] ** 2 + rho[2] ** 2 + rho[3] ** 2; t = rho[-1]
        f = lambda n_, s_, t_: edens(n_ / 2, n_ / 2, s_ / 4, s_ / 4, s_ / 4, t_ / 2, t_ / 2)
        args = [n, s, t]
        e = edens(n / 2, n / 2, s / 4, s / 4, s / 4, t / 2, t / 2)
        mask = n > 1e-12
        exc = np.where(mask, e / np.maximum(n, 1e-300), 0.0)
        if deriv == 0:
            return exc, None, None, None
        vr, vs, vt = [np.where(mask, _cs(f, args, i), 0.0) for i in range(3)]
        return exc, (vr, vs, None, vt), None, None
    ra, rb = rho[0], rho[1]
    na, nb = ra[0], rb[0]
    gaa = ra[1] ** 2 + ra[2] ** 2 + ra[3] ** 2
    gbb = rb[1] ** 2 + rb[2] ** 2 + rb[3] ** 2
    gab = ra[1] * rb[1] + ra[2] * rb[2] + ra[3] * rb[3]
    ta, tb = ra[-1], rb[-1]
    args = [na, nb, gaa, gab, gbb, ta, tb]
    e = edens(*args)
    n = na + nb
    mask = n > 1e-12
    exc = np.where(mask, e / np.maximum(n, 1e-300), 0.0)
    if deriv == 0:
        return exc, None, None, None
    d = [np.where(mask, _cs(edens, args, i), 0.0) for i in range(7)]
    vrho = np.stack([d[0], d[1]], axis=1)
    vsig = np.stack([d[2], d[3], d[4]], axis=1)
    vtau = np.stack([d[5], d[6]], axis=1)
    return exc, (vrho, vsig, None, vtau), None, None
