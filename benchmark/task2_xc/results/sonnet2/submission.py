"""SCF-SIC-free hybrid meta-GGA ("OSC-PBE0"): PBE-type exchange/correlation skeleton plus
a one-electron-spin-polarized self-correlation gate built from the iso-orbital indicator.

 e_c = e_c^PBE(n_a,n_b,grad) * [1 - zeta^2 * z^2],  z = tau_W/tau  (z=1 for any one-orbital density)
 -> e_c vanishes EXACTLY for every fully spin-polarized one-orbital density (H atom, any
    hydrogenic 1s etc.), while the uniform-gas and two-electron-singlet limits are untouched.
Exchange: spin-scaled PBE-form with kappa=0.804 (tight LO), mu = mu_PBE; 25% exact exchange
(a=1/n, n=4 from the MP4-based adiabatic connection argument). No parameter is fit to any data.
Derivatives by complex-step differentiation (exact to machine precision).
"""
import numpy as np

HYB = 0.25
META = True
KAPPA0 = 0.174
KAPPA, MU, BETA = 0.804, 0.2195149727645171, 0.06672455060314922
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
    # one-orbital (z->1) tightening: F_x <= 1.174 for any two-electron/one-orbital density
    # (exact LO-type bound, Perdew et al. SCAN constraint), smooth interpolation in z^2
    kap = KAPPA - (KAPPA - KAPPA0) * z**2
    h = kap * (1 - 1 / (1 + MU * s2 / kap))
    return ex_lda * (1 + h)


def _E(ra, rb, saa, sab, sbb, ta, tb):
    """semilocal xc energy per volume (complex-step safe); exact-exchange part removed."""
    n = ra + rb
    sig = saa + 2 * sab + sbb
    z = sig / (8 * n * (ta + tb))
    z = np.where(np.real(z) > 1, 1.0 + 0 * z, z)
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
    # one-orbital self-correlation gate
    gate = 1 - zeta**2 * z**2
    return ex * (1 - HYB) + n * (ec + H) * gate


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
