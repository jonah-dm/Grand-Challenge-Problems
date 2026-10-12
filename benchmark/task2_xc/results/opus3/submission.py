"""opus3: one-electron-consistent hybrid meta-GGA (OEC-h).  See NOTES.md.

E_xc = a E_x^HF + (1-a) E_x^SL[n_a, n_b] + E_c^G
 a = 1/4
 E_x^SL : spin-scaled, Fx(s,z) = 1 + kappa - kappa/(1 + mu(z) s^2/kappa)
          mu(z) = mu0 + (mu1 - mu0) w(z),  z = tau_W/tau (per spin channel)
          mu0 = beta*pi^2/3 (cancels the PBE-correlation gradient term in the slowly varying limit)
          mu1 = derived so that E_x^SL[H atom] = -5/16 Ha exactly (one-electron norm)
          kappa = 0.804 (Lieb-Oxford)
 E_c^G  : PBE correlation x (1 - zeta^2 z^2): vanishes for every fully polarized one-orbital density.
"""
import numpy as np
from pyscf.dft import libxc

HYB = 0.25
META = True
WPOW = 2
KAPPA = 0.804
BETA = 0.06672455060314922
MU0 = BETA * np.pi**2 / 3
MU1 = None  # set below by the hydrogen-atom norm
CX = -0.75 * (3 / np.pi) ** (1 / 3)
TINY = 1e-14


def _fx(s2, z, mu1):
    w = z**WPOW
    mu = MU0 + (mu1 - MU0) * w
    return 1 + KAPPA - KAPPA / (1 + mu * s2 / KAPPA)


def _zclip(z):
    # smooth-enough clip to [0,1] that is complex-step safe
    zr = np.real(z)
    return np.where(zr < 1.0, np.where(zr > 0.0, z, 0.0 * z), 1.0 + 0.0 * z)


def _ex_unpol(n, sig, tau, mu1):
    """energy density (per volume) of spin-unpolarized semilocal exchange."""
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    s2 = sig / (4 * kf**2 * n**2)
    z = _zclip(sig / (8 * n * tau))
    return CX * n ** (4 / 3) * _fx(s2, z, mu1)


def _hydrogen_mu1():
    from scipy.integrate import quad
    from scipy.optimize import brentq

    def ex(mu1):
        def f(r):
            n2 = 2 * np.exp(-2 * r) / np.pi
            kf = (3 * np.pi**2 * n2) ** (1 / 3)
            s2 = (2 * n2) ** 2 / (4 * kf**2 * n2**2)
            return 0.5 * CX * n2 ** (4 / 3) * _fx(s2, 1.0, mu1) * 4 * np.pi * r * r
        return quad(f, 1e-12, 60, limit=400)[0]
    return brentq(lambda m: ex(m) + 5 / 16, 0.1, 2.0)


MU1 = _hydrogen_mu1()


def _cs_grad(f, args, h=1e-30):
    """value and complex-step gradient of f wrt each arg."""
    val = np.real(f(*args))
    grads = []
    for i in range(len(args)):
        a = [x.astype(complex) for x in args]
        hi = h * np.maximum(np.abs(args[i]), 1e-30)
        a[i] = a[i] + 1j * hi
        grads.append(np.imag(f(*a)) / hi)
    return val, grads


def _xchan(n, sig, tau):
    # exchange of one spin channel: (1/2) E_x^unpol[2n]
    f = lambda n_, s_, t_: 0.5 * _ex_unpol(2 * n_, 4 * s_, 2 * t_, MU1)
    return _cs_grad(f, [n, sig, tau])


def _gate(na, nb, sig, tau):
    n = na + nb
    zeta = (na - nb) / n
    z = _zclip(sig / (8 * n * tau))
    return 1 - zeta**2 * z**2


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    rho = np.asarray(rho)
    if spin == 0:
        n = np.maximum(rho[0], TINY)
        sig = np.maximum(np.sum(rho[1:4] ** 2, axis=0), 1e-30)
        tau = np.maximum(rho[-1], sig / (8 * n) * (1 + 1e-12) + 1e-30)
        f = lambda n_, s_, t_: _ex_unpol(n_, s_, t_, MU1)
        ex, (vn, vs, vt) = _cs_grad(f, [n, sig, tau])
        ex *= (1 - HYB); vn *= (1 - HYB); vs *= (1 - HYB); vt *= (1 - HYB)
        ec, vc = libxc.eval_xc(",PBE", rho[:4], 0, 0, 1)[:2]
        exc = ex / n + ec
        vrho = vn + vc[0]
        vsigma = vs + vc[1]
        vtau = vt
        exc = np.where(rho[0] > TINY, exc, 0.0)
        return exc, (vrho, vsigma, np.zeros_like(vrho), vtau), None, None
    ra, rb = rho[0], rho[1]
    na = np.maximum(ra[0], TINY); nb = np.maximum(rb[0], TINY)
    saa = np.maximum(np.sum(ra[1:4] ** 2, axis=0), 1e-30)
    sbb = np.maximum(np.sum(rb[1:4] ** 2, axis=0), 1e-30)
    sab = np.sum(ra[1:4] * rb[1:4], axis=0)
    ta = np.maximum(ra[-1], saa / (8 * na) * (1 + 1e-12) + 1e-30)
    tb = np.maximum(rb[-1], sbb / (8 * nb) * (1 + 1e-12) + 1e-30)
    exa, (vna, vsa, vta) = _xchan(na, saa, ta)
    exb, (vnb, vsb, vtb) = _xchan(nb, sbb, tb)
    sc = 1 - HYB
    n = na + nb
    ec, vc = libxc.eval_xc(",PBE", np.array([ra[:4], rb[:4]]), 1, 0, 1)[:2]
    vcr, vcs = np.asarray(vc[0]), np.asarray(vc[1])
    stot = np.maximum(saa + 2 * sab + sbb, 1e-30)
    G, (gna, gnb, gs, gt) = _cs_grad(_gate, [na, nb, stot, ta + tb])
    Ec = n * ec
    exc = (sc * (exa + exb) + Ec * G) / n
    vrho = np.zeros((len(n), 2)); vsig = np.zeros((len(n), 3)); vtau = np.zeros((len(n), 2))
    vrho[:, 0] = sc * vna + vcr[:, 0] * G + Ec * gna
    vrho[:, 1] = sc * vnb + vcr[:, 1] * G + Ec * gnb
    vsig[:, 0] = sc * vsa + vcs[:, 0] * G + Ec * gs
    vsig[:, 1] = vcs[:, 1] * G + 2 * Ec * gs
    vsig[:, 2] = sc * vsb + vcs[:, 2] * G + Ec * gs
    vtau[:, 0] = sc * vta + Ec * gt
    vtau[:, 1] = sc * vtb + Ec * gt
    exc = np.where(n > 2 * TINY, exc, 0.0)
    return exc, (vrho, vsig, np.zeros((len(n), 2)), vtau), None, None
