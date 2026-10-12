"""haiku3 submission: 'HZX' meta-GGA exchange (one-orbital-interpolated) + PBE correlation.

Exchange (spin-scaled, per spin channel n = 2 rho_s):
    e_x = -C_x n^{1/3} * n * F_x(s, z)
    F_x = F0(s) + z * (F1(s) - F0(s)),   z = tau_W / tau  (iso-orbital indicator, per spin)
    F_k(s) = 1 + kappa_k - kappa_k / (1 + mu s^2 / kappa_k),   mu = 10/81 (2nd-order gradient
             expansion coefficient of exchange), kappa_k bounded by Lieb-Oxford (1 + kappa <= 1.804)
    kappa0 = 0.804 (Lieb-Oxford saturation, z -> 0: slowly varying / many-orbital regions)
    kappa1 = kappa0 = 0.804 (LO ceiling 1 + kappa <= 1.804 saturated at large s).
    mu1    : single constant fixed by the EXACT exchange energy of the hydrogen atom,
             E_x(H) = -5/16 Ha (one-orbital, z = 1 exactly). Solved once at import, not fitted.
             It is larger than the GE value 10/81 because one-orbital densities are far from the
             gradient-expansion regime; the GE is kept for z -> 0 only.
Correlation: PBE correlation taken unchanged from libxc (',PBE') as a building block.

Derivatives: exchange derivatives by complex-step differentiation (exact to round-off);
correlation derivatives from libxc.
"""
import numpy as np
from pyscf.dft import libxc
from scipy.optimize import brentq

HYB = 0.0
META = True

MU_GE = 10.0 / 81.0
KAPPA0 = 0.804
EPS_N = 1e-16
EPS_T = 1e-16
EPS_S = 1e-24
C_X = 0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
H_EXCHANGE_EXACT = -5.0 / 16.0


def _fx_vol(n, sig, tau, kappa1, mu1):
    """Volume exchange energy density of one spin channel written with n = 2 rho_s."""
    kf2 = (3.0 * np.pi ** 2 * n) ** (2.0 / 3.0)
    s2 = sig / (4.0 * kf2 * n * n)
    tw = sig / (8.0 * n)
    z = tw / tau
    F0 = 1.0 + KAPPA0 - KAPPA0 / (1.0 + MU_GE * s2 / KAPPA0)
    F1 = 1.0 + kappa1 - kappa1 / (1.0 + mu1 * s2 / kappa1)
    F = F0 + z * (F1 - F0)
    return -C_X * n ** (4.0 / 3.0) * F


def _ex_pol(ra, rb, saa, sbb, ta, tb, kappa1, mu1=None):  # kappa1 slot unused for MU1 (kept for signature)
    """Spin-polarized exchange energy density per volume (complex-safe, no abs/max)."""
    m1 = MU1 if mu1 is None else mu1
    return 0.5 * (_fx_vol(2.0 * ra + EPS_N, 4.0 * saa + EPS_S, 2.0 * ta + EPS_T, kappa1, m1)
                  + _fx_vol(2.0 * rb + EPS_N, 4.0 * sbb + EPS_S, 2.0 * tb + EPS_T, kappa1, m1))


def _h_atom_exchange(kappa1, mu1):
    r = np.linspace(1e-7, 40.0, 400001)
    rho = np.exp(-2.0 * r) / np.pi
    zero = np.zeros_like(r)
    ex = _ex_pol(rho, zero, 4.0 * rho ** 2, zero, rho / 2.0, zero, kappa1, mu1)
    f = 4.0 * np.pi * r ** 2 * ex
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(r)))   # trapezoid


def _solve_mu1():
    g = lambda m: _h_atom_exchange(KAPPA0, m) - H_EXCHANGE_EXACT
    return brentq(g, 0.05, 2.0, xtol=1e-12)


MU1 = _solve_mu1()
KAPPA1 = KAPPA0


def _exchange_with_derivs(ra, rb, saa, sbb, ta, tb, kappa1):
    base = [ra, rb, saa, sbb, ta, tb]
    f0 = np.real(_ex_pol(*base, kappa1))
    h = 1e-30
    grads = []
    for i in range(6):
        args = [np.asarray(b, dtype=complex) for b in base]
        args[i] = args[i] + 1j * h
        grads.append(np.imag(_ex_pol(*args, kappa1)) / h)
    return f0, grads


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    kappa1 = KAPPA1
    if spin == 0:
        rho = np.asarray(rho)
        n = rho[0]
        g = rho[1:4]
        sig = np.sum(g * g, axis=0)
        tau = rho[-1]
        ra = rb = n / 2.0
        saa = sbb = sig / 4.0
        ta = tb = tau / 2.0
        f0, gr = _exchange_with_derivs(ra, rb, saa, sbb, ta, tb, kappa1)
        ex_c, vc = libxc.eval_xc(",PBE", rho[:4], 0, relativity, 1, omega, verbose)[:2]
        nsafe = np.maximum(n, 1e-30)
        exc = (f0 + n * np.asarray(ex_c)) / nsafe
        vrho = 0.5 * (gr[0] + gr[1]) + np.asarray(vc[0])
        vsigma = 0.25 * (gr[2] + gr[3]) + np.asarray(vc[1])
        vtau = 0.5 * (gr[4] + gr[5])
        vlapl = np.zeros_like(vrho)
        return exc, (vrho, vsigma, vlapl, vtau), None, None
    rho = np.asarray(rho)
    ra_, rb_ = rho[0, 0], rho[1, 0]
    ga, gb = rho[0, 1:4], rho[1, 1:4]
    saa = np.sum(ga * ga, axis=0)
    sbb = np.sum(gb * gb, axis=0)
    sab = np.sum(ga * gb, axis=0)
    ta, tb = rho[0, -1], rho[1, -1]
    f0, gr = _exchange_with_derivs(ra_, rb_, saa, sbb, ta, tb, kappa1)
    ex_c, vc = libxc.eval_xc(",PBE", rho[:, :4], 1, relativity, 1, omega, verbose)[:2]
    ntot = ra_ + rb_
    nsafe = np.maximum(ntot, 1e-30)
    exc = (f0 + ntot * np.asarray(ex_c)) / nsafe
    vrho = np.stack([gr[0] + np.asarray(vc[0])[:, 0], gr[1] + np.asarray(vc[0])[:, 1]], axis=-1)
    vc1 = np.asarray(vc[1])
    vsigma = np.stack([gr[2] + vc1[:, 0], vc1[:, 1], gr[3] + vc1[:, 2]], axis=-1)
    vtau = np.stack([gr[4], gr[5]], axis=-1)
    vlapl = np.zeros_like(vrho)
    return exc, (vrho, vsigma, vlapl, vtau), None, None
