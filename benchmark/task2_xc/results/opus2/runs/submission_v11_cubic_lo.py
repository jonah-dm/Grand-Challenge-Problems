"""opus2: "BIG" -- bounded-indicator, integrated-gradient-expansion meta-GGA exchange
(+ r2SCAN correlation as a building block).

E_xc = HYB*E_x^HF + (1-HYB)*E_x^BIG[rho, sigma, tau] + E_c^{r2SCAN}

Exchange enhancement factor (spin scaling applied exactly):
    Fx(p, alpha) = [ F0 + (F1(p) - F0) * W(alpha) ] * g(s)
    F0      = 1 + kappa0                        (iso-orbital, alpha = 0; kappa0 = 0.174 tight 2e bound)
    F1(p)   = 1 + kappa1 - kappa1/(1 + mu1 p/kappa1)  (kappa1 = 0.489 from LO, see below)
    W       = 2*beta,  beta = ah/(1+ah) in [0,1)   (W linear in a *bounded* indicator)
    ah      = alpha^2 (1+eps)/(alpha+eps)          (numerical smoothing at alpha=0 only)
    g(s)    = 1 - exp(-a1/sqrt(s))                 (Fx ~ s^-1/2, nonuniform scaling)
Constants: mu1 is NOT a free parameter: it is fixed so that the *integrated* second-order
gradient expansion coefficient is exactly 10/81 once the linear (alpha-1) slope of W is
transformed to p by integration by parts (int rho^{4/3}(alpha-1) = -(20/27) int rho^{4/3} p).
a1 is fixed by the exact exchange energy of the hydrogen atom (-5/16 Ha).
The alpha -> infinity limit Fx -> (2 F1 - F0) g follows from linearity in beta, with no
extra parameter. See NOTES.md.
"""
import numpy as np
from pyscf.dft import libxc

HYB = 0.0           # pure meta-GGA (decided a priori; see NOTES)
META = True
ETA = 0.1

MU = 10.0 / 81.0            # second-order gradient expansion of exchange
KAPPA_LO = 0.804            # Lieb-Oxford: pointwise Fx <= 1.804
KAPPA0 = 0.174              # tight two-electron / iso-orbital bound Fx <= 1.174
# alpha->inf limit Fx -> 2 F1 - F0 must obey LO for all s:  2(1+KAPPA1) - (1+KAPPA0) <= 1+KAPPA_LO
KAPPA1 = 0.5 * (KAPPA_LO + KAPPA0)   # = 0.489, the largest value allowed (saturated)
EPS = 1.0                   # ah = 2 alpha^2/(1+alpha): smooth (quadratic) at alpha=0, ah(1)=1
MU1 = MU   # W is flat (to 2nd order) at alpha=1, so the pointwise GE coefficient 10/81 is kept

_CX = -0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
_CF = 0.3 * (3.0 * np.pi ** 2) ** (2.0 / 3.0)


def _fx(p, alpha, a1):
    s = np.sqrt(p + 1e-30)
    f1 = 1.0 + KAPPA1 - KAPPA1 / (1.0 + MU1 * p / KAPPA1)
    f0 = 1.0 + KAPPA0 + 0.0 * p
    ah = alpha * alpha * (1.0 + EPS) / (alpha + EPS)   # smooth at alpha=0 (dah/dalpha=0), ah(1)=1
    t = (ah - 1.0) / (ah + 1.0)               # bounded indicator in [-1, 1): -1 iso-orbital, 0 UEG, ->1 alpha->inf
    w = 1.0 + t * t * t                       # W(0)=0, W(1)=1, W'(1)=W''(1)=0, W(inf)=2, monotone
    g = 1.0 - np.exp(-a1 / np.sqrt(s))
    return (f0 + (f1 - f0) * w) * g


def _ex_unpol(rho, sigma, tau, a1):
    """exchange energy per volume for a spin-unpolarized density (complex-step safe)."""
    rho = rho + 1e-30
    kf2 = (3.0 * np.pi ** 2 * rho) ** (2.0 / 3.0)
    p = sigma / (4.0 * kf2 * rho * rho)
    tw = sigma / (8.0 * rho)
    tu = _CF * rho ** (5.0 / 3.0)
    alpha = (tau - tw) / (tu + ETA * tw)   # regularized iso-orbital indicator (eta=0.1, chosen for SCF stability)
    alpha = np.where(np.real(alpha) < 0, 0.0 * alpha, alpha)
    return _CX * rho ** (4.0 / 3.0) * _fx(p, alpha, a1)


def _hydrogen_norm():
    # exact exchange of H atom density rho = exp(-2r)/pi: -5/16 Ha. Solve for a1 so the
    # semilocal (100%) exchange reproduces it (iso-orbital, alpha=0 everywhere).
    r = np.linspace(1e-5, 40, 200001)
    rho = np.exp(-2 * r) / np.pi
    sig = (2 * rho) ** 2
    tau = sig / (8 * rho)
    # spin-polarized: Ex[rho,0] = 1/2 Ex_unpol[2 rho]
    def ex(a1):
        e = 0.5 * _ex_unpol(2 * rho, 4 * sig, 2 * tau, a1)
        return np.trapezoid(e * 4 * np.pi * r * r, r)
    lo, hi = 0.5, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ex(mid) < -0.3125:   # too much exchange -> reduce a1
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)



A1 = 4.948698694818636   # hydrogen-atom exact-exchange norm (from _hydrogen_norm())
H = 1e-20


def _ex_and_derivs(rho, sigma, tau):
    """returns e (per volume), d/drho, d/dsigma, d/dtau via complex step."""
    e = np.real(_ex_unpol(rho, sigma, tau, A1))
    dr = np.imag(_ex_unpol(rho + 1j * H, sigma, tau, A1)) / H
    ds = np.imag(_ex_unpol(rho, sigma + 1j * H, tau, A1)) / H
    dt = np.imag(_ex_unpol(rho, sigma, tau + 1j * H, A1)) / H
    return e, dr, ds, dt


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    cx = 1.0 - HYB
    rho = np.asarray(rho, dtype=float)
    thr = 1e-14
    if spin == 0:
        r, gx, gy, gz, tau = rho[0], rho[1], rho[2], rho[3], rho[-1]
        sig = gx * gx + gy * gy + gz * gz
        m = r > thr
        e = np.zeros_like(r); vr = np.zeros_like(r); vs = np.zeros_like(r); vt = np.zeros_like(r)
        ee, a, b, c = _ex_and_derivs(r[m], sig[m], tau[m])
        e[m], vr[m], vs[m], vt[m] = ee, a, b, c
        exc_c, vxc_c = libxc.eval_xc("MGGA_C_R2SCAN", rho, 0, 0, 1)[:2]
        exc = cx * np.where(m, e / np.maximum(r, thr), 0.0) + exc_c
        vrho = cx * vr + vxc_c[0]
        vsig = cx * vs + vxc_c[1]
        vtau = cx * vt + vxc_c[3]
        return exc, (vrho, vsig, np.zeros_like(r), vtau), None, None
    ra, rb = rho[0], rho[1]
    rt = ra[0] + rb[0]
    out_e = np.zeros_like(rt)
    vrho = np.zeros((rt.size, 2)); vsig = np.zeros((rt.size, 3)); vtau = np.zeros((rt.size, 2))
    for k, rs in enumerate((ra, rb)):
        r, gx, gy, gz, tau = rs[0], rs[1], rs[2], rs[3], rs[-1]
        sig = gx * gx + gy * gy + gz * gz
        m = r > thr
        # Ex[ra,rb] = 1/2 sum_s Ex_unpol[2 rho_s]
        ee, a, b, c = _ex_and_derivs(2 * r[m], 4 * sig[m], 2 * tau[m])
        out_e[m] += 0.5 * ee
        vrho[m, k] = a
        vsig[m, 2 * k] = 2.0 * b
        vtau[m, k] = c
    exc_c, vxc_c = libxc.eval_xc("MGGA_C_R2SCAN", rho, 1, 0, 1)[:2]
    exc = cx * np.where(rt > thr, out_e / np.maximum(rt, thr), 0.0) + exc_c
    vrho = cx * vrho + vxc_c[0]
    vsig = cx * vsig + vxc_c[1]
    vtau = cx * vtau + vxc_c[3]
    return exc, (vrho, vsig, np.zeros((rt.size, 2)), vtau), None, None
