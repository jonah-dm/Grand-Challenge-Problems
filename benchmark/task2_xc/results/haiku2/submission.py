"""haiku2 submission: spin-scaled meta-GGA with a one-orbital (iso-orbital indicator) exchange factor.

Exchange (spin-scaled): E_x = 1/2 [X(2n_up, 4 s_uu, 2 tau_up) + X(2n_dn, 4 s_dd, 2 tau_dn)], with
  X(n, s, tau) = eps_x^LDA(n) * F(p) * G(z) * n
  F(p) = 1 + kappa*(1 - exp(-mu*p/kappa)),  p = s/(4 kF^2 n^2)  (reduced gradient squared)
         mu = 10/81 (exact gradient-expansion coefficient), kappa = 0.804 (ceiling 1 + kappa).
  z    = tau_W / tau, tau_W = s/(8 n)  (iso-orbital indicator; z = 1 for a single orbital)
  G(z) = 1 - z/(1+z)  in [1/2, 1]; G = 1 in the uniform-gas limit; 1/2 at z = 1.
         Approximate one-orbital scaling (exact exchange there is -1/2 of the Hartree energy),
         not an exact pointwise identity.
Correlation: PW92 spin-interpolated LDA correlation plus the PBE-type gradient term with an
  exponential cutoff, H = gamma*phi^3 ln(1 + (beta/gamma) t^2 exp(-A t^2)); beta, gamma are the
  standard nonempirical constants. No dev-set fitting.

Derivatives are computed by complex-step differentiation of the pointwise energy density
(exact to machine precision).
"""
import numpy as np

HYB = 0.0
META = True

KAPPA = 0.804
MU = 10.0 / 81.0
BETA = 0.06672455060314922
GAMMA = (1.0 - np.log(2.0)) / np.pi ** 2
FPP0 = 1.709920934161365  # f''(0) of the PW92 spin interpolation


def _G(rs, A, a1, b1, b2, b3, b4):
    q = 2.0 * A * (b1 * rs ** 0.5 + b2 * rs + b3 * rs ** 1.5 + b4 * rs ** 2)
    return -2.0 * A * (1.0 + a1 * rs) * np.log(1.0 + 1.0 / q)


def _pw92(rs, zeta):
    ec0 = _G(rs, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
    ec1 = _G(rs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    alpha_c = -_G(rs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
    f = ((1.0 + zeta) ** (4.0 / 3.0) + (1.0 - zeta) ** (4.0 / 3.0) - 2.0) / (2.0 ** (4.0 / 3.0) - 2.0)
    z4 = zeta ** 4
    return ec0 + alpha_c * f * (1.0 - z4) / FPP0 + (ec1 - ec0) * f * z4


def _exchange_X(n, s, tau):
    """Per-volume exchange energy of a spin-resolved density, complex-safe."""
    kf = (3.0 * np.pi ** 2 * n) ** (1.0 / 3.0)
    ex_lda = -0.75 * (3.0 / np.pi) ** (1.0 / 3.0) * n ** (4.0 / 3.0)
    p = s / (4.0 * kf ** 2 * n ** 2)
    fx = 1.0 + KAPPA * (1.0 - np.exp(-MU * p / KAPPA))
    z = (s / (8.0 * n)) / tau
    gz = 1.0 - z / (1.0 + z)
    return ex_lda * fx * gz


def _E_pol(n1, n2, suu, sud, sdd, t1, t2):
    """Total exchange-correlation energy per volume for a spin-resolved density."""
    n = n1 + n2
    s_tot = suu + 2.0 * sud + sdd
    ex = 0.5 * (_exchange_X(2.0 * n1, 4.0 * suu, 2.0 * t1) + _exchange_X(2.0 * n2, 4.0 * sdd, 2.0 * t2))
    zeta = (n1 - n2) / n
    rs = (3.0 / (4.0 * np.pi * n)) ** (1.0 / 3.0)
    ec = _pw92(rs, zeta)
    phi = 0.5 * ((1.0 + zeta) ** (2.0 / 3.0) + (1.0 - zeta) ** (2.0 / 3.0))
    kf = (3.0 * np.pi ** 2 * n) ** (1.0 / 3.0)
    ks2 = 4.0 * kf / np.pi
    t2_ = s_tot / (4.0 * phi ** 2 * ks2 * n ** 2)
    g3 = GAMMA * phi ** 3
    A = (BETA / GAMMA) / (np.exp(-ec / g3) - 1.0)
    H = g3 * np.log(1.0 + (BETA / GAMMA) * t2_ * np.exp(-A * t2_))
    return ex + n * (ec + H)


_ARGS_IDX = range(7)  # n1, n2, suu, sud, sdd, t1, t2


def _derivs(args, h=1e-20):
    out = []
    for i in _ARGS_IDX:
        a = list(args)
        a[i] = a[i] + 1j * h
        out.append(np.imag(_E_pol(*a)) / h)
    return out


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None, **kwargs):
    rho = np.asarray(rho, dtype=float)
    floor = 1e-12
    if spin == 0 or rho.ndim == 2:
        n = np.maximum(rho[0], floor)
        grad = rho[1:4]
        s = np.maximum(grad[0] ** 2 + grad[1] ** 2 + grad[2] ** 2, 0.0)
        tau = np.maximum(rho[5 if rho.shape[0] >= 6 else 4], floor)
        n1 = n2 = 0.5 * n
        suu = sdd = sud = 0.25 * s
        t1 = t2 = 0.5 * tau
        args = (n1, n2, suu, sud, sdd, t1, t2)
        E = _E_pol(*args)
        d = _derivs(args)
        eps = E / n
        vrho = 0.5 * (d[0] + d[1])
        vsigma = 0.25 * (d[2] + d[3] + d[4])
        vtau = 0.5 * (d[5] + d[6])
        return eps, (vrho, vsigma, None, vtau), None, None

    # spin-polarized: rho shape (2, nvar, N)
    n1 = np.maximum(rho[0, 0], floor)
    n2 = np.maximum(rho[1, 0], floor)
    g1 = rho[0, 1:4]
    g2 = rho[1, 1:4]
    suu = np.maximum(np.sum(g1 * g1, axis=0), 0.0)
    sdd = np.maximum(np.sum(g2 * g2, axis=0), 0.0)
    sud = np.sum(g1 * g2, axis=0)
    ti = 5 if rho.shape[1] >= 6 else 4
    t1 = np.maximum(rho[0, ti], floor)
    t2 = np.maximum(rho[1, ti], floor)
    args = (n1, n2, suu, sud, sdd, t1, t2)
    n = n1 + n2
    E = _E_pol(*args)
    d = _derivs(args)
    eps = E / n
    # PySCF spin convention: each returned block is [grid, component]
    vrho = np.stack([d[0], d[1]], axis=-1)
    vsigma = np.stack([d[2], d[3], d[4]], axis=-1)
    vtau = np.stack([d[5], d[6]], axis=-1)
    return eps.reshape(-1, 1), (vrho, vsigma, None, vtau), None, None
