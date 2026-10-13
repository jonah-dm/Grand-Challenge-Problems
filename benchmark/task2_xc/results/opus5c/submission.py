"""opus5c: hybrid meta-GGA with a potential-exact one-orbital exchange branch.

E_xc = a E_x^HF + (1-a) sum_s int e_x^LDA(n_s) [(1-L_s) F0(s_s) + L_s F_H(s_s)] + E_c
  a = 1/4 (PBE0 adiabatic-connection argument; known)
  F0 = PBE-form many-orbital enhancement (kappa = 0.804 Lieb-Oxford; mu = beta pi^2/3; known)
  F_H = one-orbital branch whose exchange potential on the hydrogen density equals the exact
        exchange potential -v_H up to the constant c+ = 0.1895 (NOTES.md sec. 6b; tabulated from
        the ODE solution, see hbranch.py / make_table.py). Exact, self-consistent, for all hydrogenic ions.
  L_s = h(z_s) (W2: every one-orbital channel, exact spin scaling) or h(z_s) max(0,zeta_s)^2 (W1).
  h(z) = z^3 (4 - 3 z), z_s = tau_W,s / tau_s capped at 1.
  E_c = E_c^PBE[na,nb] - sum_s h(z_s)[1 - h(z_s')(1-zeta^2)] E_c^PBE[n_s,0]  (pair-aware Stoll
        same-spin removal; known B95/TPSS idea in the prior-art opus4 form).
Derivatives: complex-step differentiation (exact to machine precision).
Env switches exist only for the disclosed candidate comparison: O5_W in {W0,W1,W2}, O5_MU in {pbe,ge}.
"""
import os
import numpy as np
from scipy.interpolate import CubicSpline

HYB = 0.25
META = True
_W = os.environ.get("O5_W", "W2")
_KAPPA = 0.804
_BETA = 0.06672455060314922
_GAMMA = (1 - np.log(2)) / np.pi**2
_MU = _BETA * np.pi**2 / 3 if os.environ.get("O5_MU", "pbe") == "pbe" else 10 / 81
_CX = -(3 / 4) * (6 / np.pi) ** (1 / 3)
_CS2 = 4 * (6 * np.pi**2) ** (2 / 3)
_S0SQ = (6 * np.pi) ** (-2 / 3)
_TINY = 1e-14
_FH_T = np.array([0.0, 0.041673679, 0.08333402467, 0.1249943703, 0.166654716, 0.2083150617, 0.2499754074, 0.291635753, 0.3332960987, 0.3749564444, 0.41661679, 0.4582771357, 0.4999374814, 0.541597827, 0.5832581727, 0.6249185184, 0.666578864, 0.7082392097, 0.7498995554, 0.7915599011, 0.8332202467, 0.8748805924, 0.9165409381, 0.9582012837, 0.9998616294, 1.041521975, 1.083182321, 1.124842666, 1.166503012, 1.208163358, 1.249823703, 1.291484049, 1.333144395, 1.37480474, 1.416465086, 1.458125432, 1.499785777, 1.541446123, 1.583106469, 1.624766814, 1.66642716, 1.708087506, 1.749747851, 1.791408197, 1.833068543, 1.874728888, 1.916389234, 1.95804958, 1.999709925, 2.041370271, 2.083030617, 2.124690962, 2.166351308, 2.208011654, 2.249671999, 2.291332345, 2.332992691, 2.374653037, 2.416313382, 2.457973728, 2.499634074, 2.541294419, 2.582954765, 2.624615111, 2.666275456, 2.707935802, 2.749596148, 2.791256493, 2.832916839, 2.874577185, 2.91623753, 2.957897876, 2.999558222, 3.041218567, 3.082878913, 3.124539259, 3.166199604, 3.20785995, 3.249520296, 3.291180641, 3.332840987, 3.374501333, 3.416161678, 3.457822024, 3.49948237, 3.541142715, 3.582803061, 3.624463407, 3.666123752, 3.707784098, 3.749444444, 3.791104789, 3.832765135, 3.874425481, 3.916085826, 3.957746172, 3.999406518, 4.041066863, 4.082727209, 4.124387555, 4.1660479, 4.207708246, 4.249368592, 4.291028937, 4.332689283, 4.374349629, 4.416009974, 4.45767032, 4.499330666, 4.540991011, 4.582651357, 4.624311703, 4.665972048, 4.707632394, 4.74929274, 4.790953085, 4.832613431, 4.874273777, 4.915934122, 4.957594468, 4.999254814, 5.040915159, 5.082575505, 5.124235851, 5.165896196, 5.207556542, 5.249216888, 5.290877233, 5.332537579, 5.374197925, 5.41585827, 5.457518616, 5.499178962, 5.540839307, 5.582499653, 5.624159999, 5.665820344, 5.70748069, 5.749141036, 5.790801381, 5.832461727, 5.874122073, 5.915782418, 5.957442764, 5.99910311, 6.040763455, 6.082423801, 6.124084147, 6.165744492, 6.207404838, 6.249065184, 6.290725529, 6.332385875, 6.374046221, 6.415706566, 6.457366912, 6.499027258, 6.540687603, 6.582347949, 6.624008295, 6.66566864])
_FH_F = np.array([1.275573845, 1.275482828, 1.275214778, 1.274776669, 1.274175023, 1.273415928, 1.272505059, 1.271447696, 1.27024874, 1.268912733, 1.267443871, 1.265846019, 1.264122729, 1.262277254, 1.260312555, 1.258231324, 1.256035989, 1.253728728, 1.251311481, 1.248785962, 1.246153666, 1.243415884, 1.240573709, 1.237628047, 1.234579624, 1.231428998, 1.228176565, 1.224822568, 1.221367105, 1.217810134, 1.214151484, 1.210390859, 1.206527846, 1.202561923, 1.198492462, 1.194318737, 1.190039932, 1.185655141, 1.181163381, 1.176563592, 1.171854645, 1.167035347, 1.162104444, 1.15706063, 1.151902548, 1.146628797, 1.141237937, 1.135728493, 1.13009896, 1.124347807, 1.118473483, 1.112474418, 1.106349035, 1.100095744, 1.093712957, 1.087199084, 1.080552545, 1.073771767, 1.066855195, 1.059801292, 1.052608547, 1.045275475, 1.037800628, 1.030182596, 1.022420008, 1.014511547, 1.006455943, 0.9982519872, 0.9898985323, 0.9813944987, 0.9727388795, 0.9639307461, 0.954969253, 0.9458536438, 0.9365832562, 0.9271575277, 0.9175760017, 0.9078383327, 0.8979442928, 0.8878937773, 0.877686811, 0.8673235547, 0.856804311, 0.8461295314, 0.8352998227, 0.8243159534, 0.8131788612, 0.8018896594, 0.7904496443, 0.7788603022, 0.7671233172, 0.7552405778, 0.7432141852, 0.7310464608, 0.7187399536, 0.7062974483, 0.6937219734, 0.6810168087, 0.6681854939, 0.6552318366, 0.6421599201, 0.6289741122, 0.6156790732, 0.6022797639, 0.5887814541, 0.575189731, 0.5615105068, 0.5477500271, 0.533914879, 0.5200119986, 0.5060486787, 0.4920325764, 0.4779717204, 0.4638745175, 0.4497497598, 0.4356066301, 0.4214547086, 0.4073039776, 0.3931648264, 0.3790480553, 0.3649648791, 0.35092693, 0.3369462586, 0.3230353354, 0.3092070503, 0.2954747109, 0.2818520398, 0.2683531704, 0.2549926412, 0.2417853881, 0.2287467345, 0.2158923807, 0.2032383889, 0.1908011678, 0.178597453, 0.166644285, 0.1549589841, 0.1435591213, 0.1324624859, 0.1216870488, 0.1112509207, 0.1011723066, 0.09146945344, 0.08216059353, 0.07326388049, 0.06479731914, 0.05677868778, 0.04922545258, 0.04215467341, 0.03558290035, 0.02952606022, 0.02399933215, 0.01901701158, 0.01459236148, 0.01073744992, 0.007462972946, 0.004778061432, 0.002690070797, 0.001204352155, 0.0003240034547, 4.959903821e-05])

_CSP = CubicSpline(_FH_T, _FH_F, bc_type=((1, 0.0), (1, 0.0)))
_TEND = _FH_T[-1]


def _where(c, a, b):
    return np.where(c, a, b)


def _FH(s2):
    """One-orbital branch, complex-step aware: F(x + i h) = F(x) + i h F'(x)."""
    t = np.log(s2 / _S0SQ)
    tr, ti = t.real, t.imag
    trc = np.clip(tr, 0.0, _TEND)
    inside = (tr > 0.0) & (tr < _TEND)
    f = _CSP(trc)
    fp = _where(inside, _CSP(trc, 1), 0.0)
    return f + 1j * ti * fp


def _F0(s2):
    return 1 + _KAPPA - _KAPPA / (1 + _MU * s2 / _KAPPA)


def _h(z):
    zc = _where(z.real >= 1.0, 1.0 + 0j, z)
    zc = _where(zc.real <= 0.0, 0j, zc)
    return zc**3 * (4 - 3 * zc)


def _G(rs, A, a1, b1, b2, b3, b4):
    q = 2 * A * (b1 * np.sqrt(rs) + b2 * rs + b3 * rs**1.5 + b4 * rs**2)
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / q)


def _pw92(rs, zeta):
    ec0 = _G(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
    ec1 = _G(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    mac = _G(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
    fz = ((1 + zeta) ** (4 / 3) + (1 - zeta) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    z4 = zeta**4
    return ec0 - mac * fz * (1 - z4) / 1.709921 + (ec1 - ec0) * fz * z4


def _pbe_c(n, zeta, sig):
    """PBE correlation energy per volume (complex-safe). n > 0 assumed."""
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _pw92(rs, zeta)
    phi = ((1 + zeta) ** (2 / 3) + (1 - zeta) ** (2 / 3)) / 2
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi**2 * ks2 * n**2)
    phi3 = phi**3
    A = _BETA / _GAMMA / (np.exp(-ec / (_GAMMA * phi3)) - 1)
    At2 = A * t2
    H = _GAMMA * phi3 * np.log(1 + _BETA / _GAMMA * t2 * (1 + At2) / (1 + At2 + At2**2))
    return n * (ec + H)


def _pbe_c_pol(n, sig):
    """PBE correlation for a fully polarized density (zeta = 1 exactly)."""
    rs = (3 / (4 * np.pi * n)) ** (1 / 3)
    ec = _G(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    phi = 2 ** (2 / 3) / 2
    kf = (3 * np.pi**2 * n) ** (1 / 3)
    ks2 = 4 * kf / np.pi
    t2 = sig / (4 * phi**2 * ks2 * n**2)
    phi3 = phi**3
    A = _BETA / _GAMMA / (np.exp(-ec / (_GAMMA * phi3)) - 1)
    At2 = A * t2
    H = _GAMMA * phi3 * np.log(1 + _BETA / _GAMMA * t2 * (1 + At2) / (1 + At2 + At2**2))
    return n * (ec + H)


def _edens(na, nb, saa, sab, sbb, ta, tb):
    """xc energy per volume (semilocal part), complex-safe."""
    ma = na.real > _TINY
    mb = nb.real > _TINY
    na_ = _where(ma, na, 1.0)
    nb_ = _where(mb, nb, 1.0)
    saa_ = _where(ma, saa, 0.0)
    sbb_ = _where(mb, sbb, 0.0)
    ta_ = _where(ma & (ta.real > 1e-30), ta, 1e-30)
    tb_ = _where(mb & (tb.real > 1e-30), tb, 1e-30)
    za = _where(ma, saa_ / (8 * na_ * ta_), 0.0)
    zb = _where(mb, sbb_ / (8 * nb_ * tb_), 0.0)
    ha, hb = _h(za), _h(zb)
    n = na + nb
    mn = n.real > _TINY
    n_ = _where(mn, n, 1.0)
    zeta = (na - nb) / n_
    if _W == "W2":
        La, Lb = ha, hb
    elif _W == "W1":
        La = ha * _where(zeta.real > 0, zeta**2, 0.0)
        Lb = hb * _where(zeta.real < 0, zeta**2, 0.0)
    else:
        La = Lb = 0.0 * ha
    # exchange (spin channels)
    s2a = saa_ / (_CS2 * na_ ** (8 / 3))
    s2b = sbb_ / (_CS2 * nb_ ** (8 / 3))
    s2a = _where(ma, s2a, 1.0)
    s2b = _where(mb, s2b, 1.0)
    exa = _CX * na_ ** (4 / 3) * ((1 - La) * _F0(s2a) + La * _FH(s2a))
    exb = _CX * nb_ ** (4 / 3) * ((1 - Lb) * _F0(s2b) + Lb * _FH(s2b))
    ex = _where(ma, exa, 0.0) + _where(mb, exb, 0.0)
    # correlation
    eps = 1e-12
    zc = _where(zeta.real > 1 - eps, (1 - eps) + 0j, zeta)
    zc = _where(zc.real < -1 + eps, (-1 + eps) + 0j, zc)
    sig = saa + 2 * sab + sbb
    sig = _where(mn, sig, 0.0)
    ecf = _where(mn, _pbe_c(n_, zc, sig), 0.0)
    eca = _where(ma, _pbe_c_pol(na_, saa_), 0.0)
    ecb = _where(mb, _pbe_c_pol(nb_, sbb_), 0.0)
    pair = 1 - zeta**2
    Ga = ha * (1 - hb * pair)
    Gb = hb * (1 - ha * pair)
    ec = ecf - Ga * eca - Gb * ecb
    return (1 - HYB) * ex + ec


def _cstep(f, args, i):
    x = args[i]
    h = 1e-20 * np.maximum(np.abs(x), 1e-20)
    a = [np.asarray(v, dtype=complex) for v in args]
    a[i] = a[i] + 1j * h
    return f(*a).imag / h


def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, omega=None, verbose=None):
    if spin == 0:
        r = np.asarray(rho, dtype=float)
        n = r[0]
        sig = r[1] ** 2 + r[2] ** 2 + r[3] ** 2
        tau = r[-1]

        def f(n, s, t):
            return _edens(n / 2, n / 2, s / 4, s / 4, s / 4, t / 2, t / 2)

        args = [n, sig, tau]
        e = f(*[np.asarray(v, dtype=complex) for v in args]).real
        exc = np.where(n > _TINY, e / np.where(n > _TINY, n, 1.0), 0.0)
        vrho = _cstep(f, args, 0)
        vsig = _cstep(f, args, 1)
        vtau = _cstep(f, args, 2)
        return exc, (vrho, vsig, None, vtau), None, None
    ra, rb = np.asarray(rho[0], dtype=float), np.asarray(rho[1], dtype=float)
    na, nb = ra[0], rb[0]
    saa = ra[1] ** 2 + ra[2] ** 2 + ra[3] ** 2
    sbb = rb[1] ** 2 + rb[2] ** 2 + rb[3] ** 2
    sab = ra[1] * rb[1] + ra[2] * rb[2] + ra[3] * rb[3]
    ta, tb = ra[-1], rb[-1]
    args = [na, nb, saa, sab, sbb, ta, tb]
    e = _edens(*[np.asarray(v, dtype=complex) for v in args]).real
    n = na + nb
    exc = np.where(n > _TINY, e / np.where(n > _TINY, n, 1.0), 0.0)
    d = [_cstep(_edens, args, i) for i in range(7)]
    vrho = np.stack([d[0], d[1]], axis=1)
    vsig = np.stack([d[2], d[3], d[4]], axis=1)
    vtau = np.stack([d[5], d[6]], axis=1)
    return exc, (vrho, vsig, None, vtau), None, None
