"""V1: exact (N<=2 Fock space) spin-resolved Lieb functional for 1D soft-Coulomb one-orbital densities.
F[n_a,n_b] = sup_{v_a,v_b} { E_hull(v; f_a,f_b) - sum_s int v_s n_s }, computed via a smoothed (entropic, T)
grand-canonical dual, annealed T -> small. Corners: (0,0),(1,0),(0,1),(1,1) singlet-sector ground state,(2,0),(0,2).
usage: python lieb1d.py fa fb [shape]
"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.optimize import minimize
L, Ng = 10.0, 61
x = np.linspace(-L, L, Ng); h = x[1] - x[0]
w = lambda d: 1.0 / np.sqrt(d * d + 1.0)
T1 = sp.diags([np.full(Ng, 1.0 / h**2), np.full(Ng - 1, -0.5 / h**2), np.full(Ng - 1, -0.5 / h**2)], [0, 1, -1]).tocsr()
I = sp.identity(Ng, format='csr')
W2 = w(x[:, None] - x[None, :]).ravel()
# antisymmetric (2,0) basis i<j
iu, ju = np.triu_indices(Ng, 1)
P = sp.csr_matrix((np.r_[np.ones(len(iu)), -np.ones(len(iu))] / np.sqrt(2),
                   (np.r_[iu * Ng + ju, ju * Ng + iu], np.r_[np.arange(len(iu)), np.arange(len(iu))])), shape=(Ng * Ng, len(iu)))
TT = sp.kron(T1, I) + sp.kron(I, T1)

def gs1(v):
    H = (T1 + sp.diags(v)).toarray(); e, c = np.linalg.eigh(H); p = c[:, 0]**2 / h * 0 + c[:, 0]**2
    return e[0], p / h  # density normalised: sum p/h * h = 1

def gs2(va, vb, anti=False):
    H = TT + sp.diags(np.add.outer(va, vb).ravel() + W2)
    if anti:
        H = (P.T @ H @ P).tocsc()
    e, c = sla.eigsh(H, k=1, which='SA', tol=1e-12)
    c = c[:, 0]
    if anti:
        c = P @ c
    psi2 = (c**2).reshape(Ng, Ng)
    return e[0], psi2.sum(1) / h, psi2.sum(0) / h   # densities of particle 1 (spin a) and 2 (spin b)

def corners(va, vb):
    e1a, da = gs1(va); e1b, db = gs1(vb)
    e2, d2a, d2b = gs2(va, vb)
    eT, dTa1, dTa2 = gs2(va, va, True); eTb, dTb1, dTb2 = gs2(vb, vb, True)
    z = np.zeros(Ng)
    E = np.array([0.0, e1a, e1b, e2, eT, eTb])
    Na = np.array([0, 1, 0, 1, 2, 0]); Nb = np.array([0, 0, 1, 1, 0, 2])
    Da = [z, da, z, d2a, dTa1 + dTa2, z]; Db = [z, z, db, d2b, z, dTb1 + dTb2]
    return E, Na, Nb, np.array(Da), np.array(Db)

def run(fa, fb, n1, Ts=(1e-2, 2e-3, 4e-4)):
    na, nb = fa * n1, fb * n1
    v = np.r_[-1.0 / np.sqrt(x**2 + 1), -1.0 / np.sqrt(x**2 + 1)]
    for T in Ts:
        def negG(v):
            va, vb = v[:Ng], v[Ng:]
            E, Na, Nb, Da, Db = corners(va, vb)
            m = E.min(); wgt = np.exp(-(E - m) / T); Z = wgt.sum(); p = wgt / Z
            Om = m - T * np.log(Z)
            G = Om - h * (va @ na + vb @ nb)
            ga = h * (p @ Da - na); gb = h * (p @ Db - nb)
            return -G, -np.r_[ga, gb]
        r = minimize(negG, v, jac=True, method='L-BFGS-B', options={'maxiter': 1500, 'gtol': 1e-10, 'ftol': 1e-15})
        v = r.x
        print(f'  T={T} G={-r.fun:.6f} it={r.nit}', file=sys.stderr, flush=True)
    Fval = -r.fun
    # residual density error
    return Fval, r

if __name__ == '__main__':
    fa, fb = float(sys.argv[1]), float(sys.argv[2]); shape = sys.argv[3] if len(sys.argv) > 3 else 'softH'
    Zs = {'softH': 1.0, 'softHe': 2.0}[shape]
    e, n1 = gs1(-Zs / np.sqrt(x**2 + 1))
    phi = np.sqrt(n1)
    TW = h * phi @ (T1 @ phi)
    J = 0.5 * h * h * n1 @ (w(x[:, None] - x[None, :]) @ n1)
    F, r = run(fa, fb, n1)
    N = fa + fb
    Exc = F - N * TW - N * N * J
    Ex = -(fa**2 + fb**2) * J
    Ec = Exc - Ex
    lo = -2 * fa * fb * J; up = -2 * min(fa * fb, (1 - fa) * (1 - fb)) * J
    print(f"{shape} fa={fa:.3f} fb={fb:.3f} N={N:.3f} J={J:.5f} F={F:.6f} Exc/J={Exc/J:+.5f} Ec/J={Ec/J:+.5f} "
          f"lower={lo/J:+.4f} tent={up/J:+.4f} sat_up={Ec<=up+1e-6} conv={r.success} it={r.nit}", flush=True)
