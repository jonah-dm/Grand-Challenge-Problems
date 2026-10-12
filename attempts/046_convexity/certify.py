"""Exact certification of convexity at N=2 of atomic ground-state energies:

        E(Z,1) + E(Z,3) >= 2 E(Z,2)        (fixed nucleus of charge Z, nonrelativistic)

for all real Z >= 1.45.  Everything is evaluated in exact rational arithmetic.
Ingredients (proofs in README.md):

  (E1)  E(Z,1) = -Z^2/2.
  (HY)  E(Z,2) <= E_hyl(zeta,c;Z)  for every zeta>0, c>=0   (two-parameter Hylleraas function,
        closed form derived in hylleraas.py; affine and decreasing in Z).
  (HP)  E(Z,3) >= 3 [ E_sym(Z/2,2) + E_anti(Z/2,2) ]   (spin-resolved Hall-Post bound), where
        E_sym / E_anti are the ground energies of the two-electron atom of charge Z' = Z/2 over
        spatially symmetric / antisymmetric wave functions.
  (BZ)  E_sym(Z',2) and E_anti(Z',2) are bounded below by Bazley's method of intermediate problems
        with H0 = h1 + h2 and the "special choice" p_i = r12 u_i, u_i = 1s(1) nl(2) +/- nl(1) 1s(2),
        n <= K: the bound is  min( lowest eigenvalue of Z'^2 diag(e) + Z' B^-1 S ,  -Z'^2/2 (1 + 1/(K+1)^2) ).
        B, S, e are exact rationals computed by bazley_blocks.py (data/bazley_K*.json).
  (R1)  With K=1 everything is a polynomial in Z and the inequality is proved analytically for Z >= 2.054.

Usage:  python certify.py [K]      (default K=6; takes a few minutes)
"""
from fractions import Fraction as Fr
import json, sys, math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
import sympy as sp

K = int(sys.argv[1]) if len(sys.argv)>1 else 6
raw = json.load(open('data/bazley_K%d.json'%K))
blocks = {}
for name,blk in raw.items():
    blocks[(blk['sector'],blk['l'])] = ([[Fr(x) for x in row] for row in blk['B']], [Fr(x) for x in blk['S']], [Fr(x) for x in blk['e']])

# ---------- exact linear algebra ----------
def inv(M):
    n=len(M); A=[row[:]+[Fr(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        pv=A[c][c]; A[c]=[x/pv for x in A[c]]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[a-f*b for a,b in zip(A[r],A[c])]
    return [row[n:] for row in A]

def is_psd(M):
    """Exact test that the symmetric rational matrix M is positive semidefinite (Gaussian elimination)."""
    n=len(M); A=[row[:] for row in M]
    for c in range(n):
        if A[c][c] < 0: return False
        if A[c][c] == 0:
            if any(A[c][j]!=0 for j in range(c+1,n)): return False
            continue
        for r in range(c+1,n):
            f=A[r][c]/A[c][c]
            if f!=0: A[r]=[a-f*b for a,b in zip(A[r],A[c])]
    return True

def cell_lower_bound(key, a, b):
    """Rational mu such that for every Z' in [a,b] the lowest eigenvalue of the Bazley matrix
    Z'^2 diag(e) + Z' B^-1 S  (generalized eigenproblem with metric S) is >= mu.
    Uses  Z'^2 S E >= b^2 S E  (E <= 0)  and  Z' S B^-1 S >= a S B^-1 S  (B^-1 > 0)."""
    B,S,e = blocks[key]; n=len(e)
    Binv = inv(B)
    SBS = [[S[i]*Binv[i][j]*S[j] for j in range(n)] for i in range(n)]
    A = [[ (b*b*S[i]*e[i] if i==j else Fr(0)) + a*SBS[i][j] for j in range(n)] for i in range(n)]
    An = np.array([[float(x) for x in row] for row in A]); Sn=np.diag([float(s) for s in S])
    lam = eigh(An,Sn,eigvals_only=True).min()            # floating-point guess ...
    mu = Fr(math.floor(lam*10**7)-1, 10**7)               # ... rounded down ...
    M = [[A[i][j] - (mu*S[i] if i==j else 0) for j in range(n)] for i in range(n)]
    assert is_psd(M), "exact certification failed"       # ... and certified exactly.
    return mu

def persist(b):   # lowest excluded H0 level, 1s (K+1) l, evaluated at the largest Z' of the cell
    return -b*b/2*(1+Fr(1,(K+1)**2))
def lam_sector(sec,a,b):
    return min(min(cell_lower_bound(key,a,b) for key in blocks if key[0]==sec), persist(b))

# ---------- Hylleraas upper bound ----------
def E_hyl(z,c,Z):
    return z*(c**2*(35-144*Z) + 8*c*z*(-15*Z+8*c+4) + 16*z**3 + z**2*(-32*Z+50*c+10))/(2*(48*c**2+35*c*z+8*z**2))
def hyl_params(Z):
    f=lambda p: E_hyl(p[0],p[1],Z) if p[0]>0 and p[1]>=0 else 1e9
    res=min((minimize(f,[Z-d,c0],method='Nelder-Mead',options={'xatol':1e-9,'fatol':1e-12}) for d in (0.2,0.35) for c0 in (0.2,0.45)), key=lambda r:r.fun)
    return Fr(res.x[0]).limit_denominator(2000), Fr(res.x[1]).limit_denominator(2000)

def certify_interval(Z0, Z1, h, log):
    Za = Fr(Z0); worst=None; ncells=0
    while Za < Z1:
        Zb = min(Za+h, Z1); a,b = Za/2, Zb/2
        L3  = 3*(lam_sector('sym',a,b) + lam_sector('anti',a,b))   # <= E(Z,3) for all Z in the cell
        zt,ct = hyl_params(float(Za))
        U2  = min(E_hyl(zt,ct,Za), -Za*Za/2)                        # >= E(Z,2) for all Z >= Za
        RHS = 2*U2 + Zb*Zb/2                                        # >= 2E(Z,2) - E(Z,1) on the cell
        margin = L3 - RHS
        if worst is None or margin < worst[0]: worst=(margin,Za,Zb)
        if margin < 0:
            log("FAIL on [%s,%s]"%(Za,Zb)); return False
        ncells += 1; Za = Zb
    log("certified Z in [%s, %s]: %d cells of width %s, smallest margin %.6f Ha on [%.4f, %.4f]"
        % (Z0, Z1, ncells, h, float(worst[0]), float(worst[1]), float(worst[2])))
    return True

def analytic_part(log):
    """K=1: E_sym(Z',2) >= min(-Z'^2 + 16Z'/35, -5Z'^2/8), E_anti >= -5Z'^2/8, E(Z,2) <= -(Z-5/16)^2."""
    Z = sp.Symbol('Z', positive=True); Zp = Z/2
    RHS = 2*(-(Z-sp.Rational(5,16))**2) + Z**2/2
    P1 = sp.expand(3*((-Zp**2 + sp.Rational(16,35)*Zp) - sp.Rational(5,8)*Zp**2) - RHS)
    P2 = sp.expand(3*(-sp.Rational(5,8)*Zp**2 - sp.Rational(5,8)*Zp**2) - RHS)
    Z1 = sp.Rational(2054,1000)
    for name,P in (("P1",P1),("P2",P2)):
        roots = sp.solve(P,Z); rmax = max(roots, key=float)
        assert P.coeff(Z,2) > 0 and sp.simplify(rmax - Z1) < 0 and P.subs(Z,Z1) > 0
        log("%s(Z) = %s >= 0 for Z >= %s (largest root %.6f)" % (name, P, Z1, float(rmax)))
    log("analytic part: E(Z,1)+E(Z,3)-2E(Z,2) >= min(P1,P2)(Z) >= 0 for all Z >= %s" % Z1)
    return True

if __name__ == "__main__":
    lines=[]
    def log(s): print(s, flush=True); lines.append(s)
    log("K = %d (pair functions 1s nl, n <= K, l <= n-1, both symmetry sectors)" % K)
    ok = analytic_part(log)
    ok &= certify_interval(Fr(145,100), Fr(150,100), Fr(1,2000), log)
    ok &= certify_interval(Fr(150,100), Fr(160,100), Fr(1,1000), log)
    ok &= certify_interval(Fr(160,100), Fr(180,100), Fr(1,400), log)
    ok &= certify_interval(Fr(180,100), Fr(2054,1000), Fr(1,200), log)
    log("RESULT: " + ("E(Z,1)+E(Z,3) >= 2E(Z,2) certified for all Z >= 1.45" if ok else "FAILED"))
    open('data/certificate_K%d.log'%K,'w').write("\n".join(lines)+"\n")
