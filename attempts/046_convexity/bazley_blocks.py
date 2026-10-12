"""Exact Bazley ('special choice') data for the two-electron atom
      H(Z') = h1 + h2 + 1/r12,   h = -Delta/2 - Z'/r,
in the scaled variable x = Z' r (so everything is computed at Z'=1 and rescaled):
  S(Z') = Z'^-6 S1,  B(Z') = Z'^-7 B1,  eps(Z') = Z'^2 e.
Basis: unnormalized hydrogenic products u = 1s(1) nl(2) +/- nl(1) 1s(2) (pair functions),
       l=0 symmetric sector (incl. 1s1s), and antisymmetric sector for l=0..lmax.
All matrix entries are rational numbers.
B_ij = <u_i| r12 |u_j>,  S_ij = <u_i|u_j> (diagonal), e_i = -(1 + 1/n_i^2)/2.
"""
import sympy as sp, json, sys, time
r1, r2, r = sp.symbols('r1 r2 r', positive=True)

def R(n, l):
    # unnormalized hydrogenic radial function for Z=1 (rational coefficients)
    rho = 2*r/n
    return sp.expand(rho**l*sp.exp(-rho/2)*sp.assoc_laguerre(n-l-1, 2*l+1, rho))

def c_l(l, a, b):   # r_< = a, r_> = b
    return (a**l/b**(l+1))*(a**2/sp.Integer(2*l+3) - b**2/sp.Integer(2*l-1))

def dbl(F, l):
    """int int F(r1,r2) c_l(r1,r2) r1^2 r2^2 dr1 dr2 over the full quadrant."""
    F = sp.expand(F)
    I1 = sp.integrate(sp.integrate(F*c_l(l,r1,r2)*r1**2*r2**2, (r1,0,r2)), (r2,0,sp.oo))  # r1<r2
    I2 = sp.integrate(sp.integrate(F*c_l(l,r2,r1)*r1**2*r2**2, (r2,0,r1)), (r1,0,sp.oo))  # r2<r1
    return sp.nsimplify(sp.simplify(I1+I2))

def rad_int(F):
    return sp.simplify(sp.integrate(sp.expand(F)*r**2, (r,0,sp.oo)))

def block(l, ns, sigma):
    """B, S, e for pair functions 1s nl (n in ns) with symmetry sigma (+1 sym, -1 antisym)."""
    R1 = R(1,0)
    Rn = {n: R(n,l) for n in ns}
    N1 = rad_int(R1**2)
    K = len(ns)
    B = sp.zeros(K); S = sp.zeros(K)
    for i,n in enumerate(ns):
        for j,m in enumerate(ns):
            if j < i: B[i,j]=B[j,i]; S[i,j]=S[j,i]; continue
            if l == 0 and n == 1 and m == 1:
                # u = 1s1s (not symmetrized twice)
                B[i,j] = dbl(R1.subs(r,r1)**2*R1.subs(r,r2)**2, 0)
                S[i,j] = N1*N1
            elif l == 0 and (n == 1 or m == 1):
                mm = m if n == 1 else n
                # <1s1s| r12 | (1s mm)_sym > with (1s mm)_sym = 1s(1)mm(2) + mm(1)1s(2): = 2 * direct
                B[i,j] = 2*dbl(R1.subs(r,r1)**2*R1.subs(r,r2)*Rn[mm].subs(r,r2), 0)
                S[i,j] = 0
            else:
                direct = dbl(R1.subs(r,r1)**2*Rn[n].subs(r,r2)*Rn[m].subs(r,r2), 0)
                exch   = dbl(R1.subs(r,r1)*Rn[m].subs(r,r1)*Rn[n].subs(r,r2)*R1.subs(r,r2), l)/(2*l+1)
                B[i,j] = 2*(direct + sigma*exch)      # u = 1s(1)n(2) + sigma n(1)1s(2)
                S[i,j] = 2*(N1*rad_int(Rn[n]*Rn[m]) + sigma*rad_int(R1*Rn[m])*rad_int(R1*Rn[n]))
    e = [ -sp.Rational(1,2)*(1+sp.Rational(1,n**2)) for n in ns ]
    return B, S, e

if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv)>1 else 4
    out = {}
    t0=time.time()
    out[('sym',0)] = block(0, list(range(1,K+1)), +1)
    print("sym l=0 done", time.time()-t0, flush=True)
    for l in range(0, K):
        ns = [n for n in range(max(2,l+1), K+1)]
        if not ns: continue
        out[('anti',l)] = block(l, ns, -1)
        print("anti l=%d ns=%s done"%(l,ns), time.time()-t0, flush=True)
        if l >= 1:
            out[('sym',l)] = block(l, ns, +1)
            print("sym  l=%d ns=%s done"%(l,ns), time.time()-t0, flush=True)
    js = {"%s_l%d"%k: {"sector":k[0],"l":k[1],"B":[[str(x) for x in row] for row in v[0].tolist()],"S":[str(v[1][i,i]) for i in range(v[1].shape[0])],"e":[str(x) for x in v[2]]} for k,v in out.items()}
    json.dump(js, open("data/bazley_K%d.json"%K,"w"), indent=1)
    # sanity: sym block vs normalized values from before (B11 = 35/16 at Z'=1)
    B,S,e = out[('sym',0)]
    print("B11/S11 =", B[0,0]/S[0,0], "(expect 35/16)")
    print("S diag:", [S[i,i] for i in range(S.shape[0])])
    print("S offdiag zero?", all(S[i,j]==0 for i in range(S.shape[0]) for j in range(S.shape[0]) if i!=j))
    for key,(B,S,e) in out.items():
        print(key, "B=",B.tolist(), "S=",[S[i,i] for i in range(S.shape[0])], "e=",e)
