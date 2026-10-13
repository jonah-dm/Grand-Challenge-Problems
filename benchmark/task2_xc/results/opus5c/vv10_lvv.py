import numpy as np
from vv10_h import *
def parts_general(nfun, gfun, b, C, nr=300, nmu=48, rmax=30.0):
    x, wx = np.polynomial.legendre.leggauss(nr); t=0.5*(x+1); wt=0.5*wx
    r = rmax*t**3; wr = rmax*3*t**2*wt
    n = nfun(r); gn = gfun(r)
    w0 = np.sqrt(C*gn**4 + 4*np.pi*n/3)
    kap = b*1.5*np.pi*(n/(9*np.pi))**(1/6)
    beta = (1/32)*(3/b**2)**0.75
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    R2 = r[:,None,None]**2 + r[None,:,None]**2 - 2*r[:,None,None]*r[None,:,None]*mu[None,None,:]
    g = w0[:,None,None]*R2 + kap[:,None,None]; gp = w0[None,:,None]*R2 + kap[None,:,None]
    inner = np.sum(-1.5/(g*gp*(g+gp))*wmu, axis=2)
    W = 4*np.pi*r**2*wr*n
    Eb = np.sum(W*beta); Enl = 0.5*np.sum(W[:,None]*(W[None,:]/2)*inner)
    Eloc = np.sum(W*(-(3*np.pi**2/32)*n/(w0*kap)**1.5))
    return Eb, Enl, Eb+Enl, Eloc, Eb+Eloc
for Z,name in [(1.0,'H'),(1.6875,'He-like(2e, zeta=1.6875)')]:
    N = 1 if name=='H' else 2
    nf = lambda r: N*Z**3/np.pi*np.exp(-2*Z*r); gf = lambda r: 2*Z*np.ones_like(r)
    for b in [4,6,8]:
        print(name, b, ["%.5f"%v for v in parts_general(nf,gf,b,0.0093)], "(beta, nonlocal, VV10 total, local-approx nonlocal, local total)")
