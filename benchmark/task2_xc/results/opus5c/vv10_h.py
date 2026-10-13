"""Derivation numerics (no SCF): VV10 nonlocal correlation evaluated on the exact H-atom density.
Radial-symmetric double integral done by angular averaging of the kernel."""
import numpy as np
def dens(r): return np.exp(-2*r)/np.pi
def grad_over_n(r): return 2.0*np.ones_like(r)   # |grad n|/n = 2 for 1s
def vv10_parts(b, C, nr=500, nmu=64):
    # radial grid (Gauss-Laguerre-like via mapping)
    x, wx = np.polynomial.legendre.leggauss(nr)
    t = 0.5*(x+1); wt = 0.5*wx
    rmax = 30.0
    r = rmax*t**3; wr = rmax*3*t**2*wt
    n = dens(r)
    wp2 = 4*np.pi*n; wg2 = C*grad_over_n(r)**4
    w0 = np.sqrt(wg2 + wp2/3)
    kap = b*1.5*np.pi*(n/(9*np.pi))**(1/6)
    beta = (1/32)*(3/b**2)**0.75
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    E_beta = np.sum(4*np.pi*r**2*wr*n*beta)
    # double integral: int dr 4pi r^2 n(r) int dr' 2pi r'^2 n(r') int dmu Phi
    R2 = r[:,None,None]**2 + r[None,:,None]**2 - 2*r[:,None,None]*r[None,:,None]*mu[None,None,:]
    g = w0[:,None,None]*R2 + kap[:,None,None]
    gp = w0[None,:,None]*R2 + kap[None,:,None]
    Phi = -1.5/(g*gp*(g+gp))
    inner = np.sum(Phi*wmu, axis=2)  # integrate over mu
    W = 4*np.pi*r**2*wr*n
    Wp = 2*np.pi*r**2*wr*n
    E_nl = 0.5*np.sum(W[:,None]*Wp[None,:]*inner)
    C6 = 1.5*np.sum(W[:,None]*W[None,:]/(w0[:,None]*w0[None,:]*(w0[:,None]+w0[None,:])))
    return E_beta, E_nl, E_beta+E_nl, C6
if __name__ == "__main__":
    for C in [0.0089, 0.0093]:
        print("C",C,"C6(H-H)=",vv10_parts(6.0,C,nr=300,nmu=32)[3], "(exact 6.499)")
    for b in [2,4,5.9,6.0,8,10,15]:
        print(b, vv10_parts(b,0.0093,nr=300,nmu=48)[:3])
