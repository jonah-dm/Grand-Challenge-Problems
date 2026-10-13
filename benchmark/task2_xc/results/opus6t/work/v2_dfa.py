"""V2: C2 on the occupation square for exact 3D one-orbital shapes (H 1s, Gaussian), existing DFAs. No SCF.
Reports E_c/J (correlation-only codes) and E_xc/J (full DFA incl. a*exact exchange = -a(fa^2+fb^2)J)."""
import numpy as np, sys
from pyscf.dft import libxc
r = np.exp(np.linspace(np.log(1e-6), np.log(60), 6000)); wr = 4*np.pi*r**2*np.gradient(r)
def shape(name):
    if name == 'H':
        n = np.exp(-2*r)/np.pi; dn = -2*n
    else:  # Gaussian a=0.5
        a = 0.5; n = (2*a/np.pi)**1.5*np.exp(-2*a*r*r); dn = -4*a*r*n
    q = np.cumsum(n*wr); vh = q/r + (np.sum(n*wr/r) - np.cumsum(n*wr/r))
    J = 0.5*np.sum(n*vh*wr)
    return n, dn, J
def ev(code, fa, fb, n, dn):
    z = np.zeros_like(n)
    def ch(f):
        nn = np.maximum(f*n, 1e-30); g = f*dn; tw = g*g/(8*nn)
        return np.array([nn, g, z, z, z, tw])
    rho = np.array([ch(fa), ch(fb)])
    t = libxc.xc_type(code)
    if t in ('LDA',): rho = rho[:, 0]
    elif t in ('GGA', 'HF'): rho = rho[:, :4]
    exc = libxc.eval_xc(code, rho, spin=1, deriv=0)[0]
    return np.sum(exc*(fa+fb)*n*wr)
CORR = {'LDA(PW92)': ',PW', 'PBE': ',PBE', 'LYP': ',LYP', 'TPSS': ',TPSS', 'SCAN': ',SCAN', 'r2SCAN': ',R2SCAN', 'B95': ',MGGA_C_BC95', 'M06-L c': ',MGGA_C_M06_L'}
XC = {'LDA': 'LDA,PW', 'PBE': 'PBE', 'B3LYP': 'B3LYP', 'PBE0': 'PBE0', 'TPSS': 'TPSS', 'SCAN': 'SCAN', 'r2SCAN': 'R2SCAN', 'M06-L': 'M06L', 'B97M-V(no NL)': 'B97M_V', 'wB97M-V': None}
PTS = [(0.5, 0), (0.25, 0.25), (0.5, 0.5), (0.75, 0.25), (0.75, 0.75), (1.0, 0.5), (1.0, 1.0)]
for shp in ['H', 'G']:
    n, dn, J = shape(shp)
    print(f"\n### shape {shp}: J = {J:.6f}")
    hdr = "| | " + " | ".join(f"({a},{b})" for a, b in PTS) + " |"
    print("E_c/J (exact particle value -2 f_a f_b; hole: <= tent -2(1-fa)(1-fb), >= -2 fa fb)")
    print(hdr); print("|---" * (len(PTS)+1) + "|")
    ex = [(-2*a*b if a+b <= 1 else -2*(1-a)*(1-b)) for a, b in PTS]
    print("| exact / tent bound | " + " | ".join(f"{e:+.3f}" for e in ex) + " |")
    for k, c in CORR.items():
        try:
            vals = [ev(c, a, b, n, dn)/J for a, b in PTS]
            print(f"| {k} | " + " | ".join(f"{v:+.3f}" + ("!" if v > e + 1e-3 else "") for v, e in zip(vals, ex)) + " |")
        except Exception as e:
            print(f"| {k} | failed {e} |")
    print("\nE_xc/J (exact particle -N^2; hole upper bound -(1+(N-1)^2))")
    print(hdr); print("|---" * (len(PTS)+1) + "|")
    exx = [(-(a+b)**2 if a+b <= 1 else -(1+(a+b-1)**2)) for a, b in PTS]
    print("| exact / bound | " + " | ".join(f"{e:+.3f}" for e in exx) + " |")
    for k, c in XC.items():
        if c is None: continue
        try:
            hyb = libxc.hybrid_coeff(c)
            vals = [(ev(c, a, b, n, dn) - hyb*(a*a+b*b)*J)/J for a, b in PTS]
            print(f"| {k} | " + " | ".join(f"{v:+.3f}" + ("!" if v > e + 1e-3 else "") for v, e in zip(vals, exx)) + " |")
        except Exception as e:
            print(f"| {k} | failed {e} |")
# T4 check: forced-form outer-region bound for hydrogenic pair, numerically
K = 5*np.pi/(4*(3*np.pi**2)**(2/3)); sstar = (3*np.pi)**(-1/3)
for Z in [1.0, 27/16]:
    n = 2*Z**3*np.exp(-2*Z*r)/np.pi; dn = -2*Z*n; kf = (3*np.pi**2*n)**(1/3); s = np.abs(dn)/(2*kf*n)
    m = s >= sstar
    print(f"T4 hydrogenic pair Z={Z:.4f}: bound = {-K*np.sum((n**(4/3)/s**2*wr)[m]):.5f}  analytic -0.5230Z = {-0.52299*Z:.5f}")
