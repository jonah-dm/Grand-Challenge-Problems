import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import importlib
from pyscf.dft import libxc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rfun import R
# radial grid
x, w = np.polynomial.legendre.leggauss(400)
t = 0.5*(x+1); r = 30*t**3/(1-t+1e-3*0)**0 ; 
r = 40*(t/(1-t*0.999))**2; dr = np.gradient(r)  # crude; replaced below
# better: tanh-sinh-free mapping r = -ln(1-t)*c
c = 1.0
r = -c*np.log(1 - t*(1-1e-14)); jac = c*(1-1e-14)/(1 - t*(1-1e-14)) * 0.5
W = w*jac*4*np.pi*r**2
def shell(lam, kind, Z=1.0):
    if kind == 'H':
        n = lam*Z**3*np.exp(-2*Z*r)/np.pi; g = 2*Z*n
    else:  # Gaussian exponent a=Z: n=(2a/pi)^1.5 e^{-2a r^2}
        a = Z; n = lam*(2*a/np.pi)**1.5*np.exp(-2*a*r*r); g = 4*a*r*n
    return n, g
def pol(n, g, meta=True):
    z = np.zeros_like(n); tw = g*g/(8*n)
    a = np.array([n, g, z, z, z, tw]) if meta else np.array([n, g, z, z])
    b = np.zeros_like(a); b[0] = 1e-30*0
    return np.array([a, b])
def exact_sr(lam, kind, om, Z=1.0):
    if kind == 'H': return -lam**2*5/16*Z*R(Z**2/om**2), -lam**2*5/16*Z
    a = Z; full = -np.sqrt(a/np.pi); return lam**2*full*(1-om/np.sqrt(a+om**2)), lam**2*full
def E_sub(mod, lam, kind, Z=1.0):
    n, g = shell(lam, kind, Z)
    exc = mod.eval_xc('', pol(n, g), 1, 0, 1, None, 0)[0]
    om, al, be = mod.RSH; sr, full = exact_sr(lam, kind, om, Z)
    hf = al*(full - sr) + (al+be)*sr
    return np.sum(W*n*exc) + hf
def E_lib(code, lam, kind, Z=1.0, xtype='GGA'):
    n, g = shell(lam, kind, Z)
    om, al, be = libxc.rsh_coeff(code); hyb = libxc.hybrid_coeff(code)
    exc = libxc.eval_xc(code, pol(n, g, xtype=='MGGA'), 1)[0]
    if om:
        sr, full = exact_sr(lam, kind, om, Z); hf = al*(full-sr) + (al+be)*sr
    else:
        hf = hyb*(-lam**2*5/16*Z if kind=='H' else -lam**2*np.sqrt(Z/np.pi))
    return np.sum(W*n*exc) + hf
if __name__ == '__main__':
    print('norm check', np.sum(W*shell(1,'H')[0]), np.sum(W*shell(1,'G')[0]))
    rows = []
    for name in ['K0', 'K1', 'K2', 'K0b', 'K1b']:
        os.environ['O6_K'] = name[:2]
        if name.endswith('b'): os.environ['O6_OMEGA'] = '0.5767436061911404'; os.environ['O6_CAP'] = '1.174'
        else: os.environ['O6_OMEGA'] = '0.3687078276482022'; os.environ['O6_CAP'] = '1.804'
        spec = importlib.util.spec_from_file_location('sub'+name, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'submission.py'))
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        rows.append((name, (lambda mm: (lambda lam, kind='H', Z=1.0: E_sub(mm, lam, kind, Z)))(m)))
    for code, xt in [('PBE', 'GGA'), ('PBE0', 'GGA'), ('B3LYP', 'GGA'), ('HYB_GGA_XC_LC_WPBE', 'GGA'), ('CAM-B3LYP', 'GGA'), ('WB97X', 'GGA'), ('R2SCAN', 'MGGA')]:
        rows.append((code, (lambda c, x: (lambda lam, kind='H', Z=1.0: E_lib(c, lam, kind, Z, x)))(code, xt)))
    print('%-10s %9s %7s %7s %10s | %9s %7s | Gauss(a=1) E/exact  He+(Z=2) E' % ('func', 'Exc[H]', 'f(1/2)', 'f(1/4)', 'H2+inf kc', '', ''))
    for name, f in rows:
        e1 = f(1.0); e2 = f(0.5); e4 = f(0.25)
        err = (2*e2 + 5/32)*627.5095
        g = f(1.0, 'G', 1.0)/(-np.sqrt(1/np.pi))
        he = f(1.0, 'H', 2.0)
        print('%-10s %9.5f %7.4f %7.4f %10.2f | Gauss %6.3f  He+ %8.5f (exact -0.62500)  fG(1/2)/fG(1)*4 %6.3f' % (name, e1, e2/-0.3125, e4/-0.3125, err, g, he, 4*f(0.5,'G',1.0)/f(1.0,'G',1.0)))
