"""Closed-form Rayleigh quotient of the two-parameter Hylleraas trial function
       psi(s,t,u) = exp(-zeta s) (1 + c u),   s = r1+r2, t = r1-r2, u = r12,
for a two-electron atom with nuclear charge Z (fixed nucleus, atomic units).
Any zeta>0, c>=0 gives a rigorous upper bound  E(Z,2) <= E_hyl(zeta,c;Z).
The result is a rational function; it is reproduced verbatim in certify.py.
Run:  python hylleraas.py   (derives the formula with sympy and checks it)."""
import sympy as sp
s,t,u,zeta,c,Z = sp.symbols('s t u zeta c Z', positive=True)
psi = sp.exp(-zeta*s)*(1+c*u)
ps, pt, pu = sp.diff(psi,s), sp.diff(psi,t), sp.diff(psi,u)
# Hylleraas: |grad_1 psi|^2 + |grad_2 psi|^2 in (s,t,u) coordinates
grad2 = 2*(ps**2+pt**2+pu**2) + 4*pu*(ps*s*(u**2-t**2) + pt*t*(s**2-u**2))/(u*(s**2-t**2))
vol = u*(s**2-t**2)                     # d tau_1 d tau_2 = 2 pi^2 u (s^2-t^2) ds dt du  (constant cancels)
def integ(f):
    f = sp.simplify(f)
    return sp.integrate(sp.integrate(sp.integrate(f,(t,-u,u)),(u,0,s)),(s,0,sp.oo))
N  = integ(psi**2*vol)
T  = integ(sp.Rational(1,2)*grad2*vol)
Vn = integ(-Z*4*s/(s**2-t**2)*psi**2*vol)    # 1/r1 + 1/r2 = 4s/(s^2-t^2)
Ve = integ(psi**2/u*vol)
E = sp.simplify((T+Vn+Ve)/N)
print("E_hyl(zeta,c;Z) =", E)
assert sp.simplify(E.subs(c,0) - zeta*(zeta - 2*Z + sp.Rational(5,8))) == 0   # c=0: screened 1s^2
Ef = sp.lambdify((zeta,c,Z),E)
print("He, zeta=1.85, c=0.365:", Ef(1.85,0.365,2), "(literature two-parameter optimum -2.8911)")
