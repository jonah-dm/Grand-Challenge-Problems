# sonnet5c (CONSTRUCTION track) -- NOTES

(Phase-1 section, written before any SCF run. Phase-2 results appended below.)

## Phase 1: derivation

### 1. Target condition (known, on KNOWN_CONDITIONS.md): exact one-electron exchange for ALL one-orbital densities
KNOWN_CONDITIONS lists E_x[n_1] = -J[n_1] and piecewise linearity. I ask: how much of it can ANY semilocal
exchange of (n, grad n, tau) satisfy, and what does the unsatisfiable remainder cost?  Everything below was
computed by me from exact one-electron data (eigenstates of one-electron Hamiltonians with bare nuclei,
`lib.py`, no SCF): H 1s; H2+ sigma_g at R = 1, 2, 3, 4.5, 7 bohr; HeH2+-type frames. J is exact (ERI contraction).

### 2. Ensemble-amplitude scaling (derived)
Let n_1 = phi^2 (one orbital, norm 1) and consider the ensemble "lambda of an electron in phi, 1-lambda of vacuum"
(0 < lambda <= 1), density lambda n_1. Exact ensemble DFT: E[lambda] = lambda E[1] (piecewise linearity), with
E[1] = T+V because J+E_xc = 0 for one electron. Expanding E[lambda] = lambda(T+V) + lambda^2 J + E_xc[lambda n_1]:

    E_xc[lambda n_1] = E_x[lambda n_1] = - lambda^2 J[n_1],   E_c = 0.            (A)

Exact exchange is quadratic in the amplitude, a semilocal functional is not.
(A) is the one-orbital instance of piecewise linearity, but stated as a scaling relation on E_x alone it is a
condition on the exchange ENERGY FUNCTIONAL that is not on the known list in this form (the list has PWL and the
flat-plane condition in general, not the amplitude ray of one-orbital exchange or its semilocal consequences).

### 3. Obstruction theorem (derived + verified numerically)
Take any E_x = int n^{4/3} F(s) d3r on the one-orbital branch (z = 1 implies tau = tau_W, so tau adds nothing;
s is the only local scale-free variable built from n, grad n). Uniform scaling holds for any F. Under n -> lambda n,
s -> lambda^{-1/3} s, so

    E_x[lambda n_1] = lambda^{4/3} int n_1^{4/3} F(lambda^{-1/3} s_1).

(A) for all lambda and all shapes n_1 holds iff (pointwise sufficient) F(lambda^{-1/3} s) = lambda^{2/3} F(s), i.e.
F(s) = C / s^2 (e_x = -C' n^4/|grad n|^2: finite integrals at maxima, so legal). But then the shape family fails:
with C fixed by the H atom, H2+ at R=2 is off by +126% and the delocalised H2+ (R=7) by -21% (`fit5.py`).
Conversely any F that is exact for the shape family (H, H2+ R<=2) violates (A):

| F_x (z=1 branch), error of E_x vs exact | H, lambda=1 | H2+ R=3 | H2+ R=7 | H lambda=.75 | H lambda=.5 |
|---|---|---|---|---|---|
| LDA | -14.2% | -9.9% | +9.4% | +3.9% | +36% |
| PBE | -2.1% | +3.9% | +30.7% | +20.7% | +62.6% |
| H-exact one-orbital branch, s^-1/2 tail (sonnet3-type) | 0.0% | +6.5% | +34.8% | +23.8% | +67.8% |
| H-exact, s^-2 tail (the (A)-compatible tail) | 0.0% | +6.3% | +34.7% | +23.7% | +67.3% |
(`viol.py`). A 4-parameter F least-squares-fitted to shape AND lambda in {1, .75, .5} (`fit4.py`) returns errors
of -10% (H) to +11.5% (R=7) at lambda=1 and +7..+34%, +39..+75% at lambda=.75, .5: no solution.
Consequences:
 (i) Delocalisation (fractional-charge) error of semilocal exchange on one-orbital densities is not curable by
     any improvement of F(s, z): this is a no-go, not a fitting problem. Changing the tail exponent changes
     nothing (rows 3 and 4 agree to 0.2%) because the energy sits at s < 3.
 (ii) The error is SIGN-DEFINITE: semilocal exchange over-binds low-amplitude one-orbital densities.
 (iii) The escape: add the one shape invariant w = n lap(n)/|grad n|^2 (amplitude- and scale-invariant). Then
     F = s^-2 G(w) satisfies (A) exactly for every G and G is free to reproduce shapes. In this harness lap(n) is
     not passed to eval_xc (PySCF asserts `not MGGA_DENSITY_LAPL` on the Vxc path, rho has 5 rows), so the
     harness cannot use it. NUMERICAL CHECK of the escape (`fit6.py`, exact one-electron data, no SCF):
     e_x = e_x^LDA(n) s^-2 G(w), G = exp(c0 + c1 t + c2 t^2 + c3 t^3), t = w/(1+|w|), 4 constants fitted to the
     exact -J of H, H2+ (R = 1, 2, 3, 4.5, 7), HeH2+ (R = 2, 3.5): errors +2.9, -0.3, -2.7, +2.3, -8.2, 0.0, +2.2,
     +2.7 % (vs +30..35% at R = 7 for every Laplacian-free F), and (A) holds for ALL lambda by construction.
     Held-out test (fit on H, R=1, 3, 7 only; 3-parameter-shape deg-3): HeH2+ -0.7%, -0.2%, R=2 -4.5%, R=4.5 -11%.
     So the obstruction is a missing-ingredient statement: one extra scale-and-amplitude-invariant local
     variable removes it. It is a theory result, not part of the submitted (Laplacian-free) functional.

#### 3b. Parametrisation-free version (`nonparam.py`)
Ex[F] is LINEAR in the values of F on a grid, so the best any F(s) in [0, 1.804] (Lieb-Oxford cap) can do is an LP.
Minimax relative error of E_x vs exact -lambda^2 J over the library (H; H2+ R = 1..7) at amplitudes lambda:
| amplitudes included | lambda=1 | {1,.75} | {1,.75,.5} | {1,.75,.5,.25} |
|---|---|---|---|---|
| unavoidable max error of ANY F(s) | 0.0% | 5.3% | 11.6% | 20.9% |
(bonded shapes only, H and R<=3: 0% at lambda=1; with lambda in {1,.75,.5}: 7.6%.) The shape family alone is
satisfiable (a non-parametric F can be exact for all of them), the amplitude ray is what is impossible. So the
H-exact constraint can be generalised to a whole family of one-orbital shapes at no cost, and the leftover
violation is a theorem about (n, grad n, tau) functionals, not a choice of form.

### 4. What I build (CONSTRUCTION)
Since the obstruction forbids a semilocal fix, the construction goal becomes: satisfy the maximal consistent subset
of one-orbital conditions, with every constant derived, and let the global 25% exact exchange (which obeys (A)
exactly) carry the amplitude-quadratic part. Planned (details/constants in Phase 2 below).
