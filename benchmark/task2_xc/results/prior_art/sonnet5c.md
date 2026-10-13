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

### 4. What is submitted (Laplacian-free; honest summary)
    E_xc = 1/4 E_x^HF + 3/4 E_x^SL + E_c^G
    E_x^SL : per spin channel, e_x = e_x^LDA(2n_s)/2 * [ (1-w) F_gga(s) + w F_1(s) ],  w = z^2, z = tau_W/tau (total density)
    F_gga  = 1 + kappa - kappa/(1 + mu s^2/kappa), kappa = 0.804 (Lieb-Oxford), mu = beta pi^2/3 = 0.2195 (PBE GE2/xc gradient cancellation)
    F_1    = 1 + mu s^2/(1 + b s^4)   (z = 1 branch)
             * small-s coefficient mu (same as GGA branch, continuity of the GE2 coefficient through w -> 0)
             * tail s^-2: the only power compatible with the amplitude ray (A) (section 3). Numerically immaterial (rows 3, 4 of the table).
             * b = 0.037713 solved at import so that E_x[H 1s] = -5/16 Ha exactly (derived, nothing fitted).
    w = z^2: smallest integer exponent for which the weight does not alter the second-order gradient coefficient of the GGA branch
             (z ~ O(s^2) in a slowly varying gas, so w ~ s^4).
    E_c^G  : PBE correlation with the Stoll same-spin part removed in one-orbital channels and re-assigned in two-electron-singlet
             pair regions: E_c = E_c^PBE - sum_s z_s^2 [1 - z_s'^2 (1 - zeta^2)] E_c^PBE[n_s, 0]. THIS IS PRIOR ART (opus4, round 4 =
             B95/TPSS Stoll gating + pair reassignment); I did not derive it and claim nothing for it.

What is new here vs KNOWN_CONDITIONS.md and prior_art:
 1. Condition (A) as a scaling relation on E_x of one-orbital densities, the obstruction theorem (3), its parametrisation-free
    minimax numbers (3b) and the Laplacian-extended escape (3, iii). Quantitative: any (n, grad n, tau) exchange violates (A) by
    >= 5.3% / 11.6% / 20.9% at lambda down to 0.75 / 0.5 / 0.25.
 2. Joint exactness of a single two-constant one-orbital branch on six HF/hydrogenic 1s densities (H, He, Ne, Ar, Be, Li 1s):
    h = a s^2/(1+b s^4) with (a, b) = (1.984, 16.66) reproduces -J of all six to < 0.05% (`fit7.py`; the "joint H+He-exact
    iso-orbital exchange" example from the task). I implemented it (candidate HHe) and it FAILS in molecules (below): it makes
    exchange of Ne 404 mEh too negative (`ex_eval.py`), because z-weighting also activates it in shell/valence regions that are
    not single-orbital. So exactness of an isolated one-orbital density does not transfer. That is a documented negative result.
 3. The exact-exchange partition diagnostic (`diag_x.py`, `ex_eval.py`, `decomp.py`): splits the AE error of any hybrid into an
    exchange part (vs exact HF exchange on HF densities) and a correlation part (vs CCSD(T) minus HF). Finding: for the
    best baseline-like forms both parts are 20-50 kcal/mol per atomization energy and cancel to ~3 (table below). No
    constraint-only construction escapes this cancellation; it is why one-orbital exactness (which only changes the exchange part)
    moves the total by +-5..10 kcal/mol and cannot help by itself.

## Phase 2 results
### Candidate list (every candidate I compared; dev WTMAD at def2-svp; no continuous parameter tuned on these)
| id | exchange one-orbital branch | z used in exchange | weight w | corr | AE MAE | IP MAE | WTMAD svp |
|---|---|---|---|---|---|---|---|
| c1 | none (PBE-form only) | - | - | pair | 4.37 | 7.47 | 2.00 |
| c2 | H-exact | per channel (exact spin scaling) | z^2 | pair | 6.76 | 9.65 | 2.99 |
| c3 | H-exact | total | z^2 | pair | 4.29 | 8.36 | 2.02 |
| c4 | H,He,Ne,Ar,Be,Li-1s-exact (a=1.98,b=16.7) | per channel | z^2 | pair | 11.86 | 10.58 | 4.90 |
| c5 | H-exact | total | z^3 | pair | 4.16 | 8.20 | 1.96 |
| c6 | H-exact | total | z^1 | pair | 4.45 | 8.60 | 2.09 |
Submitted: c3 (w = z^2 is the GE2-derived exponent; c5 is 0.06 better but z^3 has no derivation, so not used).
Reading: the one-orbital branch is dev-neutral (c1 vs c3 within noise); per-channel z (the form that honours exact spin scaling)
is clearly worse (c2) because open-shell atoms (spin channels with 1-2 orbitals) get extra exchange; the H,He-joint branch (c4)
over-binds molecules by +10..+24 kcal/mol (same failure sonnet4 reported with a different shape).
Non-SCF exchange-energy check on non-dev atoms (`ex_eval2.py`, HF densities, mEh error of semilocal E_x): MAE over 17 atoms
PBE-form 74.4, H-exact branch (total z, w=z^2) 30.9, (per-channel) 31.3: the branch does improve absolute atomic exchange by a factor
2.4, but AE cancellation with correlation takes that away.
