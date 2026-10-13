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
     the one-orbital library errors by < 0.2% (rows 3 and 4) because the library energy sits at s < 3. (In molecules the
     branch tail DOES matter, see the QZ table below: that statement is only about the library.)
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
    E_xc = 1/4 E_x^HF + 3/4 E_x^SL + E_c^G          (HYB = 0.25, META = True; no RSH, no NLC)
    E_x^SL : per spin channel, e_x = e_x^LDA(2n_s)/2 * [ (1-w) F_gga(s) + w F_1(s) ],  w = z^2, z = tau_W/tau (TOTAL density)
    F_gga  = 1 + kappa - kappa/(1 + mu s^2/kappa), kappa = 0.804 (Lieb-Oxford), mu = beta pi^2/3 = 0.2195 (PBE-type GE2/xc gradient cancellation)
    F_1    = 1 + mu s^2/(1 + C s^3), C = 4/(5 kappa/mu)^1.25 (z = 1 branch)   -- THIS BRANCH IS THE PRIOR-ART ONE (sonnet3/opus4);
             it gives E_x[H] = -0.3098 (exact -0.3125, 0.85% off), i.e. it is NOT H-exact. Kept because it scored best (see table).
    E_c^G  : PBE correlation with Stoll same-spin part removed in one-orbital channels and re-assigned in two-electron-singlet pair
             regions: E_c = E_c^PBE - sum_s z_s^2 [1 - z_s'^2 (1 - zeta^2)] E_c^PBE[n_s, 0]   (PRIOR ART: opus4 = B95/TPSS Stoll gate + pair
             reassignment).
So the submitted FUNCTIONAL is a re-derivation of opus4's (dev WTMAD 1.08 reproduced exactly). I claim NO novelty for the functional.
The novelty I claim is the analysis around it (below). I could not find a Laplacian-free construction that beats it; the attempts and
their scores are all listed.

What I claim as new vs KNOWN_CONDITIONS.md and prior_art:
 1. Condition (A) as a scaling relation on E_x of nodeless one-orbital densities (valid because any nodeless phi is the ground state of
    some one-electron potential, so E[lambda] = lambda E[1] holds in that potential), the obstruction theorem (sec. 3), its
    parametrisation-free minimax numbers (3b) and the Laplacian-extended escape (3, iii). Quantitative: any (n, grad n, tau) exchange
    violates (A) by >= 5.3% / 11.6% / 20.9% at lambda down to 0.75 / 0.5 / 0.25 (library-restricted, discretised LP, so a numerical
    theorem, not a proof). The sign-definite over-binding (PBE +62% at lambda = 0.5) is a statement about E_x alone.
 2. Joint exactness of ONE two-constant one-orbital branch on six HF/hydrogenic 1s densities (H, He, Ne, Ar, Be, Li 1s):
    h = a s^2/(1+b s^4), (a, b) = (1.984, 16.66) reproduces -J of all six to < 0.05% (`fit7.py`) -- the "joint H+He-exact iso-orbital
    exchange" example from the task, extended to cores. Implemented as candidate c4. NEGATIVE RESULT: it makes molecules over-bind by
    +10..+24 kcal/mol and atomic exchange of Ne 404 mEh too negative (`ex_eval.py`), because the z-weight also activates the branch in
    shell and valence regions that are not single-orbital. Exactness of an isolated one-orbital density does not transfer.
 3. An exact-exchange partition diagnostic (`diag_x.py`, `diag_all.py`, `ex_eval.py`, `decomp.py`) that splits the AE error of a
    hybrid into an exchange part (vs exact HF exchange on HF densities) and a correlation part (vs CCSD(T) minus HF). For the
    PBE-form-exchange/pair-correlation functional (c1) both parts are 5..54 kcal/mol per atomization energy (N2: +48.7 / -52.3,
    H2O +19.8 / -19.4, F2 +48.0 / -53.7) and cancel to a total MAE of 3.5. No constraint-only change to exchange alone can improve the
    total unless the correlation part moves with it; that is why every one-orbital exchange variant moves AE by +-5..10 kcal/mol
    without a systematic gain.
 4. Non-SCF exchange-energy check on 17 non-dev atoms (`ex_eval2.py`, HF densities, mean |E_x^SL - E_x^HF|): PBE-form 74 mEh,
    z-blend H-exact branch (total z, w = z^2) 31 mEh, a 2.4x improvement of absolute atomic exchange that the AE cancellation (item 3)
    does not convert into better atomization energies.

## Phase 2 results
### Candidate list (every candidate compared)
Dev WTMAD. svp = def2-svp, qz = def2-qzvp. All use HYB = 0.25 and the pair-aware correlation above unless stated.
| id | one-orbital exchange branch F_1 | z in exchange | weight w | AE MAE | IP MAE | WTMAD |
|---|---|---|---|---|---|---|
| c1 | none (PBE-form only) | - | - | 4.37 | 7.47 | svp 2.00 ; qz 1.50 (AE 3.54, IP 3.93) |
| c2 | H-exact, s^-2 tail (b = 0.0377) | per channel (exact spin scaling) | z^2 | 6.76 | 9.65 | svp 2.99 |
| c3 | H-exact, s^-2 tail (b = 0.0377) | total | z^2 | 4.29 | 8.36 | svp 2.02 ; qz 1.53 (AE 3.61, IP 3.98) |
| c4 | six-1s-exact (a = 1.98, b = 16.7) | per channel | z^2 | 11.86 | 10.58 | svp 4.90 |
| c5 | H-exact, s^-2 tail | total | z^3 | 4.16 | 8.20 | svp 1.96 |
| c6 | H-exact, s^-2 tail | total | z^1 | 4.45 | 8.60 | svp 2.09 |
| c7 (SUBMITTED) | prior-art branch, C s^3 (not H-exact) | total | z^2 | - | - | **qz 1.08 (AE 2.39, IP 3.89, MCM 3.14)** |
Selection of c7 used the QZ dev score among {c1, c3, c7}; no continuous parameter was tuned. Reading: (i) per-channel z (the form that
honours exact spin scaling) is clearly worse (c2) because spin channels of open-shell atoms with 1-2 orbitals get extra exchange;
(ii) my H-exact s^-2-tail branch is no better than no branch at all (1.53 vs 1.50) while the prior-art branch with the
slower-decaying tail is 0.4 better at QZ: the large-s part of the branch carries real atomization-energy information in molecules,
which the one-orbital library (energy at s < 3) cannot see. So derived-tail (A-compatible, s^-2) loses to the inherited tail, and I
report that as a negative result, not a feature.

### Own validation set (not dev, not hidden), CCSD(T)/CBS refs from reference.py (`val_set.py`, `val_ref*.json`), def2-tzvp
| item | ref | sonnet5c err | PBE0 err |
|---|---|---|---|
| AE2:HCl | 107.20 | -1.20 | -3.06 |
| AE2:H2S | 183.20 | +0.01 | -4.12 |
| AE2:PH3 | 241.36 | +1.96 | -4.25 |
| BAR:H+H2 | 9.59 | -2.31 | -3.90 |
| PA:NH3 | 212.24 | +1.55 | +1.62 |
| PA:H2O | 172.09 | +1.51 | +1.82 |
| class MAE (AE2 / BAR / PA) | | 1.06 / 2.31 / 1.53 | 3.81 / 3.90 / 1.72 |
| MCM | | 1.63 | 3.14 |
(PBE0 is only a reference point; same basis. The validation set is tiny (3 AE2, 1 barrier, 2 PA) and was not used to choose anything.)

### Audit / derivatives
check_derivs: DERIVATIVES OK (complex-step derivatives, spin-unpolarized and polarized). audit.py: DISTINCT (closest single PBE0, R^2 = 0.970).

## Exact conditions
Satisfied: UEG limit; exact spin scaling of exchange is VIOLATED by the total-z weight (open shell only; the exact-spin-scaling variant c2 is
worse); Lieb-Oxford (F <= 1.804 in both branches, audit 0 violations); E_c = 0 for fully polarised one-orbital densities; same-spin E_c = 0 for
one-orbital channels with pair reassignment (prior art); GE2 coefficient of exchange is mu_PBE rather than 10/81 (xc-level cancellation).
Violated: ensemble-amplitude scaling (A) (provably unavoidable for any Laplacian-free form, sec. 3); exact one-electron exchange for H
(0.85% off); piecewise linearity / derivative discontinuity beyond what 25% EXX gives; no dispersion (NCI will be weak); one-electron
exchange is not exact pointwise.

## Weaknesses
Atomization energies still rest on x/c cancellation (sec. 4 item 3). NCI not addressed (no NLC; VV10 was not scored, so no claim).
IP of B +6, HF cation -5. Barrier (H+H2) and PA numbers are in the validation table above.

## What I would try with more time
(1) Implement the Laplacian-extended one-orbital branch F = s^-2 G(w) (sec. 3 iii) by supplying nabla^2 n through a custom NumInt (PySCF 2.14 does not
support the Laplacian on the Vxc path, so it needs a small patch outside this harness): it satisfies (A) exactly and fits the whole
H / H2+ / HeH2+ library to <= 3% (R = 4.5: -8%) with 4 constants, vs 30%+ for any s-only form at R = 7. (2) Re-fit the one-orbital weight
so that the exact-exchange partition (sec. 4 item 3) is balanced against a correlation functional trained to the same partition. (3) Local
hybrids with the exact one-orbital exchange density -n v_H/2, which would carry (A) exactly.
