# opus5t (THEORY track): amplitude scaling of exchange for fractional one-orbital densities

Started 00:39 UTC. Phase 1 (derivation only, no code, no SCF) 00:39-01:40. Phase 2 below.

## 0. What is used from KNOWN_CONDITIONS.md (no credit claimed)
- Uniform-gas limit, spin scaling of exchange, Lieb-Oxford (kappa=0.804), one-electron SIE
  freedom (E_x[n1]=-J[n1], E_c[n1]=0), iso-orbital indicator z=tau_W/tau, H-atom norm,
  global hybrid fraction 1/4 (adiabatic connection).
- Base functional: opus4 (round 4) = sonnet3 z-blended exchange + pair-aware Stoll correlation.
  Both are known/prior art (B95/TPSS-type ideas). They are inherited unchanged.
- Piecewise linearity (PPLB) is known. Section 1 derives the condition from the Lieb ensemble
  functional directly, but its energy-level content for N<=1 is close to PPLB. I claim novelty
  for the pointwise/amplitude form (C1), the hole statement, the no-go theorem (T2) with its unique
  enhancement factor, and the bounded construction (Section 4). I do not claim it for "fractional
  H is linear".

## 1. The condition (C1): quadratic amplitude scaling in a fractional one-orbital channel
Let phi be a normalized orbital and n1 = |phi|^2. Take a fully spin-polarized density
n = lambda n1 with 0 < lambda <= 1 (the other spin channel empty).

**Claim C1.**  E_xc[lambda n1, 0] = lambda^2 E_x[n1, 0] = -lambda^2 J[n1],  with E_c[lambda n1, 0] = 0.
Pointwise, in the conventional (hole) gauge:
  e_xc(r) = -(1/2) lambda^2 n1(r) v_H[n1](r),   h_xc(r,r') = -lambda n1(r'),
so the hole integrates to **-lambda, not -1**. The on-top value is h(r,r) = -n(r).

**Proof.** Use the Lieb (ensemble) functional F[n] = inf over ensembles Gamma -> n of Tr Gamma(T+V_ee).
If N = integral n <= 1, every admissible ensemble mixes only 0- and 1-electron states. Then V_ee = 0. For a
1-electron state with density n_k, T >= T_W[n_k], with equality for the nodeless orbital sqrt(n_k).
T_W is convex in n. So any ensemble sum_k p_k n_k = n has sum_k p_k T_W[n_k] >= T_W[n]. The bound is
attained by Gamma = (1-lambda)|0><0| + lambda|phi><phi| with phi = sqrt(n1). Hence
F[lambda n1] = T_W[lambda n1] = lambda T_W[n1]. The non-interacting functional is the same
minimisation without V_ee, so T_s[lambda n1] = lambda T_W[n1]. Therefore
E_xc = F - T_s - J = -J[lambda n1] = -lambda^2 J[n1]. The KS ensemble pair density is zero, so the
exchange hole is -n(r') and the correlation hole is zero. QED.

C1 holds for **every** one-orbital shape, not only hydrogenic densities. Its practical
meaning is this. When a single-orbital spin channel is shared between fragments (the stretched
odd-electron bond H2+, the unpaired spin in H...H...H transition states, charge transfer), each
fragment carries a fractional amplitude lambda of a one-orbital density. The fragment's
exchange must then scale as lambda^2.

## 2. Theorem T2 (no-go for semilocal exchange; new)
Write semilocal exchange (spin-scaled) as E_x = int e_x^unif(n) F(s, t), where t is any
iso-orbital indicator that is invariant under amplitude scaling (z = tau_W/tau, alpha, beta, ...).
In a one-orbital region, tau = tau_W, so z = 1 (alpha = 0) for every lambda. F is then a function
of s alone, F1(s) = F(s, z=1). Under n -> lambda n we have e_x^unif -> lambda^{4/3} e_x^unif and
s -> lambda^{-1/3} s.

**T2.** C1 holds for all lambda in (0,1] and all one-orbital n1 if and only if
F1(lambda^{-1/3} s) = lambda^{2/3} F1(s) on the realised s range, i.e. **F1(s) = c s^{-2}**.
(Proof: C1 must hold for every n1. Take the first variation of both sides with respect to a
localized change of n1 at a point with given s. Equality of the integrands follows, and the
functional equation f(lambda^{-1/3}s) = lambda^{2/3} f(s) for all lambda has only power-law
solutions with exponent -2.)
The H-atom norm fixes c: c = (5/16) / int |e_x^unif(n_H)| s^{-2}. It is then exact for every
fractional hydrogenic density lambda Z^3 e^{-2Zr}/pi, by amplitude and coordinate scaling.

Consequences:
(a) F1 = c s^{-2} is unbounded as s -> 0. So it violates the local Lieb-Oxford bound at every
    density critical point (bond midpoints, cuspless orbital centres). Analytically, for a Gaussian
    orbital density lambda (2a/pi)^{3/2} e^{-2ar^2} it gives E_x = -5 lambda^2 sqrt(a/pi), while the exact
    value is -lambda^2 sqrt(a/pi). That is 5 times too negative. **No LO-respecting (n, grad n, tau)
    meta-GGA can satisfy C1.** Every bounded F1 has amplitude degree p < 2: LDA 4/3, PBE -> 4/3,
    B88 -> 1 at large s, SCAN/r2SCAN 4/3 at small s, and sonnet3's s^{-1/2} tail 3/2.
(b) Global hybrid with a < 1 and an exact H norm: f(lambda) = E_xc[lambda n_H]/E_x[n_H]
    = a lambda^2 + (1-a) f_sl(lambda). This is not equal to lambda^2 at any lambda in (0,1), because
    f_sl(lambda) > lambda^2. Violation at lambda=1/2: (1-a)(f_sl(1/2) - 1/4).
(c) Range separation: SR semilocal exchange becomes amplitude-degree 2 where omega >> k_F,sigma(r).
    This happens because the on-top exchange hole -n_sigma(r) is exact in LSDA, and E_x^SR ->
    -(pi/2omega^2) sum_sigma int n_sigma^2. So the violation is confined to high-density regions with
    k_F,sigma > omega. Long-range HF (alpha=1) gives the interfragment part -lambda_A lambda_B/R exactly.
(d) Laplacian route: with q_L = n lap n/|grad n|^2, the local quantity y = 1/(1-q_L) equals zeta r
    for any hydrogenic density. This gives the closed form
    e_x = -(2 pi n^4/|grad n|^2)(e^{2y}-1-y)/y, which is exact pointwise for all fractional hydrogenic
    densities. It still diverges at critical points (|grad n| -> 0). I do not use it (numerically unsafe,
    and PySCF's MGGA path does not supply the Laplacian reliably).

## 3. Who violates C1 (to be verified numerically in phase 2 on the exact fractional H density)
Every semilocal functional and every global or range-separated hybrid with an SR semilocal part.
f(1/2) - 1/4 is the measure. The exact value is 0. To compute: LDA, PBE, B88, r2SCAN, SCAN, PBE0,
B3LYP, sonnet3/opus4, and this work.

## 4. Construction (bounded, partial enforcement of C1)
T2 says exact enforcement is impossible with bounded F. The best bounded choice keeps T2's
unique form where it is harmless and regularises it only where it must:
  F_q(s) = F_LO / (1 + F_LO s^2 / c_q),    F_LO = 1 + kappa = 1.804 (Lieb-Oxford, as in PBE/sonnet3)
- s -> infinity: F_q -> c_q s^{-2}, so C1 is satisfied asymptotically. This is where fractional
  fragment densities live: low density, large s.
- s -> 0: F_q -> F_LO, so the local LO bound is respected and there is no divergence at critical points.
- c_q is fixed by the H-atom norm, E_x^sl[n_H] = -5/16 exactly. The 25% HF and 75% semilocal parts
  are then each exact for H. This is not a dev-set parameter.
Where it applies: only where C1's premise holds, a single-orbital spin channel in a polarized
region. The weight is w_s = zeta^2 z_s^2 (z_s = tau_W,s/tau_s, spin-resolved).
- H atom and any fully polarized one-orbital density: zeta = 1, z_s = 1, so w = 1.
- Closed shells (RKS): zeta = 0, so w = 0 and the functional equals opus4 exactly. A singlet pair is
  not constrained by C1 (that is the flat-plane regime). There the extra binding of semilocal
  exchange is the known compensation for nondynamic correlation, which must be kept (H2 binding).
- Channel exchange: F_s = (1 - w_s) F_opus4(s_s, z) + w_s F_q(s_s).
Correlation: opus4's pair-aware Stoll form. It already satisfies the correlation half of C1,
E_c[lambda n1, 0] = 0 for every lambda (z_s = 1, zeta = 1 remove the whole polarized PBE piece).
Hybrid: HYB = 0.25 (known). No RSH and no VV10 by default.

Predicted effects (made before any run): open-shell species are raised where the unpaired spin is
iso-orbital and of low density (atomic valence tails, radicals, the TS).
- IP(B) should decrease (opus4 +6.2).
- IP(HF) should increase (opus4 -5.0).
- AE of multiply bonded molecules made from polarized atoms should increase (opus4 underbinds N2,
  CO, F2).
- C-H overbinding may get worse.
- The H+H2 barrier should rise (opus4 7.1 vs exact 9.6), because the TS has a delocalized iso-orbital spin.
Closed-shell energies and the H atom are unchanged.

Exact conditions, satisfied / violated (to confirm in phase 2):
satisfied: UEG (w=0 when z_s=0); spin scaling of the exchange weight's channel form; LO pointwise
(F <= 1.804); e_x <= 0; H atom exact (x and c); C1 for correlation at all lambda;
C1 for exchange asymptotically (large s) in polarized one-orbital channels.
violated: C1 exactly (impossible by T2); the 1.174 one-orbital bound pointwise near nuclei in
w-regions; the exchange GE2 coefficient; no long-range -1/r potential.

### Phase-1 hand results (analytic, no code)
- Pure T2 constant in spin-scaled convention: for H, s = (6 pi^2 n)^{-1/3}. So E_x = -(1/2) c_x 2^{4/3}
  (6pi^2)^{2/3} c int n_H^2 = -0.5624 c, which gives **c = 0.5557** (F1 = 0.5557 s^-2).
- Gaussian check of T2(a): the ratio of the c s^-2 energy to the exact energy is 5, for any exponent and any lambda.
- Consistency with known delocalization error: stretched H2+ is two lambda=1/2 H fragments.
  The error is 2(-5/16)(f(1/2) - 1/4). With an LDA-like amplitude degree 4/3, f(1/2) = 2^{-4/3} = 0.397,
  so the error is -0.092 Ha (-58 kcal/mol), which matches the known LSDA H2+ dissociation error
  (~50-60 kcal/mol). T2(b) says that a 25% hybrid keeps 75% of this.
- Delocalization error as amplitude-degree deficit: any functional that is exact at lambda=1 and has
  amplitude degree p gives a fractional-H error (lambda^p - lambda^2) E_x[n_H]. Only p=2 is exact.

---------------------------------------------------------------------------------------------
# Phase 2 (01:40 onward): implementation, verification, scores

## Correction to a phase-1 hand number
The "-58 kcal/mol for LDA" stretched-H2+ estimate assumed an LDA that is exact for H at lambda=1.
The computed values at the exact density, 2E_xc[n_H/2] + 5/32, are in the table below.
They are LDA -35, PBE -61, r2SCAN -57, B3LYP -54.

## Numerical verification of T2 and the C1 violations (work/norm.py, work/c1_table.py; exact analytic H density)
- T2 constant: c = 0.55556 = **5/9 exactly** (spin-scaled s). It reproduces E_x[lambda n_H] = -lambda^2 5/16 for all lambda.
- Gaussian one-orbital density: the c s^-2 functional gives 4.99999 times the exact -J (analytic prediction 5). So the
  unique C1-exact (s,z) form fails by a factor of 5 on cuspless densities, as T2(a) states.
- C1 violation on the fractional H atom (f = E_xc[lambda n_H]/(-5/16); exact f(1/2)=0.25, f(1/4)=0.0625;
  last column = delocalization error of stretched H2+ at the exact density, kcal/mol):

| functional | E_xc[H] | f(1/2) | f(1/4) | H2+(R->inf) error |
|---|---|---|---|---|
| LDA | -0.2680 | 0.340 | 0.135 | -35.4 |
| PBE | -0.3119 | 0.413 | 0.172 | -63.7 |
| B88 x | -0.3098 | 0.414 | 0.175 | -64.4 |
| TPSS x | -0.3125 | 0.406 | 0.166 | -61.0 |
| SCAN / r2SCAN | -0.3125 | 0.395 | 0.156 | -56.9 |
| PBE0 | -0.3136 | 0.373 | 0.145 | -48.4 |
| B3LYP | -0.3146 | 0.387 | 0.155 | -53.6 |
| sonnet3 / opus4 | -0.3106 | 0.372 | 0.145 | -48.0 |
| **opus5t B1 (exact-density norm)** | -0.3126 | **0.318** | **0.100** | **-26.8** |
| exact | -0.3125 | 0.250 | 0.0625 | 0 |

Every existing functional tested violates C1. The construction removes about 45% of the violation of
its own base. The rest is the price of boundedness (T2).
- Bounded variant with the small-s limit kept at sonnet3's (F = (1+mu s^2)/(1+mu s^4/c), H-normed c=20.6):
  f(1/2) = 0.419, so it gives **no** C1 improvement. Corollary of T2 checked numerically: with F(0)=1 and the
  H norm, a bounded iso-orbital F cannot reduce the fractional violation. C1 needs F(0) > 1 in polarized
  one-orbital regions, hence F_LO.

## Self-consistent H2+ (computed with the B1 constant c_q=2.2089; B2 differs by 2.7% in c_q) (aug-cc-pVTZ; exact = UHF in same basis), E(H2+) - E(H), kcal/mol (work/barrier.py)
| R (bohr) | exact | opus4 base (vtau fixed) | opus5t |
|---|---|---|---|
| 2 | -64.3 | -66.9 | -69.7 |
| 4 | -28.8 | -44.8 | -31.6 |
| 6 | -7.5 | -36.0 | -13.5 |
The stretched-bond delocalization error drops from -28.5 to -6.1 kcal/mol at R=6 bohr. That is the
physical content of C1, and the mechanism behind barrier heights and charge-transfer errors.

## HARNESS FINDING (affects every prior MGGA submission)
PySCF's define_xc_ wrapper (libxc.py `_eval_xc1`) drops `None` entries of vxc and then truncates to the
expected length. Returning vxc = (vrho, vsigma, np.zeros(N), vtau), as sonnet3/opus4 and earlier rounds did,
puts the zero "vlapl" array into the **vtau slot** and silently discards vtau. Those SCFs were therefore
not variational in tau. I found this because the C1 term (strongly tau-dependent through z_s) drove the H atom
SCF to a stationary point 5 kcal/mol above a directly minimised energy. Diagnosis: work/hdiag10.py and
hdiag11.py. The PySCF XC matrix equalled the vrho+vsigma part exactly (diff 4e-15) with no tau term.
Fix: return vlapl = None. With the fix, the Fock-matrix and finite-difference energy gradients agree to 1e-11.
For the opus4 form the effect on dev energies is small (SVP WTMAD 1.82 -> 1.83). For tau-heavy forms it is
large.

## Numerical regulator
z_s in the C1 weight uses (tau_W,s + eps tau_unif,s)/(tau_s + eps tau_unif,s) with eps = 1e-2. This removes
0/0 at points where grad n = tau = 0 (Gaussian-basis nuclei). It changes nothing where tau >> 0.01 tau_unif.
It is numerical, not fitted, and I did not vary it.

## Checks
- check_derivs: OK (check_derivs.log). audit: DISTINCT (audit.log; closest PBE0 R^2 = 0.970).
- All 24 dev species converge (QZ).

## Candidates compared (all disclosed; no continuous parameter tuned to the dev set)
"off" is the opus4 form with the vtau bug fixed. It is the base, not a submission candidate, because it contains nothing new.
| candidate | C1 weight | c_q (H norm) | SVP WTMAD | QZ WTMAD (AE / IP MAE) | H+H2 barrier, def2-QZVP (ref 9.6) |
|---|---|---|---|---|---|
| base (opus4, vtau fixed) | none | - | 1.83 | 1.10 (2.42 / 3.93) | 7.21 |
| zeta2, before the vtau fix (non-variational) | zeta^2 z_s^2 | 2.2089 | 2.16 | not run | SCF unusable |
| B1 zeta2, exact-density H norm | zeta^2 z_s^2 | 2.2089 | 1.85 | 1.31 (2.93 / 4.49) | **TS SCF did NOT converge** (score.py fallbacks, 995 s) |
| **B2 zeta2, self-consistent QZ H norm (SUBMITTED)** | zeta^2 z_s^2 | 2.1501 | not run | **1.55** (3.59 / 4.44) | **9.06** (converged) |
| bounded form with F(0)=1 | zeta^2 z_s^2 | 20.6 | not run (no C1 gain at the exact density, f(1/2)=0.419) | | |
Selection: B1 has the better dev score, and its norm is the one Section 4 defines. But B1's H3 transition state
did not converge with the harness's own fallback sequence at def2-QZVP. The README requires convergence on
every species, and the barrier is a hidden class. So I submit B2. B2 uses the same form, with the H norm
imposed self-consistently, E_SCF(H) = -0.5000000 at def2-QZVP and grid 4 (work/scnorm.py; this is still
H-atom data, not dev data). B2 converges on all dev species and on the TS. The cost is dev 1.31 -> 1.55:
removing B1's 1.1 kcal/mol H-atom overbinding moves X-H AEs up by about 1-3 kcal/mol.
This selection used dev scores and convergence. Both are disclosed.
Barrier: TS = linear H3 with r = 0.930 A, UKS, harness settings (score.run_species, def2-QZVP, grid 4).
Barrier = E(TS) - E(H) - E(H2), using the H and H2 energies from the same QZ dev runs (work/ts_conv.py).

## Honest assessment
- Dev set (first-row AEs and IPs, closed-shell molecules vs open-shell atoms): the C1 construction is
  **worse** than its base, 1.31 (B1) / 1.55 (B2, submitted) vs 1.10. Closed-shell molecules are identical by construction. The change
  comes from the atoms and radicals. Polarized atoms with large-s valence tails (N, quartet) are raised too
  much (N2 +7.7, HCN +4.2). The H atom is overbound by 1.1 kcal/mol (H2O -6.1, HF -5.8). IP(B) gets worse
  (+8.1), contrary to my phase-1 prediction, because the B atom's high-density polarized region receives F -> F_LO.
- Delocalization physics: large gains, in the direction C1 predicts. Stretched H2+ error at R = 6 bohr goes
  from -28.5 to -6.1 kcal/mol. The H+H2 barrier goes from 7.2 to about 9 (ref 9.6). These are the hidden barrier
  class and charge-transfer-type situations.
- Hidden PA and NCI items are closed shell, so they are unchanged from the base.
- Hidden second-row AEs: expect the same polarized-atom overshoot as N2 for P, S, Cl.

## Exact conditions: satisfied / violated (submitted B2)
Satisfied:
- UEG (w = 0 when z_s = 0);
- exchange spin scaling (each channel uses its own n_s, s_s, z_s and the local zeta);
- local LO (F <= 1.804);
- e_x <= 0;
- H atom exact self-consistently at def2-QZVP (E = -0.5000000; at the exact density E_x^sl is 0.8% less negative than -5/16);
- C1 for correlation at every lambda;
- C1 for exchange asymptotically (large s) in polarized one-orbital channels; at finite s, about 45% of the
  violation is removed (fractional H);
- closed-shell energies identical to the base;
- all of opus4's correlation conditions.
Violated:
- C1 exactly (impossible for bounded F, T2);
- SCAN's 1.174 one-orbital bound, pointwise near nuclei in w-regions;
- exchange GE2 (10/81);
- no -1/r potential;
- no dispersion;
- (B1 only: its relaxed H atom was 1.1 kcal/mol below exact; B2 fixes this by construction).

## What is new (summary)
1. C1: E_xc[lambda n1, 0] = -lambda^2 J[n1] for fractional fully polarized one-orbital densities. It is derived
   from the Lieb functional, and the hole normalization is -lambda. The energy content overlaps with PPLB for N <= 1.
   The pointwise amplitude-scaling statement and its use as a constraint on enhancement factors are new.
2. T2: no-go theorem. C1 for all lambda forces F(s, z=1) = c s^-2 with c = 5/9 (H norm). That form is 5x the
   exact value on Gaussians, so no LO-bounded (n, grad n, tau) functional satisfies C1. Global hybrids with a < 1
   cannot either.
3. Bounded construction F_q = F_LO/(1 + F_LO s^2/c_q), with c_q from the H norm (2.2089 at the exact density; 2.1501 self-consistent, submitted). It is weighted by
   zeta^2 z_s^2. It halves the fractional-H violation and the stretched-H2+ error, and lifts the H+H2 barrier.
   Also shown: an F(0)=1 bounded form cannot help (f(1/2) = 0.419).
4. Harness finding: the vtau slot bug in PySCF define_xc_ for vxc tuples with a zeros vlapl.

## With more time
- Restrict F_LO's high-F small-s part to low-density regions (e.g. multiply the C1 weight by a k_F,sigma/omega-type
  low-density indicator motivated by T2(c)). This removes the N-atom and B-atom overshoot while keeping the
  fragment (low-density) effect.
- A Laplacian-level C1 form (T2(d)) with a bounded regularisation at critical points.
- RSH with alpha = 1 for the interfragment part, which C1 shows is beyond any semilocal form.
- Check second-row atoms (P, S, Cl) and more barriers before trusting the hidden-set direction.

## Final dev-set score (def2-QZVP, submitted B2): WTMAD 1.55 (AE MAE 3.59, IP MAE 4.44, MCM 4.01), all 24 species converged. score.json.
This is worse than the base (1.10), sonnet3 (1.00) and B3LYP (0.97). The case for the submission is the derivation (C1, T2), the
halved delocalization error, and the H+H2 barrier (9.06 vs 7.21 for the base; ref 9.6). The base's dev advantage partly reflects
error cancellation in atoms that C1 removes.
