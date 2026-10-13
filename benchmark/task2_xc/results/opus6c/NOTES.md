# opus6c (CONSTRUCTION track, round 6): a bounded, Laplacian-free semilocal exchange that is exactly amplitude-quadratic in one-orbital channels, made possible by range separation

Started 03:13 UTC. Phase 1 (derivation only, no code, no SCF): 03:13-04:13. Phase 2 is appended below.

## 0. What is used and not claimed (from KNOWN_CONDITIONS.md and prior_art)
- UEG limit, exchange spin scaling, uniform coordinate scaling, local Lieb-Oxford cap (1+kappa = 1.804),
  one-electron SIE freedom, iso-orbital indicator z_s = tau_W,s/tau_s, H-atom norm, PBE0's 1/4 (adiabatic connection).
- Long-range-corrected hybrids with 100% LR exact exchange (LC-wPBE, wB97X, ...) are known constructions. ITYH
  attenuation of a GGA/meta-GGA enhancement factor (Iikura-Tsuneda-Yanai-Hirao 2001) is known.
- C1 (E_xc[lambda n1, 0] = -lambda^2 J[n1]) is KNOWN: Zhang & Yang, JCP 109, 2604 (1998); Mori-Sanchez, Cohen & Yang, JCP 125, 201102
  (2006); Perdew et al., PRA 76, 040501 (2007), including the open-system hole sum rule (hole integrates to -lambda). This is the organiser's
  literature check, results/round5_literature_check.md. The impossibility of a semilocal functional exact for all one-electron densities
  is known in a stronger form: Santra & Perdew, JCP 150, 174106 (2019). The only round-5 parts not found in the literature are the
  sharp Laplacian-free characterization F(s, z=1) = (5/9) s^-2 with its closed form, and the Gaussian counterexample. I use all of
  these and claim none of them. I do not use or cite the round-5 "11.6% LP floor", which is a grid artifact.
- Known related work on the remedy side: LC hybrids reduce delocalization and fractional-charge errors (Vydrov, Scuseria & Perdew
  2007; Cohen, Mori-Sanchez & Yang 2007). A Laplacian route to one-electron SIE reduction exists (Becke-Roussel 1989; Ramasamy et al.,
  PRB 112, L161112, 2025). Neither gives a bounded semilocal SR branch that is exactly amplitude-quadratic, which is what I add.
- opus5c Corollary 4 / section 1i (prior art) already states that with erfc(omega r)/r screening the one-orbital SR branch must be
  s^-2 Psi(2 omega n_s/|grad n_s|), and it gives the hydrogenic Psi = (5/9) R. That is not mine. opus5c recorded it and did not
  implement it. It also did not notice that this branch is bounded. Sections 2-5 below, and the functional, are what I add.
- Base pieces, inherited unchanged and claimed by no one here: sonnet3's z-blended exchange enhancement (PBE form with a one-orbital
  s^-1/2 tail) and opus4's pair-aware Stoll correlation (B95/TPSS-type ideas).

## 1. Setting
A spin channel holds one orbital: n_s = lambda |phi|^2, 0 < lambda <= 1, and the other channel is empty. C1 says the xc energy
must be exactly quadratic in lambda. The exact exchange hole of such a channel is h(r,r') = -n_s(r') at every reference point r.
So for any interaction kernel K(u), the exact channel exchange energy density (conventional gauge) is
    e_x^K(r) = -(1/2) n_s(r) * integral n_s(r') K(|r-r'|) dr'      (*)
This is quadratic in amplitude pointwise. Every semilocal exchange of (n, grad n, tau) has e_x = e_x^LDA(n) F(s, z) with z = 1 here,
so it scales as lambda^{4/3} F(lambda^{-1/3} s). That is the source of T2.

## 2. Short-range expansion of one-orbital SR exchange (derivation; the second-order term is Laplacian-free)
Take K = erfc(omega u)/u. Taylor-expand n_s(r+u) and average over angles:
    n + u^2 lap n/6 + u^4 lap^2 n/120 + ...
The moments are M_k = 4 pi int_0^inf u^{2k+1} erfc(omega u) du = 4 pi Gamma(k+3/2)/((2k+2) sqrt(pi) omega^{2k+2}), giving
M_0 = pi/omega^2, M_1 = 3 pi/(4 omega^4) and M_2 = 5 pi/(4 omega^6). From (*):
    e_x^SR = -(pi/(2 omega^2)) n_s^2 - (pi/(16 omega^4)) n_s lap n_s - (pi/(192 omega^6)) n_s lap^2 n_s - ...
Integrate by parts (int n lap n = -int |grad n|^2). The energy through O(omega^-4) then becomes
    E_x^SR[one orbital] = -(pi/(2 omega^2)) int n_s^2 [ 1 - b/2 ] + O(omega^-6),   b := |grad n_s|^2 / (4 omega^2 n_s^2).   (E2)
Statement: through second order in 1/omega, the exact SR exchange of ANY one-orbital channel is a Laplacian-free GGA in the single
amplitude-invariant variable b, times n_s^2. It is exactly amplitude-quadratic.
Provenance: this follows from the known Taylor expansion of the exchange hole (Becke's hole curvature with D_s = 0 for one orbital;
large-mu expansions of SR exchange by Gill et al. and Toulouse, Gori-Giorgi and Savin). What I add is the observation that for one
orbital the O(omega^-4) term becomes Laplacian-free after integration by parts, and that it is amplitude-invariant.
The O(omega^-6) term needs int (lap n)^2. That is where shape enters.

## 3. Hydrogenic closed form (re-derived; the function R is opus5c's)
Fourier transform of n_Z = Z^3 e^{-2Zr}/pi: n(k) = 16 Z^4/(4Z^2+k^2)^2. The exact SR self-exchange is then
    E_x^SR[lambda n_Z] = -lambda^2 (5/16) Z R(b),   b = Z^2/omega^2,
    R(b) = 1 - (32/(5 pi)) I(b),   I(b) = int_0^inf e^{-b u^2} (1+u^2)^{-4} du.
For a hydrogenic density |grad n|/n = 2Z everywhere, so b = |grad n_s|^2/(4 omega^2 n_s^2) is the same local variable as in (E2).
The pointwise branch below is therefore exact for every hydrogenic ion, at every lambda and every omega:
    e_x^{SR,1}(r) = -(pi/(2 omega^2)) n_s^2 * 5 Q(b),   Q(b) := R(b)/b.                  (B)
Equivalently, in enhancement-factor language (spin-scaled k_s = (6 pi^2 n_s)^{1/3}, s_s^2 = b omega^2/k_s^2):
    F^{SR,1} = (5/9) s_s^{-2} R(b) = (5/9) (k_s/omega)^2 Q(b).
Checks:
- Small b (large omega): I = 5pi/32 - (pi/32) b + (pi/64) b^2 - (5pi/192) b^3 + ..., so R = b/5 - b^2/10 + b^3/6 and
  5Q = 1 - b/2 + 5b^2/6. The O(b) term is -1/2, identical to the universal coefficient in (E2). So the hydrogenic branch is exact for
  EVERY one-orbital shape through O(omega^-4). Only from O(omega^-6) on does it use the hydrogenic value of int (lap n)^2 (I checked
  that the hydrogenic O(b^2) coefficient 5/6 reproduces -(pi/192 omega^6) int (lap n)^2 with int n^2 (1-1/r)^2 = 5 int n^2).
- Large b (tails, small omega): R -> 1, so F -> (5/9) s^-2, which is T2's form. As omega -> 0 the branch tends to the
  full-range C1 form, as it must.

## 4. Lemmas (new, subject to the scope correction below)
Scope (from the organiser's E2, which I accept): no semilocal functional, with or without the Laplacian, can be exact for all one-electron
densities even at lambda = 1. The reason is the nonlocal interfragment Hartree term (Santra-Perdew). So what follows removes only the
AMPLITUDE part of the obstruction. In this construction the interfragment part is assigned to LR exact exchange, where it is exact at
long range. The lambda = 1 shape error of the SR branch remains for non-hydrogenic shapes and carries over unchanged (times lambda^2)
to every lambda.
**L1 (exact C1 for every shape).** Any branch of the form e = n_s^2 g(b) with fixed omega is exactly amplitude-quadratic pointwise,
because n_s -> lambda n_s leaves b invariant. Hence E_x^{SR,1}[lambda n1] = lambda^2 E_x^{SR,1}[n1] for every one-orbital shape n1,
not only hydrogenic ones. The fractional-charge curvature of such a channel is then zero. The only remaining one-electron error is the
lambda = 1 shape error, which scales as lambda^2. In this way delocalization error is reduced to integer-occupation self-interaction.
**L2 (boundedness: range separation evades T2).** At fixed density, s_s x = omega/k_s with x = b^{-1/2}. So F^{SR,1} = (5/9)(k_s/omega)^2 Q(b).
R is increasing and concave with R(0) = 0 and R'(0) = 1/5, so Q is decreasing. Therefore
    sup over all gradients of F^{SR,1}(n_s, .) = (k_s/omega)^2/9,   attained at b -> 0 (density critical points),
and at that point e_x -> -(pi/(2 omega^2)) n_s^2. This is exactly the on-top (LSDA-exact) SR exchange. The T2 divergence s^-2 at
critical points is replaced by a finite, physically exact value. At large s, F -> 0.
**L3 (corollary: C1 + finiteness force 100% long-range exact exchange in one-orbital channels).** The complementary LR branch that is exact
for hydrogenic densities is (5/9) s^-2 (1 - R(b)) = (5/9)(k_s/omega)^2 (1-R(b))/b. It diverges as b -> 0, i.e. at every density
maximum, cusp-less nucleus or bond midpoint. In general a branch F = s^-2 Psi(b) is bounded at fixed density only if Psi = O(b), and
hydrogenic exactness gives Psi_LR(0) = 5/9 != 0. So no bounded semilocal LR exchange satisfies C1 even for hydrogenic densities. The
LR part must be carried by exact exchange with coefficient 1. (LC hybrids use alpha_LR = 1 for the -1/r asymptote; that is known.
The derivation of alpha_LR = 1 from C1 plus pointwise finiteness is, to my knowledge, new.)
**L4 (omega from C1 + Lieb-Oxford on the hydrogen atom; derived, not fitted).** For the exact H density, b = 1/omega^2 everywhere and
k_s is largest at the nucleus, where s_s = s0 = (6 pi)^{-1/3}. So
    max_r F^{SR,1}[n_H] = (5/9)(6 pi)^{2/3} R(1/omega^2) = 3.937 R(1/omega^2).
The pointwise exchange cap |e_x^SR| <= |e_x| <= (1+kappa)|e_x^LDA| holds because erfc <= 1, the hole is negative and LO caps the full
e_x. It requires 3.937 R(1/omega^2) <= 1.804. This gives an admissible interval omega >= omega*. omega -> infinity is full HF, which
destroys the semilocal error cancellation that molecules need. So I take the boundary
    omega* = the smallest omega with sup_r F^{SR,1}[n_H] = 1 + kappa,
which keeps the largest semilocal range consistent with C1 and LO on the one-electron reference system. The phase-1 estimate from my
R(0.3) and R(0.4) values is omega* ~ 0.37. I will solve for it exactly in phase 2. It does not depend on the SR exact-exchange fraction.
Scope: for He+-like ions (Z > 1) the bound is exceeded near the nucleus, and for cusp-less Gaussian-basis nuclei b -> 0 gives
F -> (k/omega)^2/9. A smooth cap of F at 1+kappa therefore guards it (numerical safeguard, exponent fixed a priori, disclosed).

## 5. The functional (construction)
    E_xc = E_x^{LR-HF}(omega*) + a E_x^{SR-HF}(omega*) + (1-a) E_x^{SR,sl}(omega*) + E_c^{opus4}
    RSH = (omega*, 1, -(1-a))   [alpha*HF + beta*SR-HF]
SR semilocal exchange, per spin channel (spin scaling exact):
    e^{SR,sl}_s = (1 - w_s) e^{SR,base}_s + w_s e^{SR,1}_s
- e^{SR,base}: the sonnet3 z-blended enhancement factor (prior art), ITYH-attenuated with k_eff = k_F/sqrt(F) (known).
- e^{SR,1}: branch (B), smoothly capped at F <= 1.804.
- w_s: where C1's premise holds, a single-orbital channel with the other channel absent: w_s = zeta^2 z_s^2 (opus5t's weight, known).
  With LC, the interfragment exchange (-lambda_A lambda_B/R) is already exact through LR-HF, so the branch only has to fix the
  on-fragment amplitude degree. It is exactly quadratic (L1), so no compromise in lambda is needed, unlike F_q.
- Correlation: opus4 pair-aware Stoll PBE (prior art). It gives E_c[lambda n1, 0] = 0 at every lambda (z_s = 1, zeta = 1).
Pre-registered candidates (no continuous parameter tuned on dev). All will be disclosed with scores:
  K0  LC base only (w = 0), a = 1/4                       [reference: contains no C1 branch]
  K1  LC + branch, w = zeta^2 z_s^2, a = 1/4 (PBE0 fraction)  [primary]
  K2  LC + branch, w = zeta^2 z_s^2, a = 0 (pure LC)
  K3  LC + branch, w = z_s^2 (all one-orbital channels incl. closed shells), a = 1/4  [expected to underbind X-H]
Selection rule: dev SVP WTMAD among K1-K3, subject to convergence on every species, H2+ at R = 6 bohr and the H+H2 barrier.
All of these are disclosed.

## 6. Predictions (before any run)
- Fractional H: exact for all lambda and every omega (by construction; the cap could perturb the nuclear region slightly).
- H2+ at R = 6 bohr: interfragment exchange is exact (LR-HF), and the on-fragment part is exactly quadratic. The error should come from
  SCF density relaxation and the cap only: |err| < 3 kcal/mol (opus4 base -28.5, opus5t -6.1).
- H+H2 barrier: LC alone raises it. With the branch I expect 9-10 (ref 9.6).
- Atoms: in valence regions the branch is the hydrogenic SR exchange with the local decay constant, and it falls back to the base
  where z_s < 1. The polarized-atom overshoot of F_q (F -> 1.804 at small s irrespective of omega) should be much smaller, because
  F^{SR,1} is reduced by R(b) with b ~ 7-10 in valence.
- Risk: with 100% LR-HF, closed-shell AEs may underbind (more exact exchange). K2 (a = 0) is the physically motivated hedge.

## 7. Exact conditions (planned): satisfied / violated
Satisfied: UEG (w = 0 and the ITYH base reduces to SR-LDA at s = 0 and z = 0; the LR-HF + SR-HF + SR-LDA total is LDA exchange);
exchange spin scaling; C1 exactly in fully polarized one-orbital hydrogenic channels at all lambda; C1 amplitude quadratic for every
one-orbital shape (L1); one-orbital SR exchange exact through O(omega^-4) for every shape (E2); -1/r asymptotics of v_x (LR-HF = 1);
the H atom exactly; E_c[lambda n1] = 0; local LO for the semilocal part (cap).
Violated: the shape error at lambda = 1 for non-hydrogenic one-orbital densities beyond O(omega^-4); GE2 for exchange (base);
flat-plane / two-electron singlets (not addressed); dispersion (no VV10).

## 8. Proof details added in phase 1
- L2 monotonicity: R'(b) = (32/(5pi)) int u^2 e^{-bu^2}(1+u^2)^{-4} du > 0 and R''(b) = -(32/(5pi)) int u^4 e^{-bu^2}(...) du < 0.
  So R is increasing and concave with R(0) = 0. The chord slope Q(b) = R(b)/b of a concave function through the origin is
  non-increasing. Hence sup_b Q = Q(0) = 1/5 and sup F^{SR,1} = (k_s/omega)^2/9.
- L3 uniqueness: amplitude invariance forces F = s^-2 Psi(b) (opus5t/opus5c argument with fixed omega). The hydrogenic family has a
  spatially constant b = Z^2/omega^2, so E_x = Psi(b) * (a b-independent integral). Exactness for every Z then fixes Psi(b) at every
  b > 0: Psi_SR = (5/9) R and Psi_LR = (5/9)(1 - R). b -> 0 at fixed n > 0 happens at every density critical point, so F_LR -> infinity
  there, while F_SR -> (k/omega)^2/9.
- Numerical safeguards, fixed now, before any run: the cap is F_c = F (1 + (F/1.804)^8)^{-1/8}. z_s uses opus5t's regularized form
  (tau_W,s + 0.01 tau_unif,s)/(tau_s + 0.01 tau_unif,s), with a smooth saturation at 1. b = sigma_ss/(4 omega^2 n_s^2) is used
  directly, so no square root of sigma appears and the expression stays complex-step safe. Q(b) is evaluated in closed form through
  erfcx (to be derived symbolically and checked against quadrature in phase 2).
- Gaussian test (to compute): exact SR self-exchange of a normalized Gaussian n = (2a/pi)^{3/2} e^{-2ar^2} is
  E^SR = -sqrt(a/pi) (1 - omega/sqrt(a + omega^2)) (Fourier: n(k) = e^{-k^2/(8a)}). Compare with the branch.

---------------------------------------------------------------------------------------------
# Phase 2 (04:13 onward): implementation, verification, scores

## Implementation facts
- R(b) in closed form. The integration-by-parts recurrence 2n J_{n+1} = (2n-1-2b) J_n + 2b J_{n-1}, with J_1 = (pi/2) erfcx(sqrt b) and
  J_2 = (pi/4)(1-2b) erfcx(sqrt b) + sqrt(pi b)/2, gives J_4 = I(b). It is checked against quadrature to 1e-12 for b <= 50 (work/rfun.py).
  For b < 1e-5 the series Q = 1/5 - b/10 + b^2/6 is used, and for b > 400 the asymptotic series of J_4. Both are complex-step safe.
- PySCF convention verified in the source: RSH = (omega, alpha, beta) means alpha*LR-HF + (alpha+beta)*SR-HF. So RSH = (omega, 1, a-1).
- check_derivs: DERIVATIVES OK (unpolarized and polarized). Audit (submitted K2b, audit.log): DISTINCT (closest single SCAN0, R^2 = 0.902; mix R^2 = 0.963; the first K1 run gave BHANDHLYP 0.916).
  NEG and LO are clean. The audit's UEG check reports a deviation of 0.73 (0.64 for K1). This is an artifact of the audit, which assumes a global hybrid:
  eval_xc returns only the SR semilocal part. At s = 0 and z = 0 it is (1-a) SR-LDA exchange plus PW92, and the HF parts supply the rest.
  For the UEG, LR-HF + a SR-HF equals LR-LDA + a SR-LDA, so the total is LDA.
- omega* (L4) solved exactly: 0.368708 for the cap 1+kappa = 1.804.

## A failure of the first submission form, found in phase 2: variational flattening (and its fix)
K1 (omega* = 0.3687, cap 1.804) is accurate at the exact H density (E_xc = -0.31225, fractional H nearly exact, table below). But the
self-consistent H atom at def2-QZVP came out at -0.50633 Ha, 7.5 kcal/mol too low (work/hdiag.py). The SCF density flattens around the
nucleus: n(0) drops from 0.318 to 0.228, and b drops from the hydrogenic 7.36 to about 3. This raises Q(b) and the branch up to the cap
over a finite region. The minimal-basis SVP H is fine (-0.4995), so the effect only appears in flexible bases. This is the bounded
analogue of opus5c's "charge-linear potential makes the H density non-stationary". Q(b < 1/omega^2) is not fixed by any hydrogenic
density (neutral H has the smallest b of the family). It is an extrapolation, and the 1.804 cap leaves it too much room.
Fix, on the known list: in a one-orbital channel the tight bound is the two-electron/one-orbital bound F_x <= 1.174 (exchange spin
scaling maps a one-orbital channel onto a two-electron singlet density, the class the bound holds for). So the cap is 1.174, and L4
with the same rule gives **omega* = 0.576744**. With this, H at def2-QZVP is -0.500773 (0.49 kcal/mol low), and -0.49944 at SVP.
This change came after a phase-2 run. It is disclosed as candidate K1b. It is a re-derivation with a tighter KNOWN bound, not a
tuned number (L4 is unchanged; only the bound inserted into it changed).

## One-electron tests at exact densities (work/frac.py; E_xc in Ha, f(lambda) = E_xc[lambda n_H]/(-5/16), exact f(1/2) = 0.25,
## f(1/4) = 0.0625; "H2+ inf" = 2E_xc[n_H/2] + 5/32 in kcal/mol, exact 0; Gaussian ratio = E_x/exact at lambda = 1; last column
## = 4E[n_G/2]/E[n_G], exact 1 (amplitude homogeneity on a non-hydrogenic shape))
    func          Exc[H]  f(1/2)  f(1/4)  H2+inf kc |                   | Gauss(a=1) E/exact  He+(Z=2) E
    K0          -0.30646  0.2728  0.0760      -8.96 | Gauss  0.953  He+ -0.60368 (exact -0.62500)  fG(1/2)/fG(1)*4  1.202
    K1          -0.31225  0.2510  0.0630      -0.39 | Gauss  1.140  He+ -0.61571 (exact -0.62500)  fG(1/2)/fG(1)*4  1.116
    K2          -0.31210  0.2512  0.0632      -0.49 | Gauss  1.187  He+ -0.61254 (exact -0.62500)  fG(1/2)/fG(1)*4  1.149
    K0b         -0.31036  0.2611  0.0683      -4.34 | Gauss  0.964  He+ -0.60843 (exact -0.62500)  fG(1/2)/fG(1)*4  1.119
    K1b         -0.31247  0.2506  0.0628      -0.23 | Gauss  1.055  He+ -0.60856 (exact -0.62500)  fG(1/2)/fG(1)*4  1.074
    PBE         -0.31192  0.4125  0.1717     -63.73 | Gauss  0.976  He+ -0.61811 (exact -0.62500)  fG(1/2)/fG(1)*4  1.641
    PBE0        -0.31356  0.3734  0.1449     -48.40 | Gauss  0.987  He+ -0.62139 (exact -0.62500)  fG(1/2)/fG(1)*4  1.478
    B3LYP       -0.31457  0.3867  0.1546     -53.61 | Gauss  0.986  He+ -0.62328 (exact -0.62500)  fG(1/2)/fG(1)*4  1.524
    HYB_GGA_XC_LC_WPBE  -0.31849  0.2957  0.0864     -17.91 | Gauss  0.989  He+ -0.62371 (exact -0.62500)  fG(1/2)/fG(1)*4  1.290
    CAM-B3LYP   -0.31116  0.3293  0.1143     -31.10 | Gauss  0.970  He+ -0.61444 (exact -0.62500)  fG(1/2)/fG(1)*4  1.373
    WB97X       -0.31385  0.3035  0.0952     -20.98 | Gauss  1.003  He+ -0.62316 (exact -0.62500)  fG(1/2)/fG(1)*4  1.332
    R2SCAN      -0.31250  0.3952  0.1559     -56.94 | Gauss  0.996  He+ -0.62500 (exact -0.62500)  fG(1/2)/fG(1)*4  1.582
Reading: every existing functional tested, LC ones included, is convex in lambda (LC-wPBE f(1/2) = 0.296, wB97X 0.304, CAM-B3LYP
0.329, against exact 0.25). The LC base without the branch (K0, K0b) is still off (0.273, 0.261). With the branch (K1, K1b) the
fractional H atom is exact to 0.001 in f, and the stretched-H2+ limit at the exact density is -0.4 / -0.2 kcal/mol, against -18
(LC-wPBE) and -21 (wB97X). The amplitude homogeneity on a Gaussian (last column) is broken only by the cap (1.07 for K1b; without the
cap it would be exactly 1 by L1). The lambda = 1 shape error on a Gaussian is +5.5% for K1b. This is the residual nonlocality error
that no semilocal form can remove (Santra-Perdew / organiser E2).

## Candidates compared (ALL; no continuous parameter fitted to dev data). def2-SVP dev WTMAD, plus QZ where run
| id | omega | cap | SR-HF a | weight | SVP WTMAD (AE/IP MAE) | QZ WTMAD (AE/IP MAE) | H2+ R=6 err | H2+ R=2 err | H+H2 barrier (9.6) | E_SCF(H) QZ |
|---|---|---|---|---|---|---|---|---|---|---|
| K0 (LC base, no branch) | 0.3687 | - | 1/4 | 0 | 2.14 (5.01/5.82) | - | -11.87 | -0.61 | 7.42 | -0.49438 |
| K1 | 0.3687 | 1.804 | 1/4 | zeta^2 z_s^2 | 2.31 (5.30/6.89) | killed (H unstable) | +1.29 | -10.96 | 12.56 | -0.50633 |
| K2 | 0.3687 | 1.804 | 0 | zeta^2 z_s^2 | 3.30 (8.01/6.93) | - | +2.41 | -14.78 | 14.70 | -0.51214 |
| K3 | 0.3687 | 1.804 | 1/4 | z_s^2 (all channels) | 11.12 (27.9/17.7) | - | - | - | - | - |
| K0b (LC base, no branch) | 0.5767 | - | 1/4 | 0 | - | - | -5.35 | -0.44 | 9.15 | -0.49801 |
| K1b | 0.5767 | 1.174 | 1/4 | zeta^2 z_s^2 | 2.90 (6.94/6.76) | 2.27 (5.51/4.91) | +0.04 | -3.74 | 10.38 | -0.50077 |
| **K2b (SUBMITTED)** | 0.5767 | 1.174 | 0 | zeta^2 z_s^2 | **1.73** (3.66/7.21) | **1.21** (2.08/8.27), MCM 5.17 | **+0.09** | -5.09 | **10.42** | -0.50138 |
H2+: aug-cc-pVTZ, E(H2+) - E(H) against UHF in the same basis, in kcal/mol. Barrier: linear symmetric H3 with r = 0.930 A at def2-QZVP, UKS,
grid 4, the harness's run_species (work/barrier.py). opus4 base: H2+ R=6 -28.5, barrier 7.2. opus5t: -6.1, 9.06.
Selection: the pre-registered rule (best SVP WTMAD among branch candidates that converge everywhere and pass H2+ and the barrier).
K1/K2/K3 (omega = 0.3687, cap 1.804) are excluded by the SCF flattening (H atom 4-7.5 kcal/mol too low, which inflates the barrier) and,
for K3, by the expected collapse of closed-shell channels. Between K1b and K2b, K2b wins on SVP (1.73 vs 2.90). The choice of a (1/4
vs 0) is a discrete choice between two physically motivated values (the PBE0 adiabatic-connection fraction vs pure LC). C1 does not
depend on a. The selection used dev scores and is disclosed here. K0 and K0b are references without the new branch, not submission
candidates. K3 confirms the known failure mode (exact one-orbital exchange in closed-shell channels destroys error cancellation).

## Final result (def2-QZVP, submitted K2b): dev WTMAD 1.21, AE MAE 2.08, IP MAE 8.27, MCM 5.17; all 24 species converged (score.json)
Baselines: PBE 2.95, r2SCAN 1.28, B3LYP 0.97, wB97M-V 0.56. opus5t 1.55, opus4 base 1.10, sonnet3 1.00.
- AEs: MAE 2.08, close to B3LYP (1.70) and better than r2SCAN (2.57) and opus5t (3.59). The polarized-atom overshoot of F_q is gone:
  N2 -0.5 (opus5t +7.7), HCN -0.4, F2 -2.2. The largest errors are H2CO +5.4, H2 -5.0 and CH2 +3.4.
- IPs: MAE 8.27, the weak class. IP(F) +10.4 and IP(B) +9.9 (B3LYP +6.8 / +11.5). Both cations are open-shell or closed-shell species
  whose 2p-hole spin density gets the branch and pure LC. I did not try to fix this (it would have meant tuning on the dev set).
- Delocalization: stretched H2+ (R = 6 bohr) error +0.09 kcal/mol, against -28.5 (opus4 base), -6.1 (opus5t), -5.4 (LC base K0b).
  At R = 2 the error is -5.1, the shape error of the hydrogenic branch at the bonded H2+ midpoint (b -> 0, where the on-top extrapolation
  overestimates). The H+H2 barrier is 10.42 (ref 9.6, +0.8, against 7.2 for the base). About 0.9 kcal/mol of the overshoot comes from the
  SCF H atom being 0.87 kcal/mol too low (E(H) = -0.501380).

## Exact conditions: satisfied / violated (submitted)
Satisfied: UEG (total; the audit's 0.73 "UEG deviation" is its global-hybrid assumption, explained above); exchange spin scaling; e_xc <= 0
(NEG 0); local LO (0 violations); the one-orbital bound F <= 1.174 for the branch; C1 amplitude homogeneity (degree 2) exactly in
fully polarized one-orbital channels wherever the cap is inactive, for every shape (L1); C1 exactly for every fractional hydrogenic
density with Z <= 1 at the exact density (f(1/2) = 0.2506); interfragment exchange at long range exact (LR-HF = 1); -1/r asymptotics of
v_x; one-orbital SR exchange exact through O(omega^-4) for every shape (E2); E_c[lambda n1] = 0 (Stoll gate, prior art);
H atom within 0.9 kcal/mol self-consistently.
Violated: lambda = 1 shape exactness for non-hydrogenic one-orbital densities (impossible for any semilocal form; Gaussian +5.5%, H2+ R=2
-5 kcal/mol); exactness for Z > 1 hydrogenic ions (the cap is active at the nucleus: He+ off by 0.017 Ha); exchange GE2 in the base;
flat-plane / fractional spin (not addressed); dispersion (no VV10); the SCF H density is not exactly stationary at the exact density
(slight flattening, -0.87 kcal/mol).

## What is new (precise), versus KNOWN_CONDITIONS.md, prior_art and the round-5 literature check
1. L2: range separation evades the Laplacian-free no-go for the AMPLITUDE part. The one-orbital SR branch F = (5/9)(k_s/omega)^2 Q(b),
   with b = |grad n_s|^2/(4 omega^2 n_s^2), is exactly amplitude-quadratic for every shape and is bounded at every density. Its supremum
   is the exact on-top SR exchange -(pi/(2 omega^2)) n_s^2, reached at critical points. (The function R was in opus5c's Corollary 4.
   That it is bounded, the on-top supremum, and the use of b were not.)
2. L3: C1 plus pointwise finiteness forces 100% long-range exact exchange in one-orbital channels. The LR semilocal complement diverges at
   every critical point. This is a derivation of alpha_LR = 1 from fractional-charge exactness, not from the -1/r asymptote.
3. E2 with the hydrogenic check: the O(omega^-4) term of one-orbital SR exchange is Laplacian-free after integration by parts, with
   coefficient -b/2. The hydrogenic branch reproduces it, so it is shape-exact through second order.
4. L4: a non-empirical range parameter. omega* is the smallest omega for which the C1-exact branch respects the one-orbital bound on the
   H atom: omega* = 0.5767 for F <= 1.174 (0.3687 for F <= 1.804, which proved SCF-unstable). It is fixed by the one-electron reference
   alone and does not depend on the SR exact-exchange fraction.
5. The construction: the first bounded, Laplacian-free functional in this benchmark with exact fractional-H behaviour (f(1/2) = 0.2506)
   and near-exact stretched H2+ (+0.09 kcal/mol at 6 bohr) that keeps first-row AEs near B3LYP level (MAE 2.08).
6. A documented negative result: the variational flattening instability of C1-exact branches in flexible basis sets, and that the tight
   one-orbital bound (1.174) controls it where the LO bound (1.804) does not.
Not new: C1 itself (Zhang-Yang 1998; MCY 2006), LC hybrids and ITYH, the Stoll correlation gates, opus4's correlation, sonnet3's
exchange, zeta^2 z_s^2 weighting (opus5t), the 1.174 and 1.804 bounds.

## Hidden-class expectations
Barrier heights: good (LC plus delocalization fix; H+H2 10.4). Proton affinities and noncovalent: closed-shell, so pure LC base. NCI has no
dispersion term and may be weak. Second-row AEs: the polarized atoms P, S, Cl get the branch only in one-orbital-like regions. Risk:
pure LC with omega = 0.58 underbinding heavier molecules. IPs are the weakest class.

## With more time
- Constrain Q(b) for b < 1/omega^2 (below the hydrogenic range) with an exact-hole argument instead of the cap. This should remove the SCF
  flattening and the H2+ R = 2 error, and might allow the smaller omega* = 0.369.
- Examine the IP(F) and IP(B) errors (2p-hole species) and apply the branch only to the C1 premise (n_s' -> 0) rather than zeta^2.
- VV10 with C from an exact C6 for the NCI class; second-row AEs.
