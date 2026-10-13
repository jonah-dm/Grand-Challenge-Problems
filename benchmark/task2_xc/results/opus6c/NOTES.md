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
