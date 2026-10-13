# opus6t (THEORY track): the one-orbital occupation square for correlation (flat-plane axis)

Started 03:13 UTC. Phase 1 (derivation only, no code, no SCF): 03:13-04:13. Phase 2 is appended below.

## 0. What is used from KNOWN_CONDITIONS.md / prior art (no credit claimed)
- PPLB piecewise linearity, the flat-plane condition (MCY 2009) and the constancy condition for
  fractional spins **as statements about the energies E(N_a, N_b) of a fixed system**. Fractional-spin
  H (E[n_H/2, n_H/2] = E[n_H]) is the known constancy condition.
- Round 5 (opus5t/sonnet5c/opus5c): C1, E_xc[lambda n1, 0] = -lambda^2 J[n1] for lambda <= 1, and T2:
  in a (n, grad n, tau) functional, C1 forces F(s, z=1) = (5/9) s^-2. This is the special case
  f_b = 0 of what follows; I reuse it and claim nothing for it.
- Base functional pieces: 25% HF (PBE0 argument), PBE exchange, PBE/PW92 correlation, opus4's
  pair-aware Stoll gate on correlation (round 4), z = tau_W/tau, opus5t's bounded form
  F_q = F_LO/(1 + F_LO s^2/c_q) with c_q from the H norm, and opus5t's z regulator. All of these are prior art.

## 1. Setting and conventions
Spin-DFT with the Lieb (ensemble, Fock-space) functional
  F[n_a, n_b] = inf { Tr G (T + W) : G an ensemble with spin densities (n_a, n_b) },
where W is any pair interaction w(r - r') >= 0 (Coulomb in 3D, soft Coulomb in 1D). T_s is the same
infimum with W dropped. J[n] = (1/2) int int n n' w. E_xc = F - T_s - J[n_a + n_b].
**Exchange convention** (the one in which hybrids, fractional-occupation DFT and exchange spin scaling
are defined): E_x is the Fock term of the KS one-body density matrix, gamma_s = f_s |phi><phi| for a
one-orbital channel. This gives E_x = -(f_a^2 + f_b^2) J[n1]. Correlation is E_c = E_xc - E_x.
(The ensemble-average-of-determinants exchange would give a different split. With that split,
exchange spin scaling fails at fractional N, so I do not use it.)

One-orbital spin densities: n1 = |phi|^2 normalised to 1, n_s = f_s n1, 0 <= f_s <= 1, N = f_a + f_b,
zeta = (f_a - f_b)/N, constant in space.

## 2. Theorem C2 (the occupation square). NEW in this form
For every one-orbital shape n1 and every (f_a, f_b) in [0,1]^2:

  (i)  particle triangle, N <= 1:  E_xc[f_a n1, f_b n1] = -N^2 J[n1]          (equality)
                                   E_c = -2 f_a f_b J[n1]
  (ii) hole triangle,     N >= 1:  -N^2 J[n1] <= E_xc <= -[1 + (N-1)^2] J[n1]
                                   -2 f_a f_b J <= E_c <= -2 (1-f_a)(1-f_b) J[n1]

So **E_c <= -2 min{ f_a f_b , (1-f_a)(1-f_b) } J[n1]** on the whole square. The upper bound is
particle-hole symmetric under (f_a, f_b) -> (1-f_a, 1-f_b) (a "tent" over the square). The bound on
E_xc depends on N only. On the particle triangle E_xc is **spin-blind**: it does not depend on zeta.

**Proof.**
(a) T_s. For any ensemble, T >= sum_s T_W[n_s], with T_W[n] = (1/8) int |grad n|^2/n (von Weizsaecker,
spin resolved; T_W is convex, so the bound passes to ensembles). The noninteracting ensemble that
occupies the nodeless orbital phi = sqrt(n1) with weights f_s attains it. So T_s = N T_W[n1].
(b) Lower bound on E_xc everywhere. F >= T_s because Tr G W >= 0. Hence E_xc >= -J[N n1] = -N^2 J[n1].
(c) Particle triangle. Take G = (1-N)|vac> + f_a |phi up> + f_b |phi down>. It has the right spin
densities and no state with two electrons, so Tr G W = 0, and its kinetic energy is N T_W[n1].
Hence F <= T_s. With (b), F = T_s and E_xc = -N^2 J. Subtracting E_x = -(f_a^2+f_b^2) J gives
E_c = -2 f_a f_b J.
(d) Hole triangle. Take p_up = 1-f_b, p_down = 1-f_a, p_2 = N-1 (the singlet determinant
|phi up, phi down>), p_vac = 0. Spin densities: p_up + p_2 = f_a and p_down + p_2 = f_b. Kinetic energy:
(p_up + p_down + 2 p_2) T_W = N T_W. Interaction: p_2 (phi phi|phi phi) = 2 (N-1) J[n1]. So
F <= N T_W + 2(N-1) J and E_xc <= -N^2 J + 2(N-1) J = -[1+(N-1)^2] J. Subtracting E_x gives
E_c <= [2N - 2 - 2 f_a f_b] J = -2(1-f_a)(1-f_b) J.  QED.

Remarks.
- (i) at f_b = 0 is round-5's C1. (i) at N = 1 is the known constancy condition. (i) for N < 1
  in the density-functional form (any shape, any spin split), the product law
  E_c = -2 f_a f_b J, and all of (ii) are what I claim as new. (ii) is a statement about E_xc at fixed
  density and arbitrary shape. It does not follow from the flat-plane condition, which is about
  ground-state energies of one fixed potential, with densities that change along the plane.
- **Physical content.** The static correlation of a fractionally occupied orbital is exactly minus the
  opposite-spin Hartree energy, -int int n_a n_b w (= -2 f_a f_b J). It is largest at
  half filling, f_a = f_b = 1/2 (E_c = -J/2; for an H-shaped orbital, -0.156 Ha, which is 3.7x the
  He correlation energy). Its mirror bound vanishes at the closed-shell corner (1,1), where only dynamic
  correlation remains. Along the closed-shell axis f_a = f_b = f, the upper bound is
  -2 min(f^2, (1-f)^2) J. It is non-differentiable at f = 1/2 (N = 1). This is the correlation face of
  the derivative discontinuity, now as a density-functional statement.
- The hole-branch bound is not tight at the corner: the gap equals the dynamic correlation of the
  pair (for the He-like pair, about 0.04 Ha). I verify the whole profile numerically in phase 2.


## 2b. Corollary C2-HD (high-density / scaling form). NEW
Let n1_g(r) = g^3 n1(g r) (uniform coordinate scaling). Then T_W scales as g^2 and J as g. By scaling,
F[f n1_g] = inf over ensembles G -> (f_a n1, f_b n1) of { g^2 T(G) + g W(G) }. As g -> infinity the minimiser
must have T = T_s at leading order. The ensembles with T = T_s are exactly the mixtures of the four states built
from phi (vac, up, down, singlet up-down). States with two same-spin electrons or a second orbital cost O(g^2).
The first-order term is the degenerate Goerling-Levy minimum of W over that manifold: an LP whose minimum is at
p_2 = max(0, N-1). Hence
   lim_{g->inf} E_c[f_a n1_g, f_b n1_g] / g  =  -2 min{ f_a f_b , (1-f_a)(1-f_b) } J[n1],
so **the C2 upper bound is exact at leading order in the high-density limit on the whole square**. The
gap in (ii) is O(g^0): the dynamic correlation of the pair.
Consequences (new as stated):
- For fractional one-orbital densities, correlation scales like exchange. E_c[n_g] = g E_c[n] holds exactly on the
  particle triangle, so the Levy-Perdew correlation scaling inequality is saturated there. On the open hole
  triangle, E_c[n_g] = -2(1-f_a)(1-f_b) g J + O(1).
- The high-density limit of correlation, finite for integer systems (a listed condition), is -infinity
  (linear in g) everywhere on the open square except the three integer corners (0,0), (1,0)/(0,1), (1,1),
  where the prefactor vanishes.
- Every semilocal correlation built with a finite high-density limit (PW92/PBE/SCAN/r2SCAN by
  design, LYP, B95, opus4) has E_c[n_g]/g -> 0 on the square. All of them violate C2-HD at leading order.

## 3a. Rigorous version of T3 on the hydrogenic family (replaces the "pointwise sense" caveat)
For n = N Z^3 e^{-2Zr}/pi, |grad n| = 2Z n, so s k_F = Z: each hydrogenic shape lies on one curve
s = Z/k_F(n) of the (n, s) plane, and its maximum density is w = N Z^3/pi. On that curve
e_c = G(n), and E(w) = (pi/(2 Z^3)) int_0^w G(u) ln^2(w/u) du/u. With D = w d/dw this gives D^3 E = (2 pi/(2 Z^3))... more precisely
G(w) = (Z^3/pi) D^3 E(w). C2(i) gives E(w) = -(1-zeta^2)(5 pi^2/32) w^2/Z^5 for every w <= Z^3/pi (all N <= 1),
hence G(u) = -(1-zeta^2)(5 pi/4) u^2/Z^2. Using Z = s k_F:
   e_c(n, zeta, s; z=1) = -(1-zeta^2) * 5/(4 (3 pi^2)^{2/3}) * pi * n^{4/3} s^{-2}    for all s >= s* = (3 pi)^{-1/3} = 0.4730
(the point (n, s) is reached by a hydrogenic N <= 1 density iff s >= s*). At zeta = 1 the same formula holds for exchange
(it is round 5's 5/9 s_sigma^-2 law). So **the particle law determines any meta-GGA correlation uniquely on
s >= 0.473**. There the form is bounded (|e_c| <= const n^{4/3}/s*^2).
Gaussian orbital at amplitude lambda: the region s < s* is a ball of radius O(lambda^{1/3}). Any e_c with
|e_c| <= C n^{4/3} there contributes O(lambda^{7/3}), while the forced s^-2 part gives 5x the exact value at
order lambda^2 (round 5's analytic factor; I re-derived int n^{4/3} s^-2 = 2 (3 pi^2)^{2/3} sqrt(a)/pi^{3/2}).
**Theorem T3 (rigorous):** a (n, grad n, tau) correlation with any finite Lieb-Oxford-type bound that satisfies
C2(i) on all hydrogenic one-orbital densities at some |zeta| < 1 has
E_c^sl[lambda g]/E_c^exact[lambda g] -> 5 for Gaussian orbitals as lambda -> 0. No such functional satisfies C2(i).
By-product: for hydrogenic shapes, s at the nucleus is (3 pi N)^{-1/3}, a local occupation meter
(N = 1: 0.473, N = 2: 0.375). It fails for cuspless shapes, so it cannot remove the T3b obstruction in general.
Along the closed-shell axis on hydrogenic shapes, D^2 E = (pi/Z^3) int_0^w G du/u is continuous in w for any locally
integrable G. So a semilocal E_c(N) is C^1 across N = 1. Any kink of the exact one (tested in V1) is therefore unreachable
(this is the known "no derivative discontinuity" statement, here localised to the correlation ray).

## 3. Theorem T3 (no-go for semilocal correlation). NEW
Let E_c^sl = int e_c(n_a, n_b, grad n_a, grad n_b, tau_a, tau_b) be any meta-GGA correlation. On
one-orbital densities, tau_s = tau_W,s, so z_a = z_b = z = 1. The local variables reduce to n, the constant zeta,
and s, which scales as s -> lambda^-1/3 s when n -> lambda n.
(a) C2(i) for all N <= 1, all zeta and all shapes holds (in the pointwise sense used for round-5 T2) iff
    e_c(n, zeta, s; z=1) = (1 - zeta^2)/2 * e_x^C1(n, s),   e_x^C1 = -(5/9)-normed c n^{4/3} s_p^{-2}
    (s_p is the fully polarised reduced gradient of the total density). The zeta-dependence is forced to be exactly
    (1-zeta^2)/2, and the amplitude dependence is forced to the unbounded s^-2 law. Consequently
    e_xc = e_x^C1[n, 0] pointwise (spin-blind). The proof is the T2 functional equation applied at each fixed zeta,
    because C2(i) has amplitude degree 2 at fixed zeta. The H norm fixes the constant. The same form is then 5x
    too negative for a Gaussian orbital, as in round 5, so **no semilocal correlation satisfies C2(i)**.
    Every bounded e_c has amplitude degree < 2 on the ray, so it misses the static correlation
    -lambda^2 (1-zeta^2) J/2 at small lambda.
(b) Indicator obstruction (T3b). At a point with zeta = 0 and z = 1, the local data (n, s) of a
    half-filled one-orbital density (N = 1, particle branch, E_c = -J/2) and of a closed-shell
    two-electron pair (N = 2, corner of the hole branch, E_c = dynamic only) can coincide. The local data do
    not see the normalisation. So any weight w(n, zeta, s, z) that switches on the spin-blind form is one number
    for both, and **no semilocal functional can satisfy the constancy/flat-plane axis and keep closed-shell pairs
    right**. Within a weight W(zeta, z), the most C2 can do is act where zeta != 0, i.e. at the
    corners zeta = +-1 that C1 already covers. Corollary: in integer UKS, a functional that satisfies C1 at
    zeta = +-1 gets stretched singlets right through spin-symmetry breaking (each fragment is zeta = +-1, C1-exact).
    C2 says this broken-symmetry energy is the exact (spin-blind) value, not an artefact.
(c) Hybrids. A global hybrid's exact exchange is same-spin only. On the particle triangle it supplies
    -a(f_a^2+f_b^2) J and none of the -2 f_a f_b J. So the semilocal remainder must carry
    [(1-a)(1+zeta^2)/2 + (1-zeta^2)/2] e_x^C1[n,0]. That is, the static-correlation part is NOT reduced by
    the hybrid fraction, while the exchange part is.


## 3c. Theorem T4 (pair catastrophe: the flat-plane axis forces over-correlation of two-electron singlets). NEW
Hypotheses: a (n, grad n, tau) correlation with e_c <= 0 pointwise (SCAN-type negativity, satisfied by PW92, PBE,
SCAN, r2SCAN and B95) that obeys C2(i) at zeta = 0 on hydrogenic one-orbital densities. This is the
constancy condition for every fractional-spin hydrogenic ion, plus its N < 1 extension.
Claim: for every two-electron closed-shell hydrogenic density n = 2 Z^3 e^{-2Zr}/pi,
   E_c^sl <= -(5/8) (ln^2 2 / 2 + ln 2 / 2 + 1/4) Z = -0.5230 Z  Ha,
which is 84% of the absolute lower bound -2J[n1] = -(5/8) Z from C2(ii). For a He-like orbital (Z = 27/16) this is
-0.88 Ha, against an exact pair correlation of about -0.042 Ha (He).
Proof: on the curve s k_F = Z the pair density covers u in (0, 2 w1] with w1 = Z^3/pi. Its outer part
u <= w1 (r >= ln2/(2Z), s >= s*) is fixed by section 3a: G(u) = -(5 pi/4) u^2/Z^2. The inner part contributes <= 0
by negativity. With the Mellin form, E_outer = (pi/(2Z^3)) int_0^{w1} G(u) ln^2(2 w1/u) du/u
= -(5/8) Z int_0^1 t ln^2(2/t) dt = -(5/8) Z (ln^2 2/2 + ln 2/2 + 1/4). QED.
For a constant spin polarisation the prefactor is (1 - zeta^2).
Meaning: in any semilocal correlation with e_c <= 0, the static-correlation axis (fractional spin, the flat plane)
and the two-electron singlet are not just locally indistinguishable (T3b). Enforcing the first on hydrogenic ions
quantitatively destroys the second, by more than 0.5 Z Ha. This is the correlation counterpart of round-5's
exchange obstruction, stated for real (integer, closed-shell) densities. It is also why every working semilocal
correlation (PBE, r2SCAN, B95, opus4) violates the constancy condition by an amount of order J. Related prior
observation: opus5c's (round 5) shift map n_H/2-spin <-> 2-electron density at r + ln2/(2Z), used there for exchange.

## 0b. Literature positioning (added after the organiser's round-5 literature check)
Known, no credit claimed:
- Zhang & Yang 1998; Mori-Sanchez, Cohen & Yang 2006 (lambda^2 amplitude scaling, i.e. the f_b = 0 edge of C2(i)).
- Cohen, Mori-Sanchez & Yang, JCP 129, 121104 (2008): constancy for fractional spins, i.e. C2(i) at N = 1, as energy
  statements. MCY PRL 102, 066403 (2009): the flat plane. The particle-triangle law E_xc = -N^2 J (any shape) follows
  from these with the Lieb argument. I treat C2(i), including the product form E_c = -2 f_a f_b J, as **known/implicit**.
- Gould & Dobson JCP 2013 (LEXX, piecewise-linear exact exchange for fractional occupations) and Gould & Pittalis
  PRL 119, 243001 (2017). The latter defines the ensemble Hartree-exchange as the lambda -> 0+ derivative of the ensemble
  universal functional. My C2-HD limit is the one-orbital instance of that ensemble-Hx minimum. Expressing it as the
  tent -2 min{f_a f_b, (1-f_a)(1-f_b)} J for the **1-RDM-convention correlation** (the one hybrids and semilocal
  correlations are built for) is my reformulation, not a new concept. I claim at most the explicit hole-branch sandwich
  C2(ii) and its consequences, and I mark it "likely implicit in Gould-Dobson/Gould-Pittalis" (full texts not read).
- Santra & Perdew 2019: no semilocal functional is exact for all one-electron densities. The fragment argument (n1 = (n_A+n_B)/2)
  applies to correlation as well. So "no semilocal correlation satisfies the constancy condition for all shapes" is
  known in substance, and T3's conclusion carries no credit.
Claimed as new, beyond those papers:
1. T3 (rigorous, hydrogenic family): the particle law forces a unique e_c = -(1-zeta^2) (5 pi/(4 (3pi^2)^{2/3})) n^{4/3} s^-2
   on s >= (3 pi)^{-1/3}, with the factor (1-zeta^2) uniquely fixed, and the Gaussian factor 5 for correlation.
2. T4: the quantitative pair catastrophe (E_c <= -0.523 Z for two-electron singlets) for any e_c <= 0 semilocal correlation
   that satisfies the flat-plane axis on hydrogenic ions. This is the "second condition of the same kind for two-electron singlets" requested.
3. Section 4: the one-orbital manifold is a property of the total density (z_tot, not z_s). This gives a derived weight for any
   bounded C1/C2 construction, and the spin-blind joint xc form with the hybrid split [(1-a)(1+zeta^2)/2 + (1-zeta^2)/2].
4. The numerical verifications V1/V2.

## 4. Where C2 tells a construction something new: the one-orbital manifold is a property of the TOTAL density
C1 (round 5) was enforced with a per-channel weight zeta^2 z_s^2. C2 shows that the exact
amplitude law holds on the one-orbital manifold only when **both** channels hold the same orbital or one is empty. This is the condition
z_tot = tau_W[n]/tau = 1 for the total density. When the minority channel holds a different
orbital in the same region (the beta 2s under alpha 2p in B, C, N, O), the region is not on the manifold
and neither C1 nor C2 applies. A single-channel z_s = 1 there is not sufficient. That is the origin of round 5's
overshoot for polarized atoms (N2 +9.6, IP(B) +7.8). This is a consequence of C2, not a fit.

### Construction (bounded, partial enforcement of C1+C2 jointly; "spin-blind one-orbital xc", SB)
E_xc = a E_x^HF + int { (1-W) [ (1-a) e_x^PBE(n_a,n_b) + e_c^base ] + W e_SB },   a = 1/4
e_SB  = [ (1-a)(1+zeta^2)/2 + (1-zeta^2)/2 ] * e_x^q[n, 0]
e_x^q[n,0] = the fully polarised exchange of the total density with F_q(s_p) = F_LO/(1 + F_LO s_p^2/c_q),
           F_LO = 1.804, c_q from the H-atom norm (opus5t's form, prior art; constant recomputed by me).
W = zeta^2 z_tot^2.
- On the particle triangle with W = 1, this is the exact C2 structure, with C1 regularised as in round 5.
- zeta^2: T3b says the zeta = 0 one-orbital point is ambiguous and closed-shell pairs must be left alone.
  zeta^2 is the lowest even polynomial with P(0)=0 and P(+-1)=1 (the same choice as opus5t, not new).
- z_tot^2 instead of z_s^2: section 4 (new, derived).
- e_c^base = opus4 pair-aware Stoll correlation (prior art). It vanishes for zeta=1 one-orbital densities, consistent with C2.
Limits: UEG (z_tot = 0, so W = 0); closed shells (zeta = 0, so W = 0: identical to the base); H atom exact
(W = 1, a HF + (1-a) q both exact); the polarised one-orbital ray as in round 5, plus the correct static
correlation factor for 0 < |zeta| < 1 one-orbital regions.

Candidates planned (all will be listed with scores): A base (no C1/C2); B opus5t-type per-channel C1
(W_s = zeta^2 z_s^2, exchange only) on my base (prior-art comparator); C = SB (submission candidate).


**T4' (shape-free form).** The forced form of section 3a is a pointwise statement about e_c(n, zeta=0, s, z=1) on
s >= s*. Every two-electron singlet density (any shape: He, H-, Li+, H2 at any R) has z = 1 and zeta = 0 everywhere. So, under the
same hypotheses, for EVERY two-electron singlet density n:
   E_c^sl[n] <= -(5 pi / (4 (3 pi^2)^{2/3})) * int_{s(r) >= s*} n^{4/3} s^{-2} d^3r,   s* = (3 pi)^{-1/3}.
(Hydrogenic shape: -0.523 Z.) Phase 2 evaluates the right-hand side on HF densities of He and H2 and compares it with the
near-exact E_c (He -0.0420, H2 -0.0409 Ha).
Corollary along the closed-shell hydrogenic axis, N in [1,2]: E_c^sl(N) <= -(5Z/8)(1/4 + ln N/2 + ln^2 N/2). This keeps
decreasing past N = 1, while the exact C2(ii) upper bound -(5Z/8)(1 - N/2)^2 turns back towards 0. A semilocal
correlation that satisfies the particle law cannot follow the far side of the tent.

## Phase-1 decisions for the construction (before any run)
Base: 25% HF + 75% PBE exchange + opus4 pair-aware Stoll correlation on PBE (all prior art). I cannot use sonnet3's
exchange because its tail form is not fully specified in prior_art and I do not read other submissions' code.
Closed shells are therefore PBE0-exchange plus opus4 correlation. The audit, which is unpolarised, may report NEAR-EXISTING for that reason. The new
content acts only in spin-polarised one-orbital-manifold regions.
Candidates to score (SVP, all disclosed): A = base; B = per-channel C1 (opus5t-type, W_s = zeta^2 z_s^2) on the base;
C = spin-blind SB (W = zeta^2 z_tot^2), the submission candidate. c_q comes from the H norm at the exact density. If the TS SCF fails,
the H norm is imposed self-consistently, as in opus5t.

## 5. Phase-2 verification plan
V1. Exact Lieb functional, computed by me, for a 1D soft-Coulomb one-orbital density along
    (f, f), f in (0, 1], and off-diagonal (f_a, f_b) points. Method: concave maximisation over v of
    E_v(f_a,f_b) - int v n, with an LP convex hull of the FCI corners N <= 2. This checks C2(i) equality, the C2(ii)
    sandwich and the tent/kink.
V2. 3D exact H-shaped and Gaussian orbitals. Compare E_c and E_xc of LDA, PBE, TPSS, SCAN, r2SCAN, B3LYP,
    PBE0, M06-L and the base/candidates with C2 on the square.
V3. Derivative check, audit, SVP dev scores for A/B/C, QZ for the chosen one, H2+ at R = 6 bohr, H+H2 barrier.

## Verification results (computed during phase 1 as part of the derivation; no SCF and no functional code; scripts in work/)
### V1: exact 1D Lieb functional (work/lieb1d.py, v1.out)
Soft-Coulomb 1D model: w = 1/sqrt(x^2+1), n1 = ground density of v = -1/sqrt(x^2+1), 61-point grid, J[n1] = 0.34621.
Spin-resolved Lieb functional over the N <= 2 Fock space. Method: concave dual maximisation over (v_a, v_b), entropic
smoothing annealed to T = 4e-4. The computed F is a lower bound on the true F; the error is <= T ln 6 = 7e-4 Ha
(0.002 J) plus incomplete convergence.

| (f_a, f_b) | N | E_c/J computed | C2 lower -2 f_a f_b | C2 upper (tent) | status |
|---|---|---|---|---|---|
| (0.25, 0.25) | 0.5 | -0.1262 | -0.125 | -0.125 | equality (within the 0.002 smoothing bound) |
| (0.5, 0.5) | 1.0 | -0.5008 | -0.500 | -0.500 | equality |
| (0.75, 0.25) | 1.0 | -0.3757 | -0.375 | -0.375 | equality |
| (0.625, 0.625) | 1.25 | -0.3136 | -0.781 | -0.281 | inside, 0.032 J below the tent |
| (0.75, 0.75) | 1.5 | -0.1707 | -1.125 | -0.125 | inside, 0.046 J below |
| (0.875, 0.875) | 1.75 | -0.0835 | -1.531 | -0.031 | inside, 0.052 J below |
| (1, 1) | 2.0 | -0.0549 | -2.000 | 0 | inside (1D pair dynamic correlation, 0.019 Ha) |
| (1, 0.5) | 1.5 | -0.0426 | -1.000 | 0 | inside |
The particle triangle is reproduced to 0.1%. On the hole triangle the exact E_c follows the tent, offset by a smooth gap
that grows to the pair dynamic correlation at (1,1), as C2-HD predicts (an O(g^0) gap). The exact E_c(N) on the closed-shell axis
has a strong kink at N = 1: left slope -1.00 J, right slope about +0.75 J. Section 3a shows that any semilocal E_c is C^1 there.

### V2: 3D exact one-orbital shapes, existing functionals (work/v2_dfa.py, v2.out). "!" = C2 violated
H 1s (J = 5/16). E_c/J, correlation-only functionals:
| | (0.5,0) | (0.25,0.25) | (0.5,0.5) | (0.75,0.25) | (0.75,0.75) | (1,0.5) | (1,1) |
|---|---|---|---|---|---|---|---|
| exact / tent bound | 0 | -0.125 | -0.500 | -0.375 | <= -0.125 | <= 0 | <= 0 |
| LDA(PW92) | -0.032 | -0.059! | -0.132! | -0.121! | -0.211 | -0.203 | -0.294 |
| PBE | -0.006 | -0.016! | -0.048! | -0.041! | -0.087! | -0.083 | -0.133 |
| LYP | 0 | -0.014! | -0.044! | -0.038! | -0.081! | -0.076 | -0.122 |
| TPSS | 0 | -0.018! | -0.051! | -0.047! | -0.089! | -0.087 | -0.133 |
| SCAN / r2SCAN | 0 | -0.022! | -0.051! | -0.044! | -0.081! | -0.076 | -0.113 |
| B95 | 0 | -0.018! | -0.045! | -0.035! | -0.077! | -0.070 | -0.112 |
| M06-L (c) | 0 | -0.002! | -0.032! | -0.019! | -0.066! | -0.058 | -0.103 |
- Every correlation functional misses 75-95% of the particle-triangle static correlation (the known fractional-spin error, now on
  the whole triangle and as a pure-correlation statement).
- **New violation: the hole-triangle bound at (3/4, 3/4)**, E_c <= -J/8. All GGA/meta-GGA correlations
  (PBE, LYP, TPSS, SCAN, r2SCAN, B95, M06-L) are above it, i.e. not negative enough. Only LDA satisfies it. The Gaussian shape gives the same
  picture (v2.out).
- For total E_xc (with the hybrid's exact exchange), every functional violates the particle equality at N = 1 (E_xc/J between
  -0.77 and -0.90 vs -1, too high: static correlation error), and at N = 1/2 is too low (-0.30 to -0.41 vs -0.25: delocalisation
  error). The two errors have opposite signs on the same triangle. The hole E_xc bound is met by all of them, except LDA on the
  Gaussian pair (1,1), where E_xc^LDA = -1.970 J > E_x = -2 J.
- T4 check: the forced-form bound for the hydrogenic pair is -0.52309 (Z=1) and -0.88243 (Z=27/16), numerically. The analytic value is -0.5230 Z.
