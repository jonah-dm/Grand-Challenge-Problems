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

## 5. Phase-2 verification plan
V1. Exact Lieb functional, computed by me, for a 1D soft-Coulomb one-orbital density along
    (f, f), f in (0, 1], and off-diagonal (f_a, f_b) points. Method: concave maximisation over v of
    E_v(f_a,f_b) - int v n, with an LP convex hull of the FCI corners N <= 2. This checks C2(i) equality, the C2(ii)
    sandwich and the tent/kink.
V2. 3D exact H-shaped and Gaussian orbitals. Compare E_c and E_xc of LDA, PBE, TPSS, SCAN, r2SCAN, B3LYP,
    PBE0, M06-L and the base/candidates with C2 on the square.
V3. Derivative check, audit, SVP dev scores for A/B/C, QZ for the chosen one, H2+ at R = 6 bohr, H+H2 barrier.
