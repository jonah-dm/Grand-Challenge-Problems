# opus5c — CONSTRUCTION track
# "Charge-linear one-electron exchange": a semilocal exchange that is exact on the whole
# fractional-charge line 0 <= N <= 1 for every hydrogenic one-electron system

(Phase 1 = derivation, written 00:40–01:40 before any SCF run. Phase-2 results are appended below.)

## 0. Summary of the claim (what is new, stated precisely)

Known condition used (on KNOWN_CONDITIONS.md): piecewise linearity in fractional electron
number (PPLB 1982) and E_xc[n_1] = -J[n_1] for one-electron densities. For a one-electron system
on the segment 0 <= N <= 1, the exact ensemble density is q*n_1 (q = N) and the exact functional
must give E_xc[q n_1] = -q^2 J[n_1] for every q in (0,1] (derivation in Section 1).

No semilocal or global-hybrid functional satisfies this, even for the H atom alone. All of them
give E(q) convex (delocalization error), including the ones that are exact for the H atom at
q = 1 (SCAN, opus3/haiku3-style H norms). Local hybrids with a(r) = 1 in iso-orbital regions
do satisfy it, but they use the exact-exchange energy density. That is a listed known
construction, and this harness cannot supply it.

New results (none of them is on the known list):

**Theorem 1 (uniqueness).** Take any meta-GGA built from (n_s, grad n_s, tau_s) that removes
correlation for fully polarized one-orbital densities (the known gate). In that limit
tau = tau_W, so its exchange reduces to an enhancement F_1(s), with
s = |grad n_s| / (2 (6 pi^2)^{1/3} n_s^{4/3}). The functional is exact on the whole charge line
0 < q <= 1 for the H atom if and only if

    F_1(s) = A / s^2   for all  s >= s0 = (6 pi)^{-1/3} = 0.37575,   with A = 5/9 exactly.

Below s0, F_1 is unconstrained.

**Corollary 1 (closed form).** On the hydrogenic range, the one-orbital exchange is the
parameter-free functional

    E_x^{1}[n_s] = -10 pi ∫ n_s^4 / |grad n_s|^2 d^3r .

It is homogeneous of degree 2 under density scaling n -> q n, which is exactly the scaling of
the true one-orbital exchange -J[n_s] = -q^2 J. It is homogeneous of degree 1 under coordinate
scaling, as all exchange must be. It is exact, point by point along the charge line, for every
hydrogenic ion (all Z, all 0 < q <= 1). This is a two-parameter family of exact one-electron
results from one constant, not a single H norm.

**Corollary 2 (no-go for pointwise bounds).** F_1(s0) = (5/9)(6 pi)^{2/3} = 3.937. So every
semilocal functional that imposes a pointwise iso-orbital bound F_x <= 1.174 (SCAN, the
two-electron bound) or F_x <= 1.804 (PBE/TPSS) necessarily has a convex E(q) for H. Imposing
those bounds pointwise is incompatible with piecewise linearity for one-electron systems. The
integrated Lieb–Oxford bound is still satisfied: E_x[H] = -0.3125 >= -1.68 ∫ n^{4/3} = -0.484.
In the same way, any F_1 that grows at large s (B88's s/ln s) or saturates (PBE) is excluded.
Hence **no existing meta-GGA or GGA, hybridized globally or not, can be charge-linear for H.**

**Corollary 3 (one-orbital channels in general).** For any one-orbital spin channel
n_s = q|phi|^2, exact exchange is -q^2 J[|phi|^2], which is homogeneous of degree 2 in n_s
whatever the shape of phi. The only s-dependent enhancement with that homogeneity is A/s^2.
So a semilocal exchange using F_1 has E_x^SL[q n] = q^2 E_x^SL[n] in every one-orbital
channel. Then the fractional-charge curvature of such a channel is
2(1-a)(J + E_x^SL[n]), which is proportional to its integer-occupation one-electron error.
Delocalization error becomes the same thing as one-electron self-interaction at q = 1. A
standard GGA has curvature even when it is exact at q = 1.

**Corollary 4 (range separation).** With short-range screening erfc(omega r)/r, the same
argument gives F_1 = s^{-2} Psi(|grad n_s|/(n_s omega)): the semilocal SR one-orbital
exchange must be -n^4/|grad n|^2 times a function of the local decay constant over omega.
(Recorded only; the submission uses a global hybrid.)

## 1. Derivations

### 1a. The exact condition on the charge line
Take the ensemble (1-q)|H+> + q|H>. Its density is q n_H and its energy is -q/2. The KS ensemble
has one orbital phi_H with occupation q, so T_s = q T, V_ext = q V, and the Hartree energy is
J[q n] = q^2 J. Since T + V = -1/2 (E_xc[n_H] = -J), we get E_xc[q n_H] = -q/2 - qT - qV -
q^2 J = -q^2 J[n_H]. The ensemble correlation is zero, so a functional whose correlation is
gated to zero for fully polarized one-orbital densities must give E_x^{DFA}[q n_H] = -q^2 J for
all q. The same holds for a global hybrid. There E_x = a E_x^HF + (1-a) E_x^SL, and
E_x^HF[q n] = -q^2 J is already homogeneous, so the condition falls entirely on E_x^SL.

### 1b. Proof of Theorem 1
For a fully polarized one-orbital density, tau_s = tau_W,s, so z_s = 1 and alpha = 0. Every
(n, grad n, tau) meta-GGA then has e_x = e_x^{LDA,pol}(n) F_1(s), with
e_x^{LDA,pol} = -(3/4)(6/pi)^{1/3} n^{4/3}.
For n = q e^{-2r}/pi we have |grad n|/n = 2, so s(r) = u e^{2r/3} with u = (6 pi q)^{-1/3}.
Hence u >= s0 for q <= 1, and u = s0 exactly at the nucleus when q = 1. Substitute t = s(r)
(e^{-2r} = (u/t)^3, dr = (3/2) dt/t). Then

    E_x^SL[q n_H] ∝ q^{4/3} ∫_u^∞ (u/t)^4 ln^2(t/u) F_1(t) dt/t ∝ ∫_u^∞ t^{-5} ln^2(t/u) F_1(t) dt ≡ G(u),

and the condition -q^2 J ∝ u^{-6} becomes G(u) = C u^{-6} for all u >= s0. Differentiating:
-(u/2) G' = ∫ t^{-5} ln(t/u) F_1 = 3C u^{-6}; the next derivative gives ∫_u^∞ t^{-5} F_1 = 18 C u^{-6};
one more gives F_1(u) = 108 C u^{-2}. So F_1 = A s^{-2} on [s0, ∞), and the solution is unique.
Normalization: with F_1 = A/s^2, e_x = -(3/4)(6/pi)^{1/3}·4(6pi^2)^{2/3} A n^4/|grad n|^2 = -18 pi A n^4/|grad n|^2.
For H, ∫ n^4/|grad n|^2 = ∫ n^2/4 = 1/(32 pi), so E_x = -(9/16) A. Setting this to -5/16 gives
**A = 5/9**, and e_x = -10 pi n^4/|grad n|^2. Checked numerically on an STO-6G hydrogenic density
(homog_h2.py): -0.38996 against exact -0.38750 for that zeta = 1.24 density; the 0.6% is the
missing cusp of the Gaussian fit.

### 1c. Why s0 is the right boundary
A normalized one-electron exponential, n = (a^3/8pi) e^{-a d}, with local decay a = |grad n|/n at
distance d from its centre, has a d = ln(6 pi s^3). So s >= s0 is exactly the condition that the
local data (n, |grad n|) are consistent with a hydrogenic density carrying at most one electron.
Where s < s0 (density saddles such as the H2+ / H2 bond midpoint, maxima of p-type densities,
the 2s node region) no hydrogenic interpretation exists, and Theorem 1 constrains nothing. I use
a C1 polynomial continuation there:
F_1(s<s0) = (A/s0^2) p(x), x = s/s0, p(x) = 1 + (1-x^2)(1 + d(1-x^2)).
This has p(1) = 1, p'(1) = -2 (C1 for every d) and p(0) = 2 + d. One constant, d, is left. It is
fixed by a second exact one-electron norm, **the H2+ ion at R = 2.0 bohr** (one-electron, so
HF = exact and E_x = -J[n] exactly). The norm is computed by me from the HF density in Phase 2.
H2+ is the prototype of the two-centre, s < s0 region, and it is not in the dev set.


### 1c'. Continuation revised after a derivation check (h2p_model.py; LCAO H2+ density, no SCF)
The C1 family p = 1 + y + d y^2 cannot meet the H2+ norm with F >= 0. In an LCAO H2+ density
(zeta = 1.238, R = 2.0 bohr), exact -J = -0.30780. The part with s >= s0 already gives
-0.27687; the region s < s0 (3.1% of the electron) contributes -0.0656 at d = 0, and the norm
needs d = -3.38 < -2 (negative F). The reason: C1 matching at s0 (slope -2 in x) makes F rise
above 3.94 just below s0. So the continuation must drop immediately below s0. I therefore use
the C0 power family

    F_1(s < s0) = (A/s0^2) (s/s0)^m ,

which is continuous at s0. It has a kink in F there, but the energy density is continuous and
v_xc has no delta function. Hydrogenic densities never sample s < s0, so Theorem 1 is
untouched. H2+ model: E_x + J = -0.0164 / -0.0051 / +0.0018 / +0.0098 for m = 0 / 1 / 2 / 4, so the
norm gives m ≈ 1.75. m = 2 is the "log-reflection-symmetric" form
F_1 = (A/s0^2) exp(-2|ln(s/s0)|), which misses the norm by 0.6%. m will be fixed in Phase 2 from the
HF (= exact) H2+ density.

### 1d. Where the one-orbital branch is switched on
The theorem is exact for any one-orbital channel. In a many-electron system, F_1 has to be
blended with a many-orbital branch F_0 through a weight L_s that is 1 for a lone electron.
Physically motivated candidates, to be compared in Phase 2 and disclosed:
- W1 ("lone electron"): L_s = z_s^2 · max(0, zeta_s)^2, zeta_s = (n_s - n_s')/n. It is
  1 for any fully polarized one-orbital density (H, any q), 0 in every closed-shell region,
  and density-scaling invariant for the H atom.
- W2 ("unpaired one-orbital channel", the exchange twin of the pair-aware correlation gate):
  L_s = z_s^2 [1 - z_s'^2 (1 - zeta^2)]. It is also 1 for H and 0 for He/H2 pairs, but nonzero
  (<= 1/4) in closed-shell regions with z ~ 0.7.
Both give exactly L = 1 for the H atom at every q, so both satisfy Theorem 1.
The same indicator family already gates correlation. That is the point: both one-electron
conditions are switched on in the same regions.


### 1e. Theorem 2: the charge axis and the spin axis cannot both hold semilocally (motivates W1)
Statement: let E_xc be any meta-GGA in (n_s, grad n_s, tau_s), with no Laplacian. Take hydrogenic
densities n_Z = (Z^3/pi) e^{-2Zr}. Then the fractional-spin density (n_Z/2, n_Z/2) has, at radius r,
exactly the local data (n_s, |grad n_s|, grad n_a·grad n_b, tau_s = tau_W) of the two-electron
density (n_Z, n_Z) at radius r + δ, with δ = ln2/(2Z). The reason is
n_Z(r)/2 = n_Z(r+δ) and |grad n|/n = 2Z in both. So for any energy density e <= 0,

    |E[n_Z/2, n_Z/2]| = ∫_δ^∞ 4π (r'-δ)^2 |e| dr' <= |E[n_Z, n_Z]|.

Apply this to correlation. Suppose the exchange of paired one-orbital channels is made
charge-linear too (F_1 also where zeta = 0: candidate W2). Then exchange gives exactly
-J[n_Z]/2 for the fractional-spin density, and constancy (E_xc = -J = -(5/16)Z) needs
E_c = -(5/32)Z. A correlation that is negative pointwise and correct for the Z -> ∞ two-electron
ion (density 2n_Z, E_c -> -0.0467) is bounded by |E_c| <= 0.047. So the fractional-spin error of
H is at least 0.156 - 0.047 = 0.109 Eh (68 kcal/mol), and it grows linearly with Z.
Consequence: charge-linear exchange should be applied only to lone (unpaired) electrons.
That is the W1 weight, L_s = z_s^2 · max(0, zeta_s)^2, which is zero for paired
channels (zeta = 0), where the standard many-orbital exchange keeps the smaller semilocal
fractional-spin error. (The Laplacian breaks the shift map, since ∇²n carries an explicit 2a/r,
so Laplacian-level meta-GGAs escape Theorem 2.) This is a second new result. It is a
rigorous reason for one of the design choices, not a dev-set preference.

### 1f. Implementation decisions fixed in Phase 1
- One-orbital indicator: h(z) = z^3 (4 - 3z), with z capped at 1. This is ≈ z^2 in shape,
  equals 1 at z = 1 with zero slope, so the hard cap is C1 (prior art reports SCF trouble with a
  hard cap on z^2). It is used in both the exchange weight and the correlation gate.
- s-variables enter only through s^2 (no square roots): F_1 = A/s^2 for s^2 >= s0^2.
- Derivatives: complex-step differentiation of my own numpy PW92+PBE correlation and exchange,
  which is exact to machine precision.


### 1g. Bounds respected by F_1 (checked analytically)
- Negativity: F_1 >= 0 iff d >= -2 (p = 1 + y + d y^2 with y = 1 - x^2 in [0,1]).
- Integrated one-electron Lieb–Oxford bound, E_x[n_1] >= -1.092 ∫ n_1^{4/3} (local form F <= 1.174).
  For H, E_x = -0.3125 against a bound of -0.3146: the H atom nearly saturates the bound. For fractional
  q the ratio is 0.993 q^{2/3} < 1, so the bound holds along the whole charge line.
  The pointwise form is violated (Corollary 2). Densities with large flat regions (s < s0 everywhere)
  could in principle violate the integrated bound through the continuation. That is a stated
  caveat.
- Uniform coordinate scaling: F_1 depends on s only, so the functional is exact for all Z
  once it is exact for H.
- Density scaling: degree 2 for s >= s0 (Corollary 1). It is approximate in the continuation region.


### 1h. Numerical confirmation on the exact H density (charge_line.py; radial quadrature, no SCF)
Semilocal exchange E_x^SL[q n_H] (correlation gated to 0; exact = -5q^2/16):

| q | exact | F_1 = (5/9)/s^2 | PBE form | SCAN alpha=0 branch (H-normed) |
|---|---|---|---|---|
| 0.10 | -0.003125 | -0.003125 | -0.016883 | -0.014186 |
| 0.25 | -0.019531 | -0.019531 | -0.053103 | -0.048711 |
| 0.50 | -0.078125 | -0.078125 | -0.126997 | -0.123496 |
| 0.75 | -0.175781 | -0.175781 | -0.212191 | -0.212615 |
| 1.00 | -0.312500 | -0.312500 | -0.305939 | -0.312499 |

Deviation from linearity of the Hartree + semilocal exchange at q = 1/2, [E_Hx(1/2) - E_Hx(1)/2]:
F_1: 4e-15. PBE: -0.0522 Eh. SCAN's H-normed iso-orbital branch: -0.0454 Eh (-28 kcal/mol).
Being exact at q = 1 (the H norm) does nothing for the charge line. Theorem 1 is what is needed.

Theorem 2 check with libxc PBE correlation on hydrogenic densities. E_c[n_Z/2, n_Z/2] is
-0.0149 / -0.0162 / -0.0167 for Z = 1 / 2 / 5. That is below |E_c[n_Z, n_Z]| = 0.0415 / 0.0457 /
0.0481, as the inequality says, and far from the -(5/32)Z = -0.156 / -0.313 / -0.781 that
constancy would need on top of a charge-linear exchange.


### 1i. Corollary 4 made explicit (range-separated version, not used in the submission)
With erfc(omega r)/r screening, density scaling gives F(q^{-1/3}s, q^{-1/3}w) = q^{2/3}F(s,w) with
w = omega/k_s. So F = s^{-2} Psi(w/s), and w/s = 2 omega n_s/|grad n_s| = 2 omega/a (a is the local
decay constant). Exactness for every hydrogenic ion then fixes Psi completely:

    F_1^SR(s; omega) = (5/9) s^{-2} R(2 omega n_s/|grad n_s|),
    R(x) = 1 - (32/(5 pi)) ∫_0^∞ (1+u^2)^{-4} exp(-u^2/x^2) du,

which is the short-range fraction of the hydrogenic self-exchange (R -> 1 as omega -> 0, R -> 0 as
omega -> ∞). For H (a = 2, so the argument equals omega): R(0.1) = 0.823, R(0.3) = 0.533,
R(0.4) = 0.429, R(0.5) = 0.348, R(1) = 0.144. So a long-range-corrected hybrid also has a
unique, fully determined charge-linear lone-electron branch. The only free choice left is omega
itself. Recorded as an extension; there is no time to validate it here.


**Spin-scaling vs. fractional spin (consequence of Theorem 2, recorded honestly).** Exact spin
scaling, E_x[n_a,n_b] = (E_x[2n_a] + E_x[2n_b])/2, makes each channel independent. A charge-linear
lone channel then forces the same F_1 on paired one-orbital channels (He, H2, 1s cores), which is
candidate W2. By Theorem 2, W2 has a fractional-spin error of at least 68 kcal/mol for H, i.e. for
stretched H2. W1 avoids this by letting the exchange weight depend on zeta, which violates exact
spin scaling. Equivalently, W1 can be read as a spin-scaled charge-linear exchange plus a
nondynamic-correlation term (1-a) Σ_s h(z_s)(1 - zeta_s+^2)(E_x,s^{F0} - E_x,s^{F1}) in paired
one-orbital regions. That term restores the static-correlation mimicry of semilocal exchange,
and it is not negative pointwise. The theorem shows the trade-off is unavoidable for
Laplacian-free semilocal functionals. Neither choice is presented as satisfying both conditions.
Candidates: W0 (no F_1; ablation), W1, W2 (spin-scaling exact), and W3 = h(z_s)[1 - h(z_s')(1-zeta^2)]
(the exchange twin of the pair-aware correlation gate).

**Theorem 1 when correlation is not gated.** If a functional keeps self-correlation, the theorem
applies to the combined lone-electron enhancement (1-a)F_x + e_c/e_x^{LDA}. The conclusion about
existing functionals is unchanged, since none has that combined form either.


**W2 rejected at the derivation stage (h2_model_m.py).** F_1 with the H2+-normed continuation was
evaluated on model H2 one-orbital channel densities (Wang LCAO). Exact -J against F_1:
R = 1.4 bohr: -0.3211 vs -0.3465 (8% too negative; PBE form -0.3161). R = 2.0: -0.2827 vs -0.2948.
R = 3.0: -0.2344 vs -0.2114. In H2 this would overbind by about 0.75 × 2 × 0.025 Eh ≈ 24 kcal/mol.
So the hydrogenic-exact branch transfers badly to compact two-centre paired channels, as
Theorem 2 already predicted. W2 (and W3, which mixes F_1 into closed shells) is not taken to the
dev set as a serious candidate. W1 restricts F_1 to lone electrons. These are mostly atom-centred
and close to hydrogenic: H, the Li/Na valence s electron, the tails of open-shell p shells.

## 2. The full functional (planned)
E_xc = a E_x^HF + (1-a) Σ_s ∫ e_x^{LDA}(n_s) [(1-L_s) F_0(s_s) + L_s F_1(s_s)] + E_c
- a = 1/4: PBE0 adiabatic-connection argument (KNOWN, used, no credit).
- F_0 = PBE form with kappa = 0.804 (Lieb–Oxford, KNOWN). For mu I compare two known
  candidates: mu = beta pi^2/3 (PBE xc gradient cancellation) and mu = 10/81 (exact exchange
  GE2). Both are disclosed.
- E_c = PBE correlation with Stoll same-spin removal in one-orbital channels, made pair-aware:
  E_c^PBE[na,nb] - Σ_s z_s^2 [1 - z_s'^2(1-zeta^2)] E_c^PBE[n_s,0]. This is KNOWN: the
  B95/TPSS idea in the opus4 form from prior art, with no credit claimed. For fully polarized
  one-orbital densities it gives E_c = 0 at every q, which is what Theorem 1 needs.
- Exchange spin scaling: exact by construction (each channel uses its own n_s, s_s, z_s).
- New and credited: F_1 = A/s^2 with A = 5/9 and s0 = (6pi)^{-1/3} (Theorem 1), the d-continuation
  fixed by the H2+ norm, the corollaries, and the coupling of F_1 to the lone-electron weight.

Constants: A, s0 are analytic; d comes from the exact H2+ exchange; a, kappa, mu, beta are known.
Nothing is fitted to the dev set.

## 3. Planned verification (Phase 2)
1. check_derivs.py for both spin paths (complex-step derivatives inside eval_xc).
2. Fractional-charge SCF of H with occupation q = 0.25, 0.5, 0.75, 1 (aug-cc-pVQZ): E(q) linear
   with orbital relaxation, against PBE0 and the base without F_1.
3. Predictions not used for anything: He exchange from F_1 (exact -1.0258), Li, H2 at several R,
   and the fractional-spin error of H (expected to stay violated, see Section 4).
4. Dev set: candidates W1/W2 × mu, SVP screening, then one QZ run.
5. audit.py: expect LO "pointwise" flags in the z ≈ 1 samples (Corollary 2 says this is forced).

## 4. Conditions expected to be violated
- Fractional-spin constancy (flat-plane spin axis). For H with n_a = n_b = n/2, E_xc must be -J.
  Locally that density looks like a two-electron singlet (z_s = 1, zeta = 0), and no semilocal
  correlation can tell the two apart. Corollary-style argument: the s-ranges overlap.
- Pointwise one-electron LO bound F <= 1.174 (Corollary 2: violated on purpose).
- The s < s0 continuation is normed (H2+), not derived. Charge linearity there is approximate.
- No derivative discontinuity beyond 25% HF, no -1/r potential, no dispersion (unless VV10 is
  added with derived constants; VV10 notes below).

## 5. Side derivation: VV10 one-electron self-correlation (vv10_h.py, vv10_lvv.py)
On the exact H density, VV10 (b = 6, C = 0.0093) gives +4.54 mEh, made up of β·N = +4.85 and a
nonlocal part of -0.31 mEh. He-like (zeta = 1.6875): +9.2 mEh. The β part is a constant per
electron and cancels in every reaction energy, so the chemically relevant one-electron
self-correlation of VV10 is only the nonlocal -0.3 mEh/H. I also derived a closed-form
"local VV10" (∫ n[β - (3pi^2/32) n (omega0 kappa)^{-3/2}], which vanishes for the UEG). For
finite systems it is a poor approximation. Conclusion: VV10 self-interaction is chemically
negligible. No credit is claimed here. With exact H density and C = 0.0093, C6(H–H) = 5.94
(exact 6.499), so a C6-normed C would be about 0.0078.
