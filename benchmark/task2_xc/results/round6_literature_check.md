# Round 6: literature and novelty check of the opus6t and opus6c theory claims

Reviewer: literature-check subagent. Date: 2026-10-13.
Inputs read: `results/round5_literature_check.md`, `results/opus6t/NOTES.md`, `results/opus6t/REVIEW.md`, `results/opus6c/NOTES.md`,
`KNOWN_CONDITIONS.md`, `results/opus6t/work/lieb1d.py`, `v1.out`, `t4_hf.py`.

**Access limitation (same as round 5).** WebSearch worked. WebFetch and curl were blocked for every full-text host I tried: arxiv.org, ar5iv,
export.arxiv.org and the git.irsamc mirror all failed (DNS failure or a proxy 403). So every literature verdict below rests on abstracts, search snippets and my own
recollection of the papers, and I mark which is which. "NOT FOUND" means not found in abstracts or snippets. It does not mean the full texts were checked.
The numerical re-checks (section 3) were run by me in the scratchpad. They are not saved in the repository.

## 1. Summary verdicts

| claim | verdict | closest prior statements |
|---|---|---|
| A(i) C2: occupation square. Particle triangle E_xc = −N²J; hole-triangle sandwich −2f_af_bJ ≤ E_c ≤ −2(1−f_a)(1−f_b)J; tent; C2-HD | Particle triangle: **KNOWN** (opus6t correctly claims nothing for it). Hole-triangle bound: **PARTIALLY KNOWN**. It is a one-line consequence of the convexity of the Lieb functional, and it is not tight (E3 gives a strictly stronger bound). C2-HD: **PARTIALLY KNOWN**. It is the one-orbital instance of the Gould–Pittalis ensemble Hx (λ→0⁺) functional, rewritten in the 1-RDM exchange convention | PPLB 1982; Lieb 1983 (convexity of F_L); Cohen, Mori-Sánchez & Yang 2008 (constancy); MCY 2009 (flat plane); Gould & Dobson 2013 (LEXX, piecewise-linear ensemble Hx); Gould & Pittalis 2017; Kraisler & Kronik 2013; Yang, Mori-Sánchez & Cohen 2013 |
| A(ii) T3: on hydrogenic densities the flat-plane axis forces e_c = −(1−ζ²)(5π/(4(3π²)^{2/3})) n^{4/3}s^{−2} for s ≥ (3π)^{−1/3}; Gaussian factor 5 | Conclusion (no semilocal correlation satisfies the constancy condition for all shapes): **KNOWN** in substance. The explicit unique form: **NOT FOUND** in the literature, but it is a direct corollary of round-5 internal prior art. It equals (1−ζ²)/2 times opus5c's hydrogenic-family exchange uniqueness, and the factor 5 is round 5's | MCY 2009 (smooth functionals cannot produce the flat plane); Santra & Perdew 2019 (no semilocal functional is exact for all one-electron densities); round-5 opus5c Theorem 1 and the T2 form |
| A(iii) T4: e_c ≤ 0 semilocal correlation + flat-plane axis on hydrogenic ions ⇒ E_c ≤ −0.523Z for every hydrogenic singlet pair, He ≤ −0.876 Ha | **NOT FOUND** as a quantitative statement. The qualitative incompatibility between fractional-spin exactness and correct pair/closed-shell correlation in local forms is **KNOWN** (static-correlation-error literature, B05/B13, PSTS 2008, FSLOSC). Derivation verified (E1). Scope caveats in E2 | Cohen, Mori-Sánchez & Yang JCP 2008 and Science 2008; MCY 2009; Perdew, Staroverov, Tao & Scuseria 2008; Su, Li & Yang 2018; Bajaj, Janet & Kulik 2017; round-5 opus5c "shift map" (Theorem 2), which is the mechanism T4 uses |
| A(iv) exact 1D soft-Coulomb Lieb functional over the square, kink at N = 1 | The kink is **KNOWN**. It is exactly the PPLB derivative discontinuity I − A of the one-electron potential whose ground density is n₁ (E4), and I reproduce the reported slope from that alone. A numerical Lieb-functional map of the full one-orbital spin square: **NOT FOUND**. Closest are exact-functional landscapes for the two-site Hubbard model and 1D exact-functional studies | Cohen & Mori-Sánchez PRA 2016 (two-site Hubbard landscape); Yang, Mori-Sánchez & Cohen 2013; 1D soft-Coulomb exact-KS literature (Helbig et al.; Wagner et al.; Hodgson/Godby, from recollection) |
| B opus6c: bounded, Laplacian-free SR semilocal exchange branch, exactly amplitude-quadratic in one-orbital channels, in an LC hybrid with 100% LR-HF | **PARTIALLY KNOWN.** Every ingredient is known: LC/100% LR-HF, ITYH, the large-μ on-top expansion of SR exchange, the hydrogen SR-exchange expansion, using H to fix ω, and exact one-orbital exchange via local hybrids or local RSH. NOT FOUND: the specific branch e = n_s² g(b), with g fixed by the hydrogenic family, and the argument that C1 plus pointwise finiteness forces α_LR = 1. "Bounded" is overstated, and the cap breaks the exact λ² property on compact or cuspless shapes (E6) | Iikura et al. 2001; Vydrov & Scuseria 2006; Henderson, Janesko & Scuseria 2008; Vydrov, Scuseria & Perdew 2007; Cohen, Mori-Sánchez & Yang 2007; Toulouse, Colonna & Savin 2004/2005; Gill, Adamson & Pople 1996; Modrzejewski et al. 2016; Mussard & Toulouse 2017; Jaramillo, Scuseria & Ernzerhof 2003; Brütting, Bahmann & Kümmel 2022 |

## 2. Per-reference findings

### A. Fractional spin, flat plane, ensembles
1. **Perdew, Parr, Levy & Balduz, PRL 49, 1691 (1982)** and **Lieb, IJQC 24, 243 (1983)**, doi:10.1002/qua.560240302. The Lieb functional F_L is convex in n,
   and T_W is convex. These two facts give C2(i) (the ensemble of vacuum plus one-electron states has W = 0) and the C2(ii) upper bound (convex combination of the corners). *Contains the tools; the bound is not stated in this form (recollection).*
2. **Cohen, Mori-Sánchez & Yang, J. Chem. Phys. 129, 121104 (2008)**, doi:10.1063/1.2987202 ("Fractional spins and static correlation error"). This paper proves the constancy condition
   for fractional spins: E(N=1, any spin split) = E(H). All approximate functionals violate it, and the violation is the static-correlation error of stretched bonds. **Cohen, Mori-Sánchez & Yang,
   Science 321, 792 (2008)**, doi:10.1126/science.1158722 ("Insights into current limitations of DFT"): delocalization error and static-correlation error as the two faces of the
   fractional-charge and fractional-spin conditions. *Contains C2(i) at N = 1 as an energy statement (abstract). Contains no quantitative semilocal no-go.*
3. **Mori-Sánchez, Cohen & Yang, PRL 102, 066403 (2009)**, doi:10.1103/PhysRevLett.102.066403. The flat plane over (N_α, N_β), with a line of discontinuity at the integer. Its conclusion is that "explicitly discontinuous
   functionals ... beyond these currently used smooth approximations" are needed. *This is the qualitative form of opus6t's §3a remark (semilocal E_c is C¹ across N = 1, the exact one is not) and of the spirit
   of T3/T4. The flat plane over 0 ≤ N ≤ 1 at the exact density N·n_Z is what supplies the N < 1, ζ = 0 data that T3 and T4 need.*
4. **Gould & Dobson, J. Chem. Phys. 138, 014103 (2013)**, doi:10.1063/1.4773284 (arXiv:1206.6158). They point out that the split into Hartree, exchange and correlation is not unique for ensembles. They define LEXX,
   whose pair density and energy are piecewise linear in the occupations, with explicit expressions for frontier s and p shells. They state E_EEXX ≤ E_LEXX ≤ E_EXX. *In that ensemble convention
   the opus6t "static correlation" −2f_af_bJ of the particle triangle is part of ensemble Hx and not correlation. So C2(i)/C2-HD in the 1-RDM convention is a re-bookkeeping of
   known ensemble-Hx content (abstract plus recollection; full text blocked).* opus6t already flags this possibility, and I confirm it at abstract level.
5. **Gould & Pittalis, PRL 119, 243001 (2017)**, doi:10.1103/PhysRevLett.119.243001. They define E_Hx as the right derivative at λ = 0⁺ of the ensemble universal functional, computable from
   block eigenvalues of a simple matrix. For a single orbital this derivative is the LP minimum of W over the four-state manifold. That is exactly opus6t's C2-HD (minimum p₂ = max(0, N−1), so E_Hx = 2(N−1)₊J).
   *C2-HD's leading-order term is this functional evaluated in the 1-RDM convention. PARTIALLY KNOWN.* Follow-up: Gould, Pittalis, Toulouse, Kraisler & Kronik, PCCP 21, 19805 (2019).
6. **Kraisler & Kronik, PRL 110, 126403 (2013)** (ensemble generalization; Hxc for N ∈ [0,1] carries no Hartree self-interaction). **Gidopoulos, Papaconstantinou & Gross, PRL 88, 033003 (2002)** (ensemble
   ghost interaction). Both are background for E_xc = −N²J on one-electron ensembles (recollection).
7. **Yang, Mori-Sánchez & Cohen, J. Chem. Phys. 139, 104114 (2013)**, doi:10.1063/1.4817183 (arXiv:1305.5194). A framework that extends any G₀-functional to fractional charges and spins.
   The abstract does not show a fixed-shape bound or a semilocal no-go. **Cohen & Mori-Sánchez, PRA 93, 042511 (2016)**, doi:10.1103/PhysRevA.93.042511 ("Landscape of an exact energy functional"):
   the exact functional of the asymmetric two-site Hubbard model as a function of two variables. Its gap comes from derivatives of the exact functional. *This is the closest prior "exact functional map".
   It is a lattice model, not the continuum one-orbital square.*
8. **Chai, J. Chem. Phys. 136, 154104 (2012)**, doi:10.1063/1.3703894, and **140, 18A521 (2014)**, doi:10.1063/1.4867532 (TAO-DFT). Fractional occupations come from a fictitious temperature, and static correlation is carried by
   an entropy term. *Related remedy. It contains no bound of the C2 type.*
9. **Su, Li & Yang, PNAS 115, 9678 (2018)**, doi:10.1073/pnas.1807095115 (FSLOSC). This is an explicitly nonlocal fractional-spin correction that restores the flat plane. BLYP smooths the flat-plane discontinuity into a
   continuous surface. *The qualitative counterpart of T3/T4: the fix had to be orbital-based (nonsemilocal).*
10. **Bajaj, Janet & Kulik, J. Chem. Phys. 147, 191101 (2017)**, doi:10.1063/1.5001885. Semilocal DFT shows convex fractional-charge errors and concave fractional-spin errors. DFT+U removes the first and worsens the second, as does
    HF exchange. jmDFT fits low-order corrections to the shape of the error. *This is an empirical trade-off. It contains no theorem.*
11. **Perdew, Staroverov, Tao & Scuseria, PRA 78, 052513 (2008)** (arXiv:0808.2523). Full exact exchange needs correlation with "balanced nonlocality". Open subsystems of fluctuating N are among the "abnormal" regions.
    *Qualitative statement that static correlation of fractional or open systems is nonlocal.* Becke's B05/B13 real-space models, from recollection, use the exact-exchange hole normalization to detect
    static correlation, for the same reason.
12. **Marut, Senjean, Fromager & Loos, Faraday Discuss. 224, 402 (2020)**, doi:10.1039/D0FD00059K, and the GOK-ensemble review arXiv:2109.04943: weight-dependent local xc for He/H₂ and the ghost-interaction curvature.
    *Different ensemble (excitations). Related only.*
13. **Mussard & Toulouse, Mol. Phys. 115, 161 (2017)** (arXiv:1607.03621). In RSH, the fractional-spin error of H is larger than in standard KS and only slightly smaller than in HF. Adding LR-MP2 makes it infinite. *Relevant to the
    "hybrid does not reduce the static part" remark in T3(c), and to opus6c's admitted fractional-spin gap.*
14. No paper found gives a **quantitative** incompatibility bound between fractional-spin exactness and singlet correlation for semilocal forms. Searches covered "semilocal" + "constancy condition" + "two-electron"/"singlet"/"over-correlation"/"bound",
    MCY 2009, PSTS 2008, FSLOSC 2018 and Bajaj 2017 (all abstract level).

### B. Range separation and short-range exchange
15. **Iikura, Tsuneda, Yanai & Hirao, J. Chem. Phys. 115, 3540 (2001)**, doi:10.1063/1.1383587 (LC scheme; ITYH attenuation of GGA exchange through k_F/√F). **Vydrov & Scuseria, J. Chem. Phys. 125, 234109 (2006)**,
    doi:10.1063/1.2409292 (LC-ωPBE). **Henderson, Janesko & Scuseria, J. Chem. Phys. 128, 194105 (2008)**, doi:10.1063/1.2921797 (HJS model hole reproducing the parent GGA, analytic SR). *Known constructions.
    HJS is normalized to −1 and reproduces PBE, so it is not one-electron exact (abstract).*
16. **Vydrov, Scuseria & Perdew, J. Chem. Phys. 126, 154109 (2007)**, doi:10.1063/1.2723119. LC-ωPBE is "nearly many-electron SIE-free in many cases". **Cohen, Mori-Sánchez & Yang, J. Chem. Phys. 126, 191109 (2007)**,
    doi:10.1063/1.2741248 (MCY2 / rCAM-B3LYP; 100% LR exact exchange used to minimize fractional-charge error). *Known: LC with α_LR = 1 cuts delocalization error. It is not exact.*
17. **Baer, Livshits & Salzner, Annu. Rev. Phys. Chem. 61, 85 (2010)**, doi:10.1146/annurev.physchem.012809.103321, and **Kronik, Stein, Refaely-Abramson & Baer, JCTC 8, 1515 (2012)**, doi:10.1021/ct2009363
    (optimal tuning of ω via the IP theorem, i.e. piecewise linearity). *A different route to fixing ω from a one-electron/PWL condition.*
18. **Gill, Adamson & Pople, Mol. Phys. 88, 1005 (1996)** (Coulomb-attenuated functionals). **Toulouse, Colonna & Savin, PRA 70, 062505 (2004)** and **J. Chem. Phys. 122, 014110 (2005)**: as μ → ∞ the SR xc energy becomes a local
    functional of the on-top pair density and its derivatives. A gradient expansion of SR xc is also given. *opus6c's (E2), −(π/2ω²)∫n_s²(1 − b/2), is the one-orbital specialization of this expansion. For one orbital,
    D_σ = 0, so the hole curvature is ∇²n/6, and integration by parts gives b. KNOWN in substance. opus6c itself credits these works.*
19. **Modrzejewski, Hapka, Chałasiński & Szczęśniak, JCTC 12, 3662 (2016)**, doi:10.1021/acs.jctc.6b00406 (LC-PBETPSS). They use the one-electron SIE of the H atom to choose ω (SR-PBE exchange exact for H at ω ≈ 0.33).
    They quote the exact SR exchange of H at large ω as −1/(16ω²) + 1/(32ω⁴) + …, and the snippet says LDA recovers the leading term. *This matches opus6c's R(b) series term for term
    (−(5/16)(b/5 − b²/10)). So "exact hydrogenic SR self-exchange vs ω" and "fix ω on the H atom" are PARTIALLY KNOWN. The full closed form R(b) and the LO-cap rule for ω (L4) were NOT FOUND.*
    Related (recollection, unverified): Brémond, Pérez-Jiménez, Sancho-García & Adamo, J. Chem. Phys. 150, 201102 (2019), "Range-separated hybrid density functionals made simple", fix ω from the H-atom total energy.
20. **Jaramillo, Scuseria & Ernzerhof, J. Chem. Phys. 118, 1068 (2003)**, doi:10.1063/1.1528936. Local hybrid with mixing τ_W/τ, which equals 1 in one-orbital regions. So exchange there is 100% exact, exactly λ²-homogeneous and exact for every
    shape (recollection, standard). **Brütting, Bahmann & Kümmel, J. Chem. Phys. 156, 104109 (2022)**, doi:10.1063/5.0082957: local range separation ω(r) with a one-electron SIE correction and correct
    density-scaling behaviour. A parameter-free successor (ωBT21/23) is described as one-electron-SIE-free (DFT2024 abstract, Benasque 2025 slides). **Schmidt, Kraisler, Makmal, Kronik & Kümmel, J. Chem. Phys. 140, 18A510 (2014)**
    (SIE-free local hybrid; recollection). *The physical goal of opus6c (exact amplitude-quadratic one-orbital exchange inside an RSH) is reached in prior work with exact-exchange energy densities. opus6c's
    contribution is to reach the amplitude part with a semilocal SR branch.*
21. **Janesko**: no paper found with an SR exchange hole that is exact in the one-electron limit. Related, from recollection: Janesko's "Rung 3.5" model-density-matrix functionals (J. Chem. Phys. 133, 104103, 2010) and Janesko, Scalmani & Frisch,
    J. Chem. Phys. 141, 144104 (2014) on delocalization.
22. **Amplitude (λ-) scaling of SR semilocal exchange**: no paper found that treats E_x^SR[λn₁] ∝ λ² as a design constraint. Density-scaling homogeneity is studied by Borgoo, Teale & Tozer (2012), and
    Brütting et al. (2022) impose a density-scaling constraint on ω(r). Neither is the amplitude-quadratic one-orbital condition (abstract level). The large-ω leading term −(π/2ω²)n_s² of SR-LSDA is amplitude-quadratic for
    every system. That observation is trivial but worth stating: any SR functional is λ²-exact at leading order in 1/ω.

## 3. Errors, overstatements and verifications (my re-derivations)

**E1 (verified correct): T3 and T4 algebra.** I re-derived the following:
- Hydrogenic curve: s k_F = Z, and the domain is s ≥ (3π)^{−1/3} = 0.4730 (n ≤ Z³/π ⇔ 3πs³ ≥ 1).
- Mellin form E(w) = (π/(2Z³))∫₀^w G(u) ln²(w/u) du/u, so D³E = (π/Z³)G(w).
- C2(i) gives E = −(1−ζ²)(5π²/32)w²/Z⁵. Hence G = −(1−ζ²)(5π/4)u²/Z², which equals −(1−ζ²)·5π n⁴/|∇n|². That is exactly (1−ζ²)/2 times round 5's −10πn⁴/|∇n|² law.
- The forced form gives −5Z/32 = −J/2 for the N = 1, ζ = 0 ion (checked numerically).
- T4 constant: ∫₀¹ t ln²(2/t) dt = ln²2/2 + ln2/2 + 1/4 = 0.83680, so −0.52300Z and −0.88256 Ha at Z = 27/16 (closed form and radial quadrature agree to 1e−10).
- The s at the nucleus is (3πN)^{−1/3}.
- C2(ii) algebra: −[1+(N−1)²]J + (f_a²+f_b²)J = −2(1−f_a)(1−f_b)J, and min{f_af_b, (1−f_a)(1−f_b)} = f_af_b ⇔ N ≤ 1. Both correct.

**E2 (scope of T3/T4; must be stated with the claim).**
(a) The hypothesis is *exactness of E_c at the exact hydrogenic fractional-spin densities N·n_Z, for all N ≤ 1 at ζ = 0 and for a continuum of Z*. The N = 1 constancy condition alone does not
determine G (one point per Z-curve). The N < 1 data come from the lower triangle of the MCY flat plane. That is legitimate, but T4 must say "flat plane on 0 ≤ N ≤ 1", not "constancy condition".
(b) It is a density-functional statement at fixed densities. A functional whose *self-consistent* energies satisfy the flat plane, with relaxed densities, is not covered. That is the usual meaning of "satisfies
the flat-plane condition" for an approximate functional. Round 5's T2 has the same caveat.
(c) The correlation is defined against the 1-RDM Fock exchange −(f_a²+f_b²)J. In the Gould–Dobson/Gould–Pittalis ensemble-Hx convention, the particle-triangle E_c is 0, so "the flat plane forces semilocal
*correlation* to ..." is convention-bound. The convention-free version holds for the semilocal remainder of E_xc. I checked that it survives: for a fully semilocal xc with e_xc ≤ 0 the same argument gives E_xc ≤ −1.046Z
for hydrogenic pairs (−1.77 Ha at Z = 27/16, against exact E_xc(He) ≈ −1.07). So T4's message is robust to the convention, with prefactor 2.
(d) It does not extend to Laplacian-level functionals. On the hydrogenic curve q/p depends on r, so ∇²n separates N. opus6t states this.
(e) T4 is quantitatively the round-5 opus5c "shift map" (the outer part of a hydrogenic pair density is an N = 1 density shifted by ln2/(2Z)) combined with T3. That lineage should be cited as internal prior art.

**E3 (C2(ii) is not tight, and a strictly stronger bound follows from the same argument).** The trial-ensemble proof is the convexity of F_L. Using the *exact* pair corner instead of the determinant gives
  E_c(f_a,f_b) ≤ −2(1−f_a)(1−f_b)J + (N−1)·E_c^exact[n₁,n₁]   (N ≥ 1),
because F[n₁,n₁] = 2T_W + 2J + E_c[n₁,n₁] and F_L is convex. This explains opus6t's "smooth gap that grows to the pair dynamic correlation": the gap is at least linear in N−1. Every opus6t V1 point satisfies it
(for example (3/4,3/4): bound −0.1525J, computed −0.1707J; (1,½): bound −0.0275J, computed −0.0426J). C2(ii) should be presented as the determinant-corner weakening of convexity, and only
C2-HD (its exactness at O(g)) as the content beyond Lieb 1983. That content is itself the one-orbital Gould–Pittalis E_Hx (reference 5).

**E4 (A(iv): the kink is the PPLB derivative discontinuity, quantitatively).** Along the fixed-shape ray (f,f), F'(1⁻) = T_W, and the subdifferential argument gives F'(1⁺) − F'(1⁻) = I − A of the one-electron potential v
with ground density n₁. Hence E_c'(1⁺) = −J + (I − A). For the soft-Coulomb H potential I computed I − A = 0.6101 Ha = 1.762J on the same 61-point grid (0.6094 Ha, 1.763J on 121 points;
1.762J on 161 points with L = 16). This predicts a right slope of +0.762J, and opus6t reports "about +0.75J" (a secant to N = 1.25). The numbers therefore agree, and the kink is known physics (PPLB), not a new
feature. I also note that **every V1 point is reported conv=False (L-BFGS hit 1500 iterations)**. This is benign: Ω_T ≤ Ω and any v gives a dual lower bound, so all computed F, and hence E_c, are rigorous *lower* bounds.
Consequences:
- The particle-triangle "equality to 0.1%" is meaningful.
- "Inside the tent" on the hole triangle is automatic for a lower bound and does not verify the upper bound. The upper bound is proven analytically anyway.
- The reported gaps (0.032–0.052J) are upper estimates of the true gap.

No grid artifact of the round-5 kind was found. The particle-triangle identity and the kink are grid-independent theorems, and I − A/J is stable to 0.1% under grid refinement.

**E5 (T4' numbers: basis/grid sensitivity, not an artifact).** I re-ran the T4' integral for He on HF densities:
- cc-pVTZ, grid level 3: −0.848
- aug-cc-pVQZ, grid level 3: −0.862
- aug-cc-pVQZ, grid level 6: −0.876 (opus6t's number)
- aug-cc-pV5Z, grid level 6: −0.887 Ha

The drift comes from the cuspless Gaussian-basis nucleus, which moves the s = s* boundary. Quote He ≤ −0.85 to −0.89 Ha (hydrogenic Z = 27/16: −0.883). The 20× over-correlation conclusion is unaffected.

**E6 (opus6c: "bounded" and "exact for every shape" are overstated).**
(a) L2 is correct as stated: R(0) = 0, R is increasing and concave, and Q = R/b decreases from 1/5. But F^{SR,1} = (5/9)(k_s/ω)²Q(b) ≤ (k_s/ω)²/9, which grows like n^{2/3}. It is bounded in s at fixed n,
*not* uniformly. In high-density regions (Z > 1 nuclei, cuspless Gaussian centres) it exceeds any LO-type cap. So the branch is not "bounded" in the sense the s⁻² law fails (a uniform F ≤ const). It is
finite pointwise, and a cap is still required.
(b) L1 ("exact C1 for every shape") holds for the uncapped branch only. The cap F_c = F(1+(F/F_cap)⁸)^{−1/8} is not amplitude-invariant. My check uses the branch alone on a normalized Gaussian, with 4E[n/2]/E[n] (exact 1):

| ω, cap | Gaussian a | uncapped | capped |
|---|---|---|---|
| 0.3687, 1.804 | 1 | 1.000 | 1.204 |
| 0.3687, 1.804 | 1/4 | 1.000 | 1.010 |
| 0.5767, 1.174 | 1 | 1.000 | 1.176 |
| 0.5767, 1.174 | 1/4 | 1.000 | 1.0005 |

opus6c's own table shows the same thing for the full functional (fG ratio 1.074 for K1b, 1.116 for K1, diluted by the exact-exchange parts). So C1 is exact only where the cap is inactive. NOTES §7 lists "C1 amplitude-quadratic for every
one-orbital shape (L1)" as satisfied. That entry should be qualified.
(c) The λ = 1 shape error of the branch is large for compact cuspless shapes. For a Gaussian with a = 1 at ω = 0.5767, the uncapped SR branch is 1.45× the exact SR self-exchange (1.81× at ω = 0.3687). This is consistent with opus6c's admitted
shape error, but it is large.
(d) L4/K1b: moving the cap from 1.804 to 1.174 came after a phase-2 SCF failure. opus6c discloses this. The 1.174 bound is a global two-electron bound used pointwise, as SCAN does, and pointwise validity is not proven.
I verified the numbers ω* = 0.368708 (cap 1.804) and 0.576744 (cap 1.174) by root-finding. The constant (5/9)(6π)^{2/3} is 3.935, not 3.937, which is immaterial.
(e) Verified correct:
- The moments M₀ = π/ω² and M₁ = 3π/(4ω⁴), and the (E2) coefficient −1/2 in b.
- R(b) against direct Fourier quadrature of the erfc-screened hydrogenic self-exchange at (Z, ω) = (1, 0.3), (1, 0.5768) and (2, 0.4). It agrees to 1e−14.
- The Gaussian SR self-exchange −√(a/π)(1 − ω/√(a+ω²)).
- The small-b series, which matches Modrzejewski et al.'s −1/(16ω²) + 1/(32ω⁴).
- The on-top limit −(π/2ω²)n_s², which equals the SR-LSDA leading term.

**E7 (minor).** opus6t's REVIEW.md attributes the claim "semilocal forms cannot satisfy the flat plane" to MCY 2009. At abstract level MCY 2009 says "smooth approximations" in general, which is broader than semilocal. So that attribution is fine.
opus6t's V2 "first demonstration that all GGA/meta-GGA correlations violate the hole-triangle bound at (3/4,3/4)" is new as stated (not found), but it is a fixed-density check of a weak bound. With E3's stronger
bound the violation is larger.

## 4. Bottom line for novelty credit (against KNOWN_CONDITIONS)
- **A(i)**: particle triangle, no credit (flat plane/constancy, PPLB, Lieb). Hole-triangle sandwich: little credit. It is a direct consequence of Lieb convexity, and E3 gives a tighter version. C2-HD: partial credit at most, as an explicit
  1-RDM-convention restatement of the one-orbital Gould–Pittalis E_Hx and its consequence that semilocal correlations with a finite high-density limit fail at O(g) on the open square.
- **A(ii) T3**: little credit. The conclusion is known (MCY 2009; Santra–Perdew 2019), and the unique form is a corollary of round-5 opus5c/T2 with a (1−ζ²)/2 factor. It is correct.
- **A(iii) T4/T4'**: the most defensible new item. A quantitative, proven over-correlation bound (−0.523Z for hydrogenic pairs; He ≤ about −0.85 to −0.89 Ha on HF densities) was not found in the literature. Credit it with the scope
  caveats in E2: fixed exact densities, flat plane on N ≤ 1 not just constancy, 1-RDM convention (it survives in xc form), and round-5 opus5c shift-map lineage. The qualitative incompatibility is known.
- **A(iv)**: the kink is the PPLB derivative discontinuity and reproduces I − A quantitatively (E4), so no credit for the kink. The 1D Lieb-functional map of the one-orbital spin square is a useful verification, not a new result.
  Its hole-triangle values are lower bounds (unconverged dual).
- **B (opus6c)**: partial credit. New as far as found: the semilocal SR branch e = n_s²·5Q(b) fixed by the hydrogenic family, exactly λ²-homogeneous where uncapped; pointwise finiteness at density critical points;
  and the C1 + finiteness argument for α_LR = 1. Known: LC/100% LR-HF, ITYH, the large-μ on-top expansion (E2 in their notes), the hydrogen SR-exchange expansion, fixing ω from H, and exact one-orbital exchange via local hybrids or local RSH.
  The claims "bounded" and "C1 exact for every shape" must be qualified (E6).

Sources (searched; abstracts or snippets only):
[Gould–Dobson 2013 (arXiv:1206.6158)](https://www.arxiv.org/abs/1206.6158) ·
[Gould–Pittalis 2017 (arXiv:1708.03966)](https://arxiv.org/abs/1708.03966v4) ·
[MCY 2009 (arXiv:0809.5108)](https://arxiv.org/abs/0809.5108) ·
[CMY 2008 fractional spins (arXiv:0805.1724)](https://arxiv.org/abs/0805.1724) ·
[Yang–Mori-Sánchez–Cohen 2013 (arXiv:1305.5194)](https://www.arxiv.org/abs/1305.5194) ·
[Cohen–Mori-Sánchez 2016 (arXiv:1506.02230)](https://arxiv.org/pdf/1506.02230) ·
[Bajaj–Janet–Kulik 2017](https://hjkgrp.mit.edu/publication/bajaj-communication-2017/) ·
[Su–Li–Yang 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6166809) ·
[Chai TAO-DFT](https://manual.q-chem.com/6.3/sect_tao-dft.html) ·
[Marut et al. 2020](https://pubs.rsc.org/es-co/content/articlelanding/2020/fd/d0fd00059k) ·
[Mussard–Toulouse 2017](https://hal.sorbonne-universite.fr/hal-01344404) ·
[Toulouse–Colonna–Savin 2004](https://arxiv.org/abs/physics/0410062) ·
[Modrzejewski et al. 2016 (arXiv:1608.00035)](https://arxiv.org/pdf/1608.00035) ·
[Brütting–Bahmann–Kümmel 2022](https://epub.uni-bayreuth.de/id/eprint/8471/1/v2_final.pdf) ·
[Vydrov–Scuseria–Perdew 2007](https://arxiv.org/abs/cond-mat/0702283) ·
[PSTS 2008 (arXiv:0808.2523)](https://arxiv.org/pdf/0808.2523) ·
[Range-separated hybrids made simple](https://rua.ua.es/dspace/handle/10045/92748)
