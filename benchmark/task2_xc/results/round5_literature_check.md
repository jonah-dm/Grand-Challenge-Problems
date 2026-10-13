# Round 5: literature and novelty check of the one-orbital amplitude-scaling theorem

Reviewer: literature-check subagent. Date: 2026-10-13.
Inputs read: `results/round5_report.md`, `results/opus5t/NOTES.md`, `results/opus5c/NOTES.md` (theory sections),
`results/sonnet5c/NOTES.md` (sections 2–3b), `KNOWN_CONDITIONS.md`, and `results/sonnet5c/nonparam.py` (re-run).

**Access limitation.** WebSearch worked. WebFetch and curl to the publishers and to arxiv.org were blocked: DNS failure, and the proxy answered 403.
So every entry below rests on abstracts, search snippets, secondary citations, and my own knowledge of the papers. I mark
which is which. "NOT FOUND" means not found in abstracts or snippets. It does not mean the full texts were checked line by line.

## Summary verdicts

| claim | verdict | closest prior statements |
|---|---|---|
| (1) E_x[λn₁] = −λ²J[n₁] (hole integrates to −λ) | **KNOWN** | PPLB 1982 + Perdew–Zunger 1981 (both on KNOWN_CONDITIONS); Zhang & Yang 1998 (a "scaling relation necessary for … self-interaction free"); Mori-Sánchez, Cohen & Yang 2006; Perdew et al. PRA 2007 (exact-exchange hole sum rule for open systems of fluctuating N) |
| (2) Laplacian-free semilocal: amplitude ray ⇔ F₁ = c s⁻² (c = 5/9); shape-dependent ⇒ no (n,∇n,τ) functional is free of one-electron delocalization error | **PARTIALLY KNOWN.** The conclusion is known in a stronger form. The sharp "iff c s⁻²" uniqueness, the c = 5/9 constant and the closed form −10π∫n⁴/\|∇n\|² were NOT FOUND | Santra & Perdew 2019 ("semilocal functionals cannot be exact for all one-electron densities", because the Hartree energy is fully nonlocal); Perdew et al. PRA 2007 (semilocal holes normalized to −1 vs exact open-system holes); Ruzsinszky et al. 2006, 2007; MCY 2006/2008; Vydrov–Scuseria–Perdew 2007; Borgoo–Teale–Tozer 2012 (density-scaling homogeneity ≈ 4/3) |
| (3) LP minimax floor ≥ 11.6 % | **NOT FOUND, and NOT VALID as stated.** The number is a discretization artifact (see E1) | Hait & Head-Gordon 2018 (delocalization error ≈ quadratic in occupation; empirical) |
| (4) Laplacian invariant w = n∇²n/\|∇n\|² (= q/p) allows F = s⁻²G(w), which removes the obstruction | **PARTIALLY KNOWN, and overstated** (see E2) | Becke & Roussel 1989 (Laplacian + τ hole, exact for hydrogenic densities at integer occupation); Ramasamy, …, Sun, PRB 112, L161112 (2025): Laplacian meta-GGA that cuts one-electron SIE in H₂⁺; Perdew & Constantin 2007; Kaplan & Perdew 2022; Mejía-Rodríguez & Trickey 2017 |

## Per-reference findings

1. **Perdew & Zunger, PRB 23, 5048 (1981)**, and **Perdew, Parr, Levy & Balduz, PRL 49, 1691 (1982)**. Both are on KNOWN_CONDITIONS.
   Combined, E[λ] = λE[1] for 0 < λ ≤ 1, with T_s and V linear in λ and J[λn₁] = λ²J[n₁], gives (1) in two lines. sonnet5c's derivation is exactly this.
   opus5t's Lieb-functional route (convexity of T_W) is a clean, self-contained proof of the same energy statement. *Contains (1) implicitly.*
2. **Zhang & Yang, J. Chem. Phys. 109, 2604 (1998)**, "A challenge for density functionals: Self-interaction error increases for systems with a noninteger
   number of electrons". The abstract says the paper derives "a scaling relation necessary for an exchange-correlation functional to be
   self-interaction free". It discusses fractional-electron fragments (stretched radicals, transition states, CT complexes).
   From my recollection, not verified here because the full text was blocked: the relation is E_xc[δρ₁] = −J[δρ₁] = δ²E_xc[ρ₁] for 0 < δ ≤ 1, which is (1).
   *Contains (1) (high confidence). Contains no s⁻² uniqueness or semilocal no-go as far as the abstract shows.*
3. **Mori-Sánchez, Cohen & Yang, J. Chem. Phys. 125, 201102 (2006)**, doi:10.1063/1.2403848, "Many-electron self-interaction error". Defines SIE-freedom
   through linearity of E(N), and shows that one-electron-SIE-free functionals still fail. Also: Cohen, Mori-Sánchez & Yang, J. Chem. Phys. 126,
   191109 (2007), doi:10.1063/1.2741248. Mori-Sánchez, Cohen & Yang, PRL 100, 146401 (2008) (localization/delocalization error). Cohen, Mori-Sánchez &
   Yang, Science 321, 792 (2008). These papers carry the fractional-H E(N) curves and the convexity of semilocal functionals.
   *Related: they contain (1) as a condition and numerical violations. They contain no uniqueness theorem.*
4. **Ruzsinszky, Perdew, Csonka, Vydrov & Scuseria, J. Chem. Phys. 125, 194112 (2006)**, doi:10.1063/1.2387954 (spurious fractional charge on
   dissociated atoms), and **J. Chem. Phys. 126, 104102 (2007)**, doi:10.1063/1.2566637 ("one- and two-electron SIE-free, but not always many-electron
   SIE-free", H₂⁺, He₂⁺, LiH⁺, Ne₂⁺). These papers treat stretched H₂⁺ as the standard semilocal failure, and the many-electron-SIE concept as the
   generalization. *Related, empirical. No theorem about enhancement-factor form.*
5. **Perdew, Ruzsinszky, Csonka, Vydrov, Scuseria, Staroverov & Tao, PRA 76, 040501(R) (2007)**, doi:10.1103/PhysRevA.76.040501, "Exchange and
   correlation in open systems of fluctuating electron number". The abstract says semilocal functionals fail for separated open systems such as stretched A₂⁺. It
   develops "a theory of exact exchange" and "an exchange-hole sum rule" that explain the failure. The secondary source arXiv:0808.2523 (Perdew,
   Staroverov, Tao & Scuseria 2008) describes this as the exact hole of open subsystems integrating to a value > −1, whereas semilocal holes are normalized to −1.
   *Contains the hole-normalization half of (1) (−λ, not −1), and the qualitative mechanism behind (2). It contains no s⁻² uniqueness. Its remedy is a
   local hybrid, not a Laplacian.*
6. **Vydrov, Scuseria & Perdew, J. Chem. Phys. 126, 154109 (2007)** (tests on fractionally charged atoms). Semilocal functionals have the largest
   many-electron SIE and global hybrids reduce it only slightly. *This is the empirical counterpart of opus5t's T2(b), the global-hybrid statement.*
7. **Santra & Perdew, J. Chem. Phys. 150, 174106 (2019)**, doi:10.1063/1.5090534, and related Perdew-group SIC papers. These state that, because the Hartree energy is fully
   nonlocal, no (semi)local approximation can be exact for all one-electron densities. PZ-SIC achieves exactness for all collections of separated one-electron
   densities. *This contains the CONCLUSION of (2) in a stronger form: it covers every semilocal functional, including Laplacian-level ones. It is stated as an
   argument, not as an enhancement-factor theorem.* Perdew et al.'s Jacob's-ladder papers (e.g. JCP 123, 062201, 2005) put full one-electron
   SIE-freedom at the hyper-GGA rung for the same reason. That last attribution is from recollection.
8. **Becke & Roussel, PRA 39, 3761 (1989)**. The model hole uses n, ∇n, ∇²n and τ. It is exact for hydrogenic densities (at integer occupation) and
   normalized to −1 by construction. My own check: the BR89 parameter equation has RHS ∝ n^{5/3}/Q, and for a one-orbital density Q ∝ ∇²n. So under n → λn
   the hole parameter x changes as λ^{2/3}, and the hole still integrates to −1. **BR89 therefore violates (1) for λ < 1**. It is exact only at λ = 1. *This anticipates the
   ingredient in (4): the Laplacian lets local data identify the hydrogenic decay constant and the distance to the nucleus, the same mechanism as opus5t's
   y = 1/(1−w) = ζr. It does not anticipate the degree-2 amplitude-homogeneous form.*
9. **Ramasamy, Hou, Vega Bazantes, Irons, Wibowo-Teale, Lebeda & Sun, PRB 112, L161112 (2025)**, doi:10.1103/dkf7-jkb4 (arXiv:2410.08887, "How
   semilocal are semilocal DFAs? / Tackling the one-electron SIE within the semilocal framework"). This is a nonempirical exchange-only meta-GGA that uses the
   Laplacian and greatly reduces SIE for H₂⁺. Its binding curve is exact at equilibrium and better than PBE and SCAN over a range of R. **This is the closest prior work to (4).**
   I could not read the full text, so I cannot say whether it discusses amplitude (density) scaling, s⁻², or q/p. The round-6 authors must read it before claiming (4).
10. **Perdew & Constantin, PRB 75, 155109 (2007)** (Laplacian-level KE density and xc); **Kaplan & Perdew, Phys. Rev. Materials 6, 083803 (2022)**
    (Laplacian-level meta-GGA for metals); **Mejía-Rodríguez & Trickey, PRA 96, 052512 (2017)** and PRB 98, 115161 (2018) (deorbitalized SCAN-L:
    τ replaced by a Laplacian-level KEDF). These use the Laplacian to replace τ or to describe metals, not for fractional-charge exactness. The variable w
    is the standard ratio q/p of the reduced Laplacian and the squared reduced gradient. *They contain the ingredient, but nothing on (4)'s purpose.*
11. **Levy & Perdew, PRA 32, 2010 (1985)** (coordinate scaling) and the Levy-Perdew scaling literature deal with coordinate scaling, not amplitude scaling.
    The amplitude (density) scaling n → λn was studied as "effective homogeneity" by **Borgoo, Teale & Tozer, J. Chem. Phys. 136, 034101 (2012)**.
    They found the xc homogeneity averaged across the integer is ≈ 4/3. Tozer, J. Chem. Phys. 159, 244102 (2023) covers Fermi–Amaldi forms. This matches opus5t's
    "amplitude degree p" language. *Related: the concept is known, but not as a constraint that pins F.* Aside: the Fermi–Amaldi term −J/N is exact for all integer one-electron
    densities but gives −λJ on λn₁. That is an explicit nonlocal example that integer exactness does not imply (1).
12. **Kaplan, Levy & Perdew, Annu. Rev. Phys. Chem. 74, 193 (2023)**, doi:10.1146/annurev-physchem-062422-013259. This reviews about 20 exact constraints,
    including one-electron SIE-freedom and PWL, and the limits of the semilocal rungs. As far as the abstract and snippets show, it contains no amplitude-ray uniqueness theorem.
13. **Hait & Head-Gordon, JPCL 9, 6280 (2018)**, arXiv:1808.02770: delocalization error is essentially quadratic in fractional occupation, assessed for 47 functionals,
    with the H atom as a test. **Bryenton, Adeleke, Dale & Johnson, WIREs CMS 13, e1631 (2023)**, doi:10.1002/wcms.1631: a review of delocalization error as
    the greatest outstanding challenge, which includes one-electron SIE. **Lonsdale & Goerigk, PCCP 22, 15805 (2020)**, doi:10.1039/D0CP01275K (74 DFAs on
    hydrogenic mono- and dinuclear one-electron systems), and J. Chem. Phys. 158, 044102 (2023). **Kraisler & Kronik, PRL 110, 126403 (2013)**:
    ensemble generalization restores PWL for any approximate functional, which is a nonsemilocal escape. *These are all related or empirical. None contains the theorem.*
14. A meta-GGA review returned by search says it "seems unlikely that a meta-GGA will generally be able to meet the straight-line condition". This is a qualitative
    expectation, not a proof. The source was not identified precisely.

## Errors and overstatements found in the agents' derivations

**E1 (serious): the "LP minimax floor ≥ 11.6 %" is a discretization artifact, not a floor.** `sonnet5c/nonparam.py` restricts F to
piecewise-linear functions on a 70-node log-s grid on [0.01, 60]. A restricted class gives an *upper* bound on the true infimum, not a lower bound.
I re-ran the identical LP (same library, H plus H₂⁺ at R = 1–7, same amplitudes) with finer grids:

| grid nodes K | cap F ≤ 1.804: λ ∈ {1,.75} / {1,.75,.5} / {1,.75,.5,.25} | no cap (F ≤ 10⁴) |
|---|---|---|
| 70 (as submitted) | 5.3 / **11.6** / 20.9 % | 3.3 / 8.7 / 14.0 % |
| 300 | 3.9 / 9.8 / 18.5 % | 0 / 0 / 7.8 % |
| 1000 | 0.6 / 6.0 / 12.8 % | 0 / 0 / 0 % |
| 2000 | 0 / 3.0 / 8.7 % | — |
| 5000 | 0 / 0 / 2.6 % | — |

The "floor" goes to zero as the grid is refined, even with the Lieb–Oxford cap. The reason is that a finite library gives only 6 shapes × a few λ linear constraints, and
a bounded, oscillatory F(s) can satisfy any finite set of them. A genuine quantitative floor needs either infinitely many shapes or a regularity
or smoothness class for F, stated as part of the claim. The 11.6 % and 20.9 % numbers must not be reported as theorem-level bounds.
The qualitative no-go in (2) is unaffected, because it uses all shapes.
(Script: scratchpad `lpchk.py`, `lpchk2.py`. The scratchpad lies outside the repository, so these scripts were not saved there.)

**E2 (substantive): (4) does not "remove the obstruction"; it removes only the amplitude part.** For any semilocal functional of any finite
derivative order, including Laplacian and τ, exactness for all one-orbital densities is impossible even at λ = 1. Proof: take n₁ = ½(n_A + n_B), where n_A and n_B are smooth densities with disjoint
compact supports. The semilocal value E[½n_A] + E[½n_B] does not depend on the separation R. The exact value −J[n₁] contains the cross term
−½∫∫n_A n_B/|r−r'|, which does depend on R. This is the Santra–Perdew nonlocality argument made concrete. So:
- Laplacian-free: degree-2 amplitude homogeneity forces F = c s⁻², which is shape-wrong. So even the R → ∞ fragment limit (a sum of independent
  fractional fragments) cannot be right for both λ = 1 and λ < 1. That is the genuine content of (2).
- With w: F = s⁻²G(w) is homogeneous of degree 2 for every G. The fragment limit can then be right *if* G makes single-centre shapes exact. sonnet5c's fitted G
  reaches only 2–11 % on held-out shapes, and the 1/(2R)-type interfragment exchange of stretched H₂⁺ is still missed at any finite R.
  sonnet5c's statement "(A) holds for ALL λ by construction" is true only for the homogeneity. The λ = 1 error carries over unchanged to every λ.
  The round-5 report's sentence "the Laplacian invariant removes the obstruction" should read "removes the amplitude-scaling obstruction".

**E3 (minor, presentation): the T2 proof is informal, and a shorter route exists.** opus5t justifies "integrands equal" with a first-variation argument.
A cleaner chain: (1) for all λ and all shapes is equivalent to E_x^sl[m] = −J[m] for every density m with ∫m ≤ 1. Equate first functional derivatives on that set:
the semilocal side is local, while v_H is nonlocal, which is an immediate contradiction. The s⁻² uniqueness is then needed only for the weaker statement that is
genuinely informative: "degree-2 amplitude homogeneity on the z = 1 branch ⇔ F₁ = c s⁻²". Note also that for a semilocal functional, integer-N exactness on
separated multi-fragment one-electron densities already *implies* (1) at the fragment amplitudes (p and 1−p). So (1) is not an independent constraint on semilocal
forms. This is the old stretched-H₂⁺ argument (Perdew 1990 Adv. Quantum Chem. 21, 113 and Zhang–Yang 1998, from recollection).

**E4 (verified correct):**
- opus5c's hydrogenic-family uniqueness, Theorem 1: the triple differentiation G → F = 108C u⁻² checks, on s ≥ s₀ = (6π)^{−1/3}.
- c = A = 5/9 and e_x = −18πA n⁴/\|∇n\|² = −10π n⁴/\|∇n\|²: (6/π)^{1/3}(6π²)^{2/3} = 6π, and for H, ∫n⁴/\|∇n\|² = 1/(32π).
- The Gaussian ratio of 5: E = −5√(a/π) against exact −J = −√(a/π).
- The Lieb-functional proof of (1): convexity of T_W, and V_ee = 0 for N ≤ 1 ensembles.
- Spin-scaled s₀ at the H nucleus.
Opus5c's "Theorem 2" (shift map: fractional-spin H ↔ two-electron density at r + ln2/2Z) is a neat observation. It was not part of the check list and was not found in prior work.

**E5 (minor):** the round-5 report lists "Perdew–Ruzsinszky–Csonka 2007 density scaling". The 2007 Ruzsinszky/Perdew/Csonka papers concern fractional charge and
many-electron SIE. Density (amplitude) scaling homogeneity is Borgoo–Teale–Tozer 2012 and the earlier Liu/Parr/Chan–Handy work, which I did not verify here.

## Bottom line for novelty credit (against KNOWN_CONDITIONS)
- (1): no credit. PPLB and one-electron SIE-freedom are listed, and Zhang–Yang 1998, MCY 2006 and Perdew et al. 2007 state it, including the open-system hole sum rule.
- (2): the impossibility of a semilocal functional free of one-electron delocalization error is known in a stronger, all-semilocal form. Possible credit: the precise
  uniqueness result. Amplitude homogeneity on the iso-orbital branch forces F₁ = (5/9)s⁻², with closed form −10π∫n⁴/\|∇n\|², hydrogenic-family uniqueness
  for s ≥ (6π)^{−1/3}, and the Gaussian factor of 5. None of these was found in the literature searched. A full-text check of Zhang–Yang 1998, Perdew et al. PRA 2007 and
  Ramasamy et al. 2025 is still needed.
- (3): no credit. The number is a grid artifact (E1).
- (4): partial at most. A Laplacian route to one-electron SIE reduction exists (BR89; Ramasamy…Sun 2025). The specific degree-2-homogeneous family s⁻²G(q/p) was
  not found, but it must be described as removing only the amplitude obstruction (E2).

Sources (searched; abstracts or snippets only):
[Zhang–Yang 1998](https://pubs.aip.org/aip/jcp/article-abstract/109/7/2604/182940/A-challenge-for-density-functionals-Self) ·
[Perdew et al. PRA 2007](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.76.040501) ·
[arXiv:0808.2523](https://arxiv.org/abs/0808.2523) ·
[Santra–Perdew 2019](https://pubs.aip.org/aip/jcp/article/150/17/174106/198327/Perdew-Zunger-self-interaction-correction-How) ·
[Ramasamy et al. PRB 2025](https://journals.aps.org/prb/abstract/10.1103/dkf7-jkb4) / [arXiv:2410.08887](https://arxiv.org/abs/2410.08887) ·
[Becke–Roussel 1989](https://link.aps.org/doi/10.1103/PhysRevA.39.3761) ·
[Kaplan–Perdew 2022](https://arxiv.org/abs/2203.09403v2) ·
[Mejía-Rodríguez–Trickey 2017](https://link.aps.org/accepted/10.1103/PhysRevA.96.052512) ·
[Kaplan–Levy–Perdew 2023](https://arxiv.org/pdf/2207.03855) ·
[Borgoo–Teale–Tozer](https://www.lct.jussieu.fr/pagesperso/reinh/labo/seminaires/resume_borgoo.html) ·
[Hait–Head-Gordon 2018](https://arxiv.org/abs/1808.02770) ·
[Bryenton et al. 2023](https://research-repository.griffith.edu.au/items/8adac253-9bf3-4599-9308-c634554359ae) ·
[Lonsdale–Goerigk 2020](https://pubs.rsc.org/mk/content/articlelanding/2020/cp/d0cp01275k) ·
[Kraisler–Kronik 2013](https://arxiv.org/pdf/1211.5950) ·
[Perdew–Constantin 2007](https://arxiv.org/html/cond-mat/0612430)
