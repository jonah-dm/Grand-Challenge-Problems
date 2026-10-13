# Round 5 — theory and construction tracks, 3 h, derivation-first, KNOWN_CONDITIONS gate

| Agent | Track | Central result | Dev WTMAD |
|---|---|---|---|
| opus5t | theory | C1 amplitude scaling (Lieb ensemble proof); T2 no-go (F = c s⁻² unique, unbounded); bounded construction F_q; H₂⁺ and barrier fixed; found harness vtau bug | 1.55 |
| opus5c | construction | Same uniqueness theorem via fractional H; closed form E_x = −10π∫n⁴/|∇n|²; LO-cap incompatibility; variational-instability lemma | 1.94 |
| sonnet5c | construction | Same condition from PPLB; obstruction theorem with numerical table; LP minimax floor (≥11.6%); Laplacian escape with held-out test | 1.53 |

See REVIEW.md in each folder. Organiser verified: closed form exact for H (−5/16), N² scaling,
uniqueness argument, Gaussian counterexample (2.5×), and the PySCF vtau bug.

## Assessment, revised after the literature check (results/round5_literature_check.md)
- (1) λ² amplitude scaling: KNOWN (follows from PPLB 1982 + PZ 1981; Zhang & Yang 1998; MCY 2006).
- (2) No-go for Laplacian-free semilocal exchange: PARTIALLY KNOWN. Santra & Perdew 2019 state the
  stronger result (no semilocal functional is exact for all one-electron densities; Hartree is
  nonlocal). Creditable and not found in print: the sharp characterization F(s, z=1) = (5/9)s⁻²,
  the closed form E_x = −10π∫n⁴/|∇n|², the hydrogenic-family uniqueness, and the Gaussian
  counterexample.
- (3) LP floor ≥ 11.6%: INVALID. Grid artifact; finer grids give 6% (1000 pts) and 0% (5000 pts).
  The qualitative theorem stands; the number does not.
- (4) Laplacian escape: PARTIALLY KNOWN (Becke–Roussel 1989; Ramasamy et al., PRB 112, L161112
  (2025), a Laplacian meta-GGA that greatly reduces H₂⁺ SIE) and OVERSTATED: the Laplacian fixes
  only the amplitude scaling; the nonlocality obstruction at λ = 1 (separated fragments) remains.
- Derivations checked correct: ensemble proof of (1), uniqueness proof, c = 5/9, −10π prefactor,
  Gaussian factor, s₀. Misattribution corrected: density scaling is Borgoo–Teale–Tozer, JCP 136,
  034101 (2012).
Net: a correct, well-verified re-derivation with a few sharp new corollaries; not a publishable
theorem on its own. The benchmark's reviewer step (literature check) was necessary: without it the
round would have been recorded as a new result. No functional beat B3LYP.
