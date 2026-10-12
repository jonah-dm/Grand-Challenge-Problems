# Chemical Accuracy at Scale

## Background

Density functional theory (DFT) is the workhorse for computing the properties of molecules and materials: formation and reaction energies, structures, band gaps, phonons, magnetism and dynamics. It is exact in principle. In practice, two limits remain:

- **Accuracy.** All of the unknown physics is packed into the exchange–correlation (XC) functional. Today's functionals are typically off by 2–10 kcal/mol and fail qualitatively for band gaps, strong correlation and charge transfer. Chemical accuracy (≈1 kcal/mol) still requires methods like CCSD(T), which scale as O(N⁷).
- **Cost.** Standard Kohn–Sham DFT scales as O(N³), and hybrid functionals are worse. That caps routine calculations at roughly 10³ atoms, well below what interfaces, disordered materials and devices need.

These problems test whether AI systems can make real progress on either front. Each task is posed so that its answers can be **scored automatically and re-run as models improve**.

## Tasks

**1. Efficiency.**
Find approximations to DFT that scale linearly with system size without losing chemical accuracy. Implement this in [code TBD].
*Scored by:* measured scaling exponent, and error versus the standard O(N³) calculation (target ≤ 1 kcal/mol) on a fixed suite of insulators, semiconductors and metals.

**2. Accuracy.**
Design an XC functional, by parametrization or any other means, that outperforms [reference functional] on [properties] across [benchmark set].
*Scored by:* error on a hidden test set against CCSD(T), QMC or experimental reference data, with no property class allowed to get worse than a set tolerance.

**3. Theory.**
- (a) Find tighter bounds or higher-order expansions of known properties of the exact functional (e.g. the Lieb–Oxford bound, the gradient expansion, scaling relations).
- (b) Find new exact sum rules the exact XC functional should satisfy, and show that enforcing them gives better ground-state properties.

*Scored by:* proof checking (formal or expert review), numerical tests against exact or near-exact densities, and improvement on the Task 2 benchmark when the new constraint is enforced.

**4. Extensions.**
Extend the above results to generalizations of DFT, such as generalized Kohn–Sham (GKS), real-time and linear-response TDDFT, relativistic DFT and superconducting DFT.
*Scored by:* the analogous accuracy or cost metrics for that generalization (e.g. excitation energies, spin–orbit splittings, superconducting T_c).

## Benchmarking

- **Fixed harness:** standard inputs, reference data and scoring scripts, so results stay comparable over time.
- **Hidden test set:** refreshed periodically to avoid contamination.
- **Leaderboard per task:** tracks the best AI result against the best human/literature baseline (e.g. r²SCAN, ωB97M-V, existing linear-scaling codes).

## Open decisions
- Target code for Task 1
- Reference functional, properties and benchmark set for Task 2
- How Task 3 proofs are verified
