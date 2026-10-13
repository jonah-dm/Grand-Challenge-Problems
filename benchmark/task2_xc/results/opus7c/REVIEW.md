# Organiser review of opus7c (round 7, construction track)

Submission S1 = opus6c with per-channel τ_W,s/τ_s in the base exchange (restores exact spin scaling;
known condition, no novelty credit). Dev WTMAD 1.14 (AE 1.95, best of the benchmark; IP 7.89).
Delocalization fix intact: H₂⁺ at 6 bohr +0.09, f(½) = 0.2507, barrier 10.38. Audit DISTINCT.

Diagnosis (non-self-consistent decomposition on UHF densities, reproduces QZ to 1–2 kcal/mol):
IP(F) exchange error +21.3 of which the one-orbital branch +1.3; IP(HF) +20.6 / +0.2; IP(B) mostly
correlation (+5.8; Be-like static correlation in B⁺). Spatially: core–valence exchange inside
0.6 bohr grows by −21.7 kcal/mol on ionization; semilocal SR exchange changes by −1.0. Conclusion:
a single-ω RSH cannot fix this; needs exact exchange localized to the core–valence region
(local hybrid / local range separation, cf. Brütting–Bahmann–Kümmel 2022). Sound and useful.

Levers tested and rejected with scores (disclosed): SR-HF 1/4 fixes IPs (MAE ~4.7–5) but collapses
AEs at every derived ω; branch weight and correlation gate move IPs ≤ 1 kcal/mol.

New but not adopted: "bathtub" pointwise bound F ≤ (4π²/18)^(1/3) = 1.299 at a one-orbital density
maximum (attained by a uniform ball) → ω* = 0.5256; gauge caveat that opus6c's ω rule is a
regularization, not a theorem; hydrogenic floor on b (S8) halving the H₂⁺ 2-bohr error.

Selection: pre-registered rule picked S5; agent overrode after QZ on physical grounds (H atom,
H₂⁺, barrier) and recorded it. S8's rule was set before results and missed by 0.05 on the barrier.

Verdict: no novelty credit for the functional; diagnosis is the result; conduct exemplary.
