# Round 5 — theory and construction tracks, 3 h, derivation-first, KNOWN_CONDITIONS gate

| Agent | Track | Central result | Dev WTMAD |
|---|---|---|---|
| opus5t | theory | C1 amplitude scaling (Lieb ensemble proof); T2 no-go (F = c s⁻² unique, unbounded); bounded construction F_q; H₂⁺ and barrier fixed; found harness vtau bug | 1.55 |
| opus5c | construction | Same uniqueness theorem via fractional H; closed form E_x = −10π∫n⁴/|∇n|²; LO-cap incompatibility; variational-instability lemma | 1.94 |
| sonnet5c | construction | Same condition from PPLB; obstruction theorem with numerical table; LP minimax floor (≥11.6%); Laplacian escape with held-out test | 1.53 |

See REVIEW.md in each folder. Organiser verified: closed form exact for H (−5/16), N² scaling,
uniqueness argument, Gaussian counterexample (2.5×), and the PySCF vtau bug.

## Assessment
First publication-adjacent output of the benchmark: a clean no-go theorem for Laplacian-free
semilocal exchange on one-electron delocalization, with a quantified floor and a constructive
escape. Three independent derivations in one round. Novelty vs. literature (Becke–Roussel 1989
Laplacian hole model; Perdew–Ruzsinszky–Csonka 2007 density scaling; MCY 2008) is unverified and
is the subject of the round-6 literature check. No functional beat B3LYP; the theorem explains why
the semilocal attempts of rounds 2–5 plateaued on self-interaction-sensitive quantities.
