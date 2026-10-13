# Organiser review of opus5c (round 5, construction track)

## Independent verification (organiser, numerical)
- E_x = −10π ∫ n⁴/|∇n|² on the exact H density: −0.3125000 Ha (exact −5/16). ✔
- Fractional charge n → N n: E_x(N) = N²·(−5/16) at N = 0.5 and 0.25 to 1e-9, as piecewise linearity
  of E(N) requires (E_x must cancel J ∝ N²). ✔
- F = (5/9)/s² with spin-scaled exchange reproduces the same value. ✔
- Uniqueness argument: under n → N n, s → N^(−1/3) s and α is invariant, so n^(4/3) F(s, α=0)
  can scale as N² only if F ∝ s^(−2). Sound for Laplacian-free forms (q = ∇²n/(4k_F²n) would scale
  differently, hence the qualifier). ✔
- F(s0) = 3.93 at s0 = (6π)^(−1/3): incompatible with pointwise caps 1.174 and 1.804. ✔

## Classification
Theory claims 1–2 (uniqueness of charge-linear one-orbital exchange; incompatibility with
pointwise Lieb–Oxford-type caps) are verified and not in KNOWN_CONDITIONS.md as stated. Claims 3–4
(charge/spin incompatibility of ≥68 kcal/mol; variational instability of the charge-linear
potential) are plausible and documented but not independently verified here. Novelty relative to
the literature (Perdew–Ruzsinszky–Csonka 2007 on density scaling and many-electron SIE; Mori-Sánchez,
Cohen, Yang 2008 on delocalization) needs an expert check: the N^(−1/3) scaling of s is standard;
the closed form and the iff statement are not something the organiser has seen written down.

## Functional
Dev WTMAD 1.94; the cap at F ≤ 1.804 (needed for SCF stability) removes the exactness. Audit
DISTINCT. Candidate selection disclosed; c fixed by the Ne norm, not the dev set. Minor rule slip
(SVP runs started 4 min early, killed, rerun) disclosed.

Verdict: THEORY RESULT WORTH EXPERT REVIEW; functional not competitive.
