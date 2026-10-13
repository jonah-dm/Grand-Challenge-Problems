# opus4: pair-reassigned Stoll correlation (25% hybrid meta-GGA)

**Final dev-set WTMAD (def2-QZVP): 1.08** (AE MAE 2.39, IP MAE 3.89 kcal/mol; all 24 species converge).
`score.json` holds that run. Baselines: PBE 2.95, r2SCAN 1.28, PBE0 1.31, B3LYP 0.97, wB97M-V 0.56,
sonnet3 1.00. **On the dev set this does not beat sonnet3 or B3LYP.** It does beat them on the two
self-generated checks I could make for the hidden classes (H+H2 barrier, H2/LiH X-H bonds; see below).
Audit: DISTINCT (closest single PBE0 R^2 = 0.970, mix R^2 = 0.998). NEG 0, LO 0, UEG 1.6e-15.
check_derivs: OK for both spin paths (complex-step derivatives).

## What is new
I diagnosed the frontier functional (sonnet3) on data I generated myself. I used HF densities at
aug-cc-pVQZ and a level-6 grid (`corr_check.py`, `he_part.py`):

| E_c (Ha) | H | He | H2 (0.741 A) |
|---|---|---|---|
| exact | 0 | -0.0420 | -0.0409 |
| PBE | -0.0060 | -0.0420 | -0.0436 |
| sonnet3 (Stoll z_s^2 removal) | -0.0002 | **-0.0313** | **-0.0307** |
| **opus4 (this work)** | -0.0002 | -0.0417 | -0.0432 |

sonnet3 removes the Stoll "same-spin" piece E_c^PBE[n_s,0] wherever channel s holds one orbital.
In a two-electron singlet (He, H2, the 1s cores, X-H bond pairs) **both** channels hold one orbital,
so it removes 26% of the correlation energy. That correlation is real. In a singlet pair all
correlation is opposite-spin. PBE gets He right only because its opposite-spin piece (-0.0309) is
too small, and the Stoll same-spin piece (-0.0111) supplies the missing part. In a pair region the
"same-spin" label on that piece is wrong: it is opposite-spin correlation that PBE has
partitioned into the wrong slot. This 0.010 Ha deficit explains sonnet3's two largest AE errors
almost exactly: H2 -6.0 and LiH -5.8 kcal/mol (0.0096 Ha). It also lowers the H+H2 barrier.

**The new condition.** Same-spin correlation vanishes in a spin channel holding one orbital
(sonnet3). It is removed only to the extent that the region is **not** an opposite-spin pair.
Where both channels are one-orbital and both spins are present, that piece is reassigned to
opposite-spin correlation:

    E_c = E_c^PBE[na,nb] - sum_s z_s^2 [1 - z_s'^2 (1 - zeta^2)] E_c^PBE[n_s, 0]

z_s = tau_W,s/tau_s and s' is the other spin. The factor 1 - zeta^2 = 4 n_a n_b / n^2 is the
uncorrelated on-top pair-density indicator: it is 1 for a closed-shell pair and 0 where one spin is
absent. This is the "per-spin same-spin removal paired with an opposite-spin correlation built
for it" that prior rounds named but did not try. The opposite-spin correlation is
E_os = E_PBE - sum_s E_ss,s + sum_s z_s^2 z_s'^2 (1-zeta^2) E_ss,s. It has **no new constants**.
The He norm is met to 1% as a consequence of PBE's own He accuracy. It is not fitted.

Limits:
- One electron (H, any fully polarised one-orbital density): zeta = 1, so the full Stoll removal applies and E_c = 0.
- Two-electron singlet: E_c = E_c^PBE, which is within 1% of exact for He.
- Uniform gas: z_s = 0, so the result is PBE (the audit's UEG check passes).
- Open-shell atoms: in a channel with a single valence orbital (for example the beta 2s of C or N), the removal stays where the other channel is multi-orbital. The genuine opposite-spin correlation of open-shell atoms is never gated (this avoids the opus2/opus3 failure).

Exchange (25% HF and the z-blended PBE-form / s^-1/2-tail iso-orbital exchange, mu = beta pi^2/3,
kappa = 0.804) is **sonnet3's, unchanged**. My contribution is the correlation.

## Exact conditions
Satisfied:
- uniform-gas limit;
- Lieb-Oxford bound (pointwise, audit);
- e_xc <= 0;
- E_c ~ 0 for one-electron densities (2e-4 Ha residual from the smooth z saturation);
- no same-spin correlation in single-orbital channels, except in opposite-spin pair regions;
- singlet two-electron correlation is entirely opposite-spin and He-normed to 1%;
- PBE high-density limit and gradient-expansion cancellation (mu, beta);
- correct H-atom energy to 0.9 kcal/mol.

Violated or not imposed:
- **exact spin scaling of exchange.** The inherited exchange uses the total-density z in each spin channel. I tested the spin-scaled version (candidate E below). It is exact but scores worse, so I did not adopt it.
- the exchange gradient expansion coefficient (10/81);
- no derivative discontinuity or fractional-charge linearity beyond 25% HF;
- no long-range -1/r exchange potential;
- no dispersion (no VV10). For noncovalent interactions the functional relies on PBE0-like behaviour;
- H2 correlation is overestimated by 0.0023 Ha, as in PBE.

## Self-generated validation (not dev set, not used for any parameter)
H + H2 -> H2 + H barrier (linear TS, r = 0.930 A, aug-cc-pVTZ; exact 9.6 kcal/mol; `barrier_check.py`):
sonnet3 5.91, **opus4 7.13**, pairt (candidate D) 7.57. H2 De at aug-cc-pVTZ: sonnet3 103.5, opus4 111.1 (exact 109.5).
The barrier is one of the hidden classes. Here the new correlation helps by 1.2 kcal/mol, because
the TS has an extended two-electron bonding pair whose correlation sonnet3 removed.

## Candidates compared (all disclosed; no continuous parameter was tuned)
def2-SVP dev WTMAD (AE MAE / IP MAE):
| | correlation gate | SVP WTMAD | AE | IP |
|---|---|---|---|---|
| A | sonnet3 Stoll z_s^2 (reproduced) | 1.99 | 4.30 | 7.68 |
| B | pair, weight z_s'^2 without the (1-zeta^2) factor | Li atom did not converge (wandered to -7.33 Ha); run abandoned. z_s' is numerically ill-defined where n_s' -> 0, so I added the on-top factor | | |
| C1 | **pairz, z_s'^2 (1-zeta^2) (SUBMITTED)** | **1.82** | 3.82 | 7.86 |
| C2 | pairz, z_s' (1-zeta^2) (power 1) | 1.87 | 3.97 | 7.71 |
| D | pairt, z_tot^2 (1-zeta^2): total-density single-orbital indicator | 1.85 | 3.95 | 7.60 |
| E | C1 + spin-scaled z_s inside exchange (exact spin scaling) | 2.10 | 4.49 | 8.58 |
def2-QZVP: C1 **1.08** (AE 2.39, IP 3.89); sonnet3's A (published) 1.00 (AE 2.22, IP 3.56).
Only C1 was run at QZ. C1 was picked over D (both physically equivalent statements of "pair
region") on the SVP dev score, and over C2 because power 2 matches the z_s^2 of the removal term.
I also tried a long-range-corrected variant (RSH = (omega, 1, -0.75) with ITYH-attenuated
semilocal exchange; code kept behind O4_OMEGA, default off). I tried to fix omega non-empirically
by the IP theorem (`tune_h.py`). For H, -eps_HOMO is 0.328 / 0.447 / 0.488 Ha at omega = 0 / 0.4 / 1.0,
so it reaches 0.5 only as omega -> infinity (pure HF). For He, the tuned omega is about 0.85. No
finite, size-transferable omega follows from this, so I did not use RSH and did not score it.

## Dev-set errors (QZ, kcal/mol) vs sonnet3
H2 +1.6 (-6.0), LiH -0.4 (-5.8), CH4 +4.1 (-1.8), NH3 +1.4, H2O -3.0, HF -3.9, N2 -3.0 (-0.6),
CO -4.1 (-1.4), F2 -4.8, C2H2 +1.2, HCN -1.2, H2CO +0.9, OH 0.0, CH2 +3.9; IP B +6.2, F -0.5, HF -5.0.
The reassignment fixes the H-H and Li-H pairs. It overshoots C-H bonds (CH4, CH2) and slightly
underbinds multiple bonds. sonnet3's better dev score comes partly from cancellation: the 26%
correlation deficit in X-H pairs offsets the C-H overbinding of its exchange.

## What I would try with more time
- Fix the exchange side of the C-H / O-H imbalance with an exchange that obeys exact spin scaling (candidate E shows the inherited exchange is tuned to the total-z form).
- A He-and-H2-normed opposite-spin form, so that the H2 overcorrelation of PBE is not inherited.
- A system-size-aware range separation, since the IP-theorem tuning on small atoms gives omega that is too large.
- VV10 with b derived from the semilocal correlation's own gradient-damping length, for noncovalent interactions.
