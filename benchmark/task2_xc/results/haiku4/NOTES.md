# haiku4 NOTES

Submission: `results/haiku4/submission.py` (= `cand_B_hnorm.py`), hybrid meta-GGA, HYB = 0.25, META = True.
Final dev-set WTMAD (def2-qzvp): **QZ_WTMAD_PLACEHOLDER** (QZ AE/IP MAE in `score.json`).
def2-svp WTMAD: 1.98 (AE MAE 4.33, IP MAE 7.34 kcal/mol). Baselines: PBE 2.95, r2SCAN 1.28, B3LYP 0.97, wB97M-V 0.56.
Audit (`audit_B.log`): UEG rel. dev 1.6e-15, LO violations 0, negativity 0, mix R2 0.9999, closest single PBE0
(R2 0.9997), verdict **NEAR-EXISTING**. check_derivs: OK (complex-step derivatives, max rel. err 1e-5).

Honest summary: this is a recombination of two prior routes (sonnet3's Stoll same-spin gate and haiku3's
H-atom-normed one-orbital exchange). Its dev score is within 0.01 of sonnet3 (1.00 at QZ) and it does not beat
the frontier. The audit and the dev score both say it is close to PBE0-type behaviour. Reviewer should weigh it as
a combination of existing constructions with one derived constant.

## What is new
1. **One-orbital exchange normed to the exact H atom.** Spin-scaled exchange, z = tau_W/tau (total, as in sonnet3).
   e_x = sum_s (1/2) e_x^unif(2 n_s) [1 + (1 - z^2) h_GGA(s) + z^2 h_orb(s)],
   h_GGA(s) = kappa (1 - 1/(1 + mu0 s^2/kappa)),  mu0 = beta pi^2/3 = 0.21951 (PBE value, z -> 0 branch),
   h_orb(s) = kappa (1 - 1/(1 + mu1 s^2/kappa)),  kappa = 0.804 (Lieb-Oxford: F_x <= 1.804 saturated at large s).
   mu1 = 0.2758318 is the unique root of E_x^SL[H atom] = -5/16 Ha (one constant, one exact datum, brentq at import;
   solved value confirmed to 1e-14). For H, z = 1 exactly, so this branch is the one that enters.
2. **Stoll same-spin gate on correlation (from sonnet3).** E_c = E_c^PBE[n_a,n_b] - sum_s z_s^2 E_c^PBE[n_s,0].
   A spin channel holding one orbital has no same-spin correlation. For H, E_c = 0 exactly.
   Consequence: the H-atom xc energy is exact by construction (exchange norm + zero correlation).
3. 25% exact exchange, semilocal exchange scaled by 0.75 (adiabatic-connection 1/4).

## Exact conditions
Satisfied: uniform-gas limit; spin scaling of exchange; Lieb-Oxford bound on F_x (audit: 0 violations); negativity;
exact H-atom exchange (by construction); exact H-atom correlation (zero); same-spin correlation zero in every
one-orbital channel; GE2 cancellation coefficient mu0 for the z -> 0 branch (PBE value).
Violated / not imposed:
- exchange is exact only in the integrated sense for H; not pointwise, and not for general one-electron densities;
- He exchange norm not imposed (opus3 found that no PBE-form iso-orbital exchange with two constants meets both H and He norms);
- no derivative discontinuity, no fractional-charge/spin linearity beyond what 25% EXX gives;
- no dispersion: noncovalent (NCI) classes will be weak;
- z^2 weight is an interpolation choice (not derived); total z used in both branches.

## Candidates compared (all disclosed; def2-svp dev-set scores; selection by dev score)
| id | one-orbital branch h_orb | mu1 | svp WTMAD | svp AE MAE | svp IP MAE | svp MCM | H-atom E_x (Ha) | status |
|---|---|---|---|---|---|---|---|---|
| A (`cand_A_tail.py`) | sonnet3 s^-1/2 tail: mu_PBE s^2/(1+C s^5/2) | - | 1.99 | 4.30 | 7.68 | 5.99 | -0.3098 | not selected |
| B (`cand_B_hnorm.py`) | LO-saturating PBE-form, mu1 solved from H norm | 0.27583 | **1.98** | 4.33 | 7.34 | 5.84 | -0.3125 | **selected** |

Logs: `score_svp_A_tail.log/json`, `score_svp_B_hnorm.log/json`. The difference is 0.01 in WTMAD and is within
noise of the choice; B was selected because it is the one that satisfies the H-atom exchange norm exactly, which is
the constraint the candidate is built on, and because its MCM is lower. Only B was run at def2-qzvp.
Sonnet3's own submitted form scores 1.99 at def2-svp in its notes; candidate A reproduces that (1.99), which is a
consistency check on this reimplementation.
No continuous parameter was fitted to the dev set. kappa, mu0, beta are fixed by LO / GE / PBE; mu1 by the H atom.

## Known weaknesses (svp per-item errors, kcal/mol)
Underbinding of X-H species: H2 -8.7, NH3 -8.3, H2O -11.3, HF -10.2, LiH -4.9. IP:HF -11.6, IP:B +8.0.
Second-row, barriers, NCI, PA not tested (no hidden-class data used).

## What I would try with more time
1. The X-H underbinding persists in every one-orbital-normed variant (haiku3, sonnet3, this one). The one-orbital
   branch needs a pointwise, not just integrated, constraint (for example the exact one-electron exchange hole
   radial moments), which would change the shape not the single constant.
2. Replace the z^2 weight by an exponent derived from a He/H2 model density (sonnet3 flagged this).
3. Add a range-separated exchange (RSH) for IPs and barriers, with the Stoll gate kept.
4. Add a non-local correlation term (VV10 or a derived dispersion kernel) for NCI.
