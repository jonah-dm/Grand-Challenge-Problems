# Round 3 — novelty rules + prior art + hidden-class disclosure

Rules as round 2, plus: round-2 NOTES given as prior art; hidden-set-B classes named (second-row
AEs, barrier, noncovalent, proton affinity); NEAR-EXISTING goes to review; disclosed form
selection allowed; check_derivs.py required; 90 min.

| Model | Submission | Audit | Dev WTMAD | AE MAE | IP MAE |
|---|---|---|---|---|---|
| haiku3 | meta-GGA: F_x = F0(s) + z(F1(s)−F0(s)), μ₁ solved from exact H-atom exchange; PBE correlation | NEAR-EXISTING (PBE, R²=0.9996) | 1.94 | 4.52 | 5.42 |
| sonnet3 | 25% hybrid meta-GGA; exchange blend from sonnet2; correlation E_c^PBE[n_a,n_b] − Σ_s z_s² E_c^PBE[n_s,0] (Stoll same-spin removal in one-orbital channels) | DISTINCT (PBE0, R²=0.97) | **1.00** | **2.22** | 3.56 |
| opus3 | "OEC-h": 25% hybrid; μ(z) exchange with μ₁ from H norm, μ₀ from xc gradient cancellation; correlation × (1−ζ⁴z²); H atom exact by construction | NEAR-EXISTING (PBE0, R²=0.9997) | 1.59 | 3.63 | 4.87 |

Baselines: PBE 2.95, r²SCAN 1.28, PBE0 1.31, B3LYP 0.97, ωB97M-V 0.56. Round 2 best: sonnet2 1.35.

## What happened
- **Cumulative progress.** Sonnet built on its own round-2 idea and went 1.35 → 1.00, now beating
  r²SCAN and within 0.03 of B3LYP, with no fitted parameter. Prior art as input works.
- **Convergent physics.** All three independently landed on the same exact condition, one-electron
  self-interaction/self-correlation via z = τ_W/τ, and the two strongest both chose 25% exact
  exchange from the adiabatic-connection argument. Opus made the sharpest formal point: in a hybrid
  the one-electron condition constrains the *sum* (1−a)E_x^SL + E_c, so fixing exchange or
  correlation alone breaks PBE0's error cancellation. Its functional is exact for the H atom.
- **Sonnet's Stoll partition is the best single idea so far**: removing only same-spin correlation
  in singly-occupied channels is an exact condition that also fixes the He/H₂ singlet case the
  round-2 gate got wrong.
- **Haiku is now competent** (derivative checker fixed its round-2 failure) but conservative.
- **Honest disclosure continued**: all three listed every candidate compared and its score.
- **Harness limits named by contestants**: no range separation (wanted for barriers/IPs), no
  nonlocal correlation (wanted for noncovalent). Both added to score.py for round 4 (RSH and
  VV10 tested against built-in ωB97X-V to machine precision).
- **Hidden-set-B metric**: the mean-|reference| scaling makes the 9.6 kcal/mol barrier and the
  5.4 kcal/mol dimer dominate; PBE's WTMAD there is 8.8. Ranking on set B uses MCM (plain mean of
  class MAEs) from now on.

## Rule / harness changes for round 4
1. RSH = (ω, α, β) and NLC = (b, C) supported in submissions.
2. Round-3 NOTES added to prior art; prompt names the frontier (Stoll same-spin removal; joint
   one-electron condition on the hybrid sum) and what each failed on.
3. MCM reported alongside WTMAD; hidden set B ranked on MCM.
