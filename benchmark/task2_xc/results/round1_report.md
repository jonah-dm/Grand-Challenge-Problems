# Round 1 (pilot) — loose rules

Rules: "design an XC functional that beats the baselines"; libxc building blocks allowed; 45 min,
1 thread, hidden set announced but not described. Harness bug present during the round: custom
functionals with HYB > 0 silently lost their exact exchange (fixed in score.py afterwards).

| Model | Submission | Audit verdict | Dev WTMAD (QZ) | Notes |
|---|---|---|---|---|
| Haiku | B88 exchange + LYP correlation (= BLYP) | WRAPPER (BLYP, R²=1.000) | 1.68 | Tried PBE0-type hybrids; got nonsense scores because of the HYB bug; fell back to BLYP. |
| Sonnet | 0.4 MN15 + 0.4 M06-2X + 0.2 PW6B95 | MIX (recovered weights exactly) | n/a (no valid run; rescored by organisers) | Screened 7 libxc hybrids on the dev set at TZ, picked the best blend. Hit the HYB bug, could not diagnose it. |
| Opus | r²SCAN unchanged | WRAPPER (R2SCAN, R²=1.000) | 1.28 (= baseline) | Screened hybrids at TZ, planned 0.5 M06-2X + 0.5 MN15; found and correctly diagnosed the HYB bug; declined to monkeypatch the scorer; fell back to r²SCAN. |

Baselines: PBE 2.95, r²SCAN 1.28, B3LYP 0.97, ωB97M-V 0.56.

## What happened
- **All three did functional shopping.** Every agent's first move was to screen existing libxc
  functionals against the dev set and submit the best one or a blend. None attempted a new form,
  a new constraint, or any derivation. Under the loose rules this is the rational strategy.
- **None beat the best baseline.** ωB97M-V (0.56) is already a strong functional on these small
  molecules; blends of Minnesota hybrids were estimated at ~0.67–0.72 at best.
- **The harness bug dominated the round.** All three hit it. Only Opus diagnosed it correctly
  (PySCF decides whether to build exact exchange from the name in `mf.xc`, not from the `hyb`
  passed to `define_xc_`). Haiku noticed "hyb has no effect" but did not find why. Sonnet
  misattributed it to array layout.
- **Audit tool works** once its library is broad enough: it identified BLYP and r²SCAN exactly
  and recovered Sonnet's 0.4/0.4/0.2 mixture weights. It also flagged BLYP's 19% deviation from
  the uniform-gas limit, which is the known LYP defect.
- **Compute contention mattered.** Three agents plus two reference runs on 4 cores; BLAS threads
  were not pinned, so load averages hit 16–19 and each QZ dev run took ~10 min. Round 2 pins
  OMP/OpenBLAS threads.

## Rule changes for round 2
1. Prompt rewritten around "new constraints / new forms / ML with constraints"; wrappers, mixes
   and re-fits are disqualified explicitly, and the dev set is declared a smoke test only.
2. audit.py is given to contestants and its verdict gates scoring (only DISTINCT is scored).
3. Hidden set B (second-row AEs, barrier height, noncovalent, proton affinities) decides ranking.
4. NOTES.md must state what is new, the physical argument, and which exact conditions are
   satisfied/violated; a reviewer reads it.
5. Harness: HYB bug fixed; thread pinning in the test command; 60-minute budget.
