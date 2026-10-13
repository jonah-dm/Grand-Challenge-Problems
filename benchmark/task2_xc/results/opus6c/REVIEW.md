# Organiser review of opus6c (round 6, construction track)

Functional: LC hybrid, 100% LR exact exchange, 0% SR exact exchange; SR semilocal exchange with a new
one-orbital branch e = −(π/2ω²) n_s² · 5Q(b), b = |∇n_s|²/(4ω²n_s²), Q from the exact SR self-exchange
fraction of hydrogenic ions (closed form via erfcx), weight ζ²z_s²; rest ITYH-attenuated sonnet3
exchange; opus4 correlation. ω* = 0.5767 from the one-orbital bound F ≤ 1.174 on H (derived).

Dev WTMAD 1.21 (AE 2.08, best of any new functional; IP 8.27, worst class). All converged. Audit
DISTINCT (UEG-deviation flag is an artifact of the audit assuming a global hybrid; noted).

Targeted quantities (agent-reported, not independently recomputed):
- H₂⁺ at R = 6 bohr, self-consistent: +0.09 kcal/mol (base −28.5; opus5t −6.1).
- Fractional H f(½) = 0.2506 (exact 0.25; LC-ωPBE 0.296; ωB97X 0.304).
- H+H₂ barrier 10.42 (exact 9.6). N₂ −0.5, HCN −0.4, F₂ −2.2 (opus5t: N₂ +7.7).
- Remaining: H₂⁺ at R = 2 bohr −5.1 (shape error); Gaussian λ=1 shape error +5.5%; IP(B) +9.9, IP(F) +10.4.

Claims (to literature check, results/round6_literature_check.md): (a) with erfc screening the
one-orbital SR branch is amplitude-invariant in b, hence exactly λ² and bounded (max = on-top SR
exchange); (b) finiteness + amplitude condition force 100% LR exact exchange in one-orbital
channels; (c) the O(ω⁻⁴) term needs no Laplacian after integration by parts; (d) ω* from the
one-orbital bound. The R function is credited to opus5c's Corollary 4.

Disclosures: seven candidates listed with scores; SR-HF fraction (1/4 vs 0) chosen on dev SVP
score; cap changed from 1.804 to 1.174 after phase-2 runs because of H-density flattening.

Verdict: best targeted construction of the benchmark (delocalization error on one-electron
systems essentially removed with derived constants); overall accuracy not competitive because of
IPs; novelty pending literature check.
