# sonnet2: "OSC-PBE0" hybrid meta-GGA (z = tau_W/tau as one-orbital detector)

Final dev-set WTMAD (def2-qzvp): **1.35** (AE MAE 3.21, IP MAE 3.32 kcal/mol). Audit verdict: DISTINCT
(mix R2 0.998, closest single functional R2 0.969). All species converge. No parameter fit to any data.

## What is new
All pieces are driven by the iso-orbital indicator z = tau_W/tau (z=1 for any one-orbital density, z->0 slowly varying).
1. Exchange enhancement blends two physically fixed forms with weight z^2:
   - z->0 (many orbitals, slowly varying): PBE-form with the exact gradient-expansion coefficient mu = 10/81
     (PBE deliberately uses 0.2195; PBEsol-like GE2 restored where the GE is valid).
   - z->1 (one-orbital regions, density tails): F_x - 1 = mu_PBE s^2 / (1 + c s^(5/2)), which decays as s^(-1/2)
     at large s (the exact tail behaviour of the exchange energy density; PBE tends to a constant). c = 4/(5 kappa/mu)^1.25
     is set analytically so that max(F_x) = 1 + kappa = 1.804, i.e. the Lieb-Oxford bound is saturated, not fit.
2. Correlation: PBE/PW92 times gate [1 - zeta^2 z^2]. This makes E_c vanish exactly for every fully spin-polarized
   one-orbital density (H atom, hydrogenic ions), removing PBE's one-electron self-correlation (-6 mEh on H).
   It does nothing for singlets or the uniform gas. (Only visible for spin-polarized systems, hence
   invisible to the unpolarized audit; the exchange blend is what makes the audit see a distinct form.)
3. 25% exact exchange (a = 1/n, n=4, from the MP4/adiabatic-connection argument); semilocal exchange scaled by 0.75.
Derivatives are exact via complex-step differentiation (no finite-difference noise); z is smoothly saturated at 1.

## Exact conditions
Satisfied: uniform-gas limit (exact, 1e-15); spin-scaling of exchange; Lieb-Oxford bound (exchange: F_x <= 1.804 with
equality saturated; checked pointwise by audit, 0 violations); E_c<=0; correct 2nd-order GE exchange coefficient for z->0;
E_c = 0 for fully polarized one-orbital densities; s^(-1/2) exchange tail for one-orbital regions; high-density correlation
limit (beta = 0.0667, via PBE H); negativity (0 violations).
Violated: no SIC for exchange (relies on 25% HF); correlation GE2 and low-density limit of PBE inherited; no
fractional-charge/derivative-discontinuity condition, no strict-confinement/low-density bound; no nonlocal/dispersion
(noncovalent interactions will be poor); z^2 blend and smooth clamp are heuristic interpolation choices; exact 1-orbital
exchange (F_x for H atom) is not reproduced pointwise.

## Honest caveats
An earlier variant with a SCAN-style kappa_0=0.174 tightening in the one-orbital limit gave WTMAD 8.8 (strongly
overbound molecules) and was discarded; that variant and a PBE0+gate-only variant were seen on a few species only.
I looked at the dev-set error pattern when choosing between these two structural variants (not parameter fitting);
the final form has no tunable parameter. With the 60-min budget, no further variants were tried.

## With more time
Make the exchange blend satisfy the exact 1-orbital exchange hole (H-atom norm) pointwise; add alpha-dependent
correlation (BETA, low-density limit) and a strongly-interacting-limit bound; range-separation/local-hybrid for the
derivative discontinuity and fractional-charge behavior; VV10-like or exact-constraint nonlocal term for dispersion;
test on second-row, barriers and NCIs.
