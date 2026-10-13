# Prior art index (what each round found, and what was already known)
- sonnet2: correlation × [1−ζ²z²] (= TPSS-style one-electron SIE-free gate, 2003). Dev 1.35.
- opus2: constraint-only exchange with integrated GE2 (= SCAN's handling). Dev 3.23.
- haiku3: z-interpolated exchange, H-atom norm (= SCAN norm idea). Dev 1.94.
- sonnet3: Stoll same-spin removal in one-orbital channels (= B95/M06 correlation idea, 1996). Dev 1.00. BEST.
- opus3: joint one-electron condition on the hybrid sum (sound argument; SCAN-style norms). Dev 1.59.
  Lemma: no PBE-shaped iso-orbital F_x meets both H and He norms.
- haiku4: recombination of sonnet3 + haiku3. Dev 1.30.
- sonnet4: sonnet3 + He-normed opposite-spin rescale. Dev 1.30. Found a two-norm (H+He) exchange
  with s^-1/2 tail (different shape from opus3's lemma) that overbinds ~20 kcal/mol.
- opus4: pair-aware same-spin removal, constant-free. Dev 1.08. Diagnosed sonnet3's 26% missing
  correlation in two-electron singlets from self-computed He/H2 correlation energies.
None of the above is new relative to KNOWN_CONDITIONS.md. Hidden-set results: results/hiddenA, hiddenB.

## Hidden set B (unseen classes) — mean of class MAEs, kcal/mol (lower is better)
wB97M-V 1.41 | opus2 1.80 | B3LYP 2.12 | sonnet2 2.22 | sonnet3 2.22 (1 species unconverged) | opus3 2.62 | r2SCAN 2.63 | haiku3 3.20
Per class (AE2 / barrier / noncovalent / proton affinity): wB97M-V 2.71/0.70/0.17/2.06; opus2 4.49/1.11/0.67/0.95;
B3LYP 1.92/5.30/0.62/0.63; sonnet3 4.32/3.65/0.48/0.42.
LESSON: dev-set rank did not transfer. The constraint-only exchange (opus2), worst on the dev set, is the
best new functional on unseen classes; forms selected on the dev set regressed there. Second-row AEs are
where every new functional loses to B3LYP; barrier heights are where self-interaction work wins.

## Round 5 (theory / construction tracks, 3 h)
- opus5t, opus5c, sonnet5c independently derived the same result: for a fractional one-orbital
  density λn₁, exact E_x = −λ²J[n₁]; a Laplacian-free semilocal exchange satisfies this iff
  F(s, z=1) = c·s⁻² (c = 5/9 from H); that form is unbounded and shape-dependent (2.5× off for a
  Gaussian, +126% for H₂⁺ at 2 bohr), so NO (n, ∇n, τ) functional can be one-electron
  delocalization-free; sonnet5c's LP gives an unavoidable ≥11.6% exchange error at λ∈{1,.75,.5}
  for any bounded F(s). Escape: the Laplacian invariant w = n∇²n/|∇n|² makes F = s⁻²G(w) exact on
  the amplitude ray for any G (sonnet5c: few-% errors on H, H₂⁺, HeH²⁺ exact data; opus5t: closed
  form exact for all fractional hydrogenic densities). Harness cannot pass the Laplacian.
- Bounded partial enforcement (opus5t, F_q = F_LO/(1+F_LO s²/c_q), weight ζ²z_s²): H₂⁺ at 6 bohr
  −28.5 → −6.1 kcal/mol, H+H₂ barrier 7.2 → 9.06 (exact 9.6); dev 1.55 (polarized atoms raised).
- LITERATURE CHECK (results/round5_literature_check.md): (1) is KNOWN (Zhang & Yang 1998; MCY 2006);
  (2) PARTIALLY KNOWN, stronger form in Santra & Perdew, JCP 150, 174106 (2019); the s⁻² closed form
  and Gaussian counterexample are the new parts; (3) the 11.6% LP floor is a GRID ARTIFACT (→0% on fine
  grids), do not cite it; (4) Laplacian escape PARTIALLY KNOWN (Becke–Roussel 1989; Ramasamy et al.
  PRB 112, L161112 (2025)) and only fixes amplitude scaling, not the λ=1 nonlocality obstruction.
- Harness bug found by opus5t (zeros vlapl displacing vtau); fixed; all META scores being redone.
- Dev scores: opus5t 1.55, sonnet5c 1.53, opus5c 1.94.

## Round 6
- opus6t (theory): occupation-square hole-triangle bounds; T3 (flat-plane axis forces s⁻² correlation
  on hydrogenic densities); T4 (nowhere-positive semilocal correlation + flat-plane axis ⇒ 10–33×
  over-correlation of two-electron singlets), verified against an exact 1D Lieb-functional model.
  Direct enforcement wrecks open-shell atoms (negative result). Functional B2 (prior-art pieces) dev 1.37.
- opus6c (construction): LC hybrid, 100% LR-HF, 0% SR-HF; SR one-orbital branch exactly λ² and
  bounded under erfc screening (b = |∇n_s|²/(4ω²n_s²) is amplitude-invariant); ω* = 0.5767 from the
  one-orbital bound. H₂⁺ at 6 bohr +0.09 kcal/mol (base −28.5); f(½) = 0.2506; N₂ −0.5. Dev 1.21:
  AE MAE 2.08 (best new) but IP MAE 8.27 (IP(B) +9.9, IP(F) +10.4). Literature check pending.
