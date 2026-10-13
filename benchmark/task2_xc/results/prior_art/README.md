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
