# Chemical Accuracy at Scale

> **In one line:** Make density functional theory both *accurate enough* to trust (≈1 kcal/mol) and *cheap enough* to apply to the systems that matter (10⁴–10⁶ atoms), and build the theory that tells us why it works.

---

## Background

Density functional theory (DFT) is the workhorse of computational chemistry and materials science. It is used in tens of thousands of papers a year and underpins nearly every large materials database (Materials Project, OQMD, AFLOW, NOMAD) and most of the training data for machine-learned interatomic potentials. DFT is how we routinely predict:

- **Energetics** – formation enthalpies, reaction energies and barriers, adsorption energies, phase stability, defect formation energies
- **Structure** – lattice constants, bond lengths, elastic constants, phonons
- **Electronic properties** – band structures, band gaps, densities of states, work functions, magnetic order
- **Response and dynamics** – dielectric and optical response, ab initio molecular dynamics, transport

DFT is exact in principle: Hohenberg–Kohn and Kohn–Sham showed that the ground-state energy is a functional of the electron density alone. In practice, everything we don't know is pushed into a single term, the **exchange–correlation (XC) functional**, which has to be approximated. That leaves two bottlenecks:

1. **Accuracy.** Widely used functionals (LDA, PBE, SCAN/r²SCAN, B3LYP, HSE) typically have errors of 2–10 kcal/mol on thermochemistry and can fail qualitatively on band gaps, charge-transfer states, strongly correlated systems (transition-metal oxides, f-electron materials), dispersion-bound systems and reaction barriers. Most of these failures come from known pathologies: self-interaction/delocalization error, static correlation error, and a missing derivative discontinuity. **Chemical accuracy**, usually taken as **1 kcal/mol (≈ 43 meV per formula unit / reaction)**, is the level where predictions can reliably rank competing structures, mechanisms and candidate materials. Today it generally needs wavefunction methods such as CCSD(T), which scale as O(N⁷).

2. **Cost.** Standard Kohn–Sham DFT scales as **O(N³)** with system size because of diagonalization and orbital orthogonalization. Hybrid functionals, which are often required for acceptable accuracy, scale as O(N⁴) without special treatment. In practice this limits routine calculations to a few hundred to a few thousand atoms, well short of what grain boundaries, amorphous materials, electrochemical interfaces, disordered alloys, battery cathode–electrolyte interphases and biomolecular systems need.

Kohn's **nearsightedness principle** says local electronic properties depend only weakly on distant perturbations, so linear scaling should be achievable in principle. Linear-scaling methods such as ONETEP, CONQUEST, BigDFT, SIESTA/order-N, divide-and-conquer, density-matrix purification and Fermi-operator expansions do exist. In practice, though, they tend to give up accuracy, become unreliable for metals and small-gap systems (where the density matrix decays slowly), or are limited to semilocal functionals.

Meanwhile, **machine-learned functionals** (e.g. DM21, Skala) and **constraint-satisfying non-empirical functionals** (e.g. SCAN) have shown that substantial accuracy gains are possible. Still, no functional is yet at chemical accuracy across main-group chemistry, transition metals and solids at the same time, and the relationship between data-driven fitting and exact physical constraints remains poorly understood.

**The grand challenge is to close both gaps at once:** chemical accuracy at a cost that scales linearly with system size, grounded in theory strong enough to tell us when and why it can be trusted.

---

## Challenge Tracks

### Track 1 — Efficiency: Linear scaling without loss of accuracy

**Prompt.** Find approximations to Kohn–Sham DFT that scale linearly, O(N), with system size without losing chemical accuracy relative to the reference cubic-scaling calculation with the same functional. Implement this in an open-source electronic-structure code *(see Open decisions)*.

**Why it is hard**
- Metals and narrow-gap semiconductors: at zero temperature the density matrix decays only algebraically, so locality-based truncation fails.
- Hybrid and higher-rung functionals: exact exchange is nonlocal and dominates cost.
- Prefactors: many O(N) methods only beat O(N³) above a crossover of thousands of atoms.
- Controlled error: truncation errors must be systematically convergeable and accompanied by an error estimate, not tuned case by case.

**Success criteria**
- Demonstrated O(N) (or O(N log N)) wall-time and memory scaling over at least two decades of system size, ideally up to ≥ 10⁵ atoms.
- Total-energy differences (reaction/formation/relative-phase energies) within **≤ 1 kcal/mol** (or ≤ 1 meV/atom for solids) of the conventional O(N³) result for the same functional.
- Forces accurate enough for stable, energy-conserving molecular dynamics.
- Works on a test set covering **insulators, semiconductors and metals**, including at least one hybrid functional.
- Crossover with conventional DFT at ≤ ~1,000 atoms on standard hardware (CPU or GPU).

**Possible directions:** localized/Wannier-like representations with adaptive truncation, finite-temperature smearing combined with Chebyshev/rational Fermi-operator expansions, subsystem/embedding DFT, orbital-free DFT with machine-learned kinetic energy functionals, low-rank/tensor-hypercontraction treatments of exact exchange, and ML-accelerated density-matrix prediction.

---

### Track 2 — Accuracy: A better exchange–correlation functional

**Prompt.** Design an exchange–correlation functional, whether parametrized, machine-learned, derived from first principles or obtained by other means, that performs better than **[reference functional(s) — see Open decisions]** on **[target properties]** across **[benchmark set]**.

**Target properties (suggested)**
| Domain | Properties |
|---|---|
| Molecular thermochemistry | Atomization energies, reaction energies, barrier heights, non-covalent interactions |
| Solids | Formation enthalpies, lattice constants, bulk moduli, cohesive energies, relative phase stability |
| Electronic structure | Band gaps, ionization potentials / electron affinities, fundamental vs. optical gaps |
| Hard cases | Transition-metal complexes and oxides, magnetic ordering energies, stretched bonds, charge-transfer systems |

**Suggested benchmark sets:** GMTKN55 and MGCDB84 (molecules), W4-11 (high-accuracy thermochemistry), TMC/transition-metal sets, and solid-state sets of experimental or CCSD(T)/QMC-quality lattice constants, cohesive energies, formation enthalpies and band gaps. A **held-out test set** that functional developers have not seen should be used for final evaluation.

**Success criteria**
- Weighted mean absolute error (e.g. WTMAD-2 on GMTKN55) below the best existing functional on the same rung of Jacob's ladder, **and** no worse than a stated tolerance on any individual property class (no "robbing Peter to pay Paul").
- Transferability: the improvement holds on held-out materials and chemistries absent from training.
- Self-consistent, numerically stable, with smooth potentials and practical grid requirements; usable in at least one standard code.
- Cost no more than a stated multiple (e.g. ≤ 2×) of the reference functional of the same rung.
- Bonus: satisfies as many known exact constraints as possible (see Track 3), and reports which constraints it violates.

---

### Track 3 — Theory: Exact conditions on the functional

The exact XC functional is unknown, but many of its properties are known: coordinate-scaling relations, the uniform-electron-gas limit, the gradient expansion, the Lieb–Oxford bound, exchange- and correlation-hole sum rules, piecewise linearity in fractional electron number and the flat-plane condition, the derivative discontinuity, and the asymptotic behaviour of the XC potential. Non-empirical functionals such as PBE and SCAN are built entirely from these constraints, and they are the main reason such functionals transfer well.

**Prompt 3a — Tighter bounds and higher-order expansions.** Find tighter bounds or higher-order expansions of known properties of the exact functional. Examples:
- Close the gap between the best rigorous Lieb–Oxford constant and the value conjectured to be tight; establish tight bounds for specific classes of densities (e.g. two-electron, spin-polarized, low-dimensional).
- Derive higher-order terms in the gradient expansion of exchange and correlation, or in the semiclassical (large-Z) expansion of atomic energies.
- Strengthen bounds relating the kinetic, exchange and correlation energies (e.g. Lieb–Thirring-type inequalities, the T_c ≤ −E_c relationships, adiabatic-connection convexity).

**Prompt 3b — New exact sum rules.** Find new exact sum rules or constraints that the exact XC functional must satisfy, and show that enforcing them in an approximate functional improves ground-state properties.

**Success criteria**
- A rigorous proof, or a conjecture backed by convincing numerical evidence from exact or near-exact calculations (e.g. FCI, DMRG, QMC, Lieb-optimized Kohn–Sham inversions).
- Demonstration that common approximations violate the new condition, with an analysis of which errors that explains.
- For 3b: an approximate functional modified to satisfy the new constraint that shows measurable, transferable improvement on Track 2 benchmarks.

---

### Track 4 — Extensions: Beyond ground-state Kohn–Sham DFT

**Prompt.** Extend the results of Tracks 1–3 (linear-scaling algorithms, improved functionals, exact constraints) to generalizations of DFT, including but not limited to:

| Generalization | Example open problems |
|---|---|
| **Generalized Kohn–Sham (GKS)** | Rigorous meaning of GKS band gaps; constraints for hybrids and range-separated hybrids; optimal tuning without system-specific fitting |
| **Time-dependent DFT (real-time and linear-response)** | Memory (non-adiabatic) kernels; charge-transfer and double excitations; excitons in solids; linear-scaling real-time propagation |
| **Relativistic DFT** | Relativistic corrections to XC; spin–orbit coupling in heavy-element and topological materials; four-component vs. two-component consistency |
| **Superconducting DFT (SCDFT)** | Predictive, parameter-free T_c; better XC kernels for the anomalous density; treatment of unconventional pairing |
| **Spin and current DFT** | Noncollinear magnetism; spin-current functionals; orbital magnetism |
| **Finite-temperature / ensemble DFT** | Temperature-dependent XC for warm dense matter; ensemble functionals for excited states and fractional occupation |
| **Embedding and multicomponent DFT** | Systematic DFT-in-DFT/WF-in-DFT embedding; nuclear–electronic (NEO) DFT for proton quantum effects |

**Success criteria:** analogous to Tracks 1–3 for the chosen generalization, namely a demonstrable improvement in cost or accuracy against an accepted reference, or a rigorous new exact condition with evidence that it matters in practice.

---

## Evaluation and Rules (draft)

- **Reference data:** CCSD(T)/CBS, W4-quality composite methods, diffusion Monte Carlo, or carefully curated experiment (with zero-point and thermal corrections applied consistently).
- **Reproducibility:** submissions include code, inputs and outputs; results must be reproducible by the organizers.
- **Hidden test set:** final ranking for Track 2 (and accuracy checks for Track 1) uses a held-out set.
- **Data hygiene:** submitters disclose all training data; no training on the hidden test set or close duplicates.
- **Compute reporting:** wall time, memory and hardware are reported for all Track 1 and Track 2 results.

---

## Open Decisions (to be filled in)

- [ ] **Track 1:** Which code(s) must the implementation target? Candidates include PySCF, CP2K, Quantum ESPRESSO, GPAW, ONETEP, BigDFT, or "any open-source code".
- [ ] **Track 1:** Which functional(s) are required: semilocal only, or must hybrids be supported?
- [ ] **Track 2:** Reference functional(s) to beat (the "XXX"), e.g. ωB97M-V for molecules, r²SCAN or HSE06 for solids, or the best functional per rung.
- [ ] **Track 2:** Exact property list and benchmark/held-out sets, and how the scoring metric weights them.
- [ ] **Track 2:** Are ML functionals and non-local (e.g. full density-matrix) functionals in scope? What are the cost limits?
- [ ] **All tracks:** Prize structure, milestones (e.g. partial prizes for crossing intermediate thresholds) and timelines.

---

## Key References

- P. Hohenberg & W. Kohn, *Phys. Rev.* **136**, B864 (1964); W. Kohn & L. J. Sham, *Phys. Rev.* **140**, A1133 (1965).
- W. Kohn, "Density functional and density matrix method scaling linearly with the number of atoms," *Phys. Rev. Lett.* **76**, 3168 (1996).
- E. Prodan & W. Kohn, "Nearsightedness of electronic matter," *PNAS* **102**, 11635 (2005).
- J. P. Perdew & K. Schmidt, "Jacob's ladder of density functional approximations," *AIP Conf. Proc.* **577**, 1 (2001).
- J. Sun, A. Ruzsinszky & J. P. Perdew, "Strongly constrained and appropriately normed semilocal density functional," *Phys. Rev. Lett.* **115**, 036402 (2015).
- A. J. Cohen, P. Mori-Sánchez & W. Yang, "Insights into current limitations of density functional theory," *Science* **321**, 792 (2008).
- L. Goerigk et al., "A look at the density functional theory zoo with the advanced GMTKN55 database," *PCCP* **19**, 32184 (2017).
- N. Mardirossian & M. Head-Gordon, "Thirty years of density functional theory in computational chemistry," *Mol. Phys.* **115**, 2315 (2017).
- J. Kirkpatrick et al., "Pushing the frontiers of density functionals by solving the fractional electron problem," *Science* **374**, 1385 (2021).
- E. Runge & E. K. U. Gross, *Phys. Rev. Lett.* **52**, 997 (1984); M. Lüders et al. / M. A. L. Marques et al., *Phys. Rev. B* **72**, 024545 / 024546 (2005).
