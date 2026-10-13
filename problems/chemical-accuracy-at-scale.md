Problem Statement: Chemical Accuracy at Scale for Density Functional Theory

Prompt: Develop methods that bring the accuracy of high-level quantum chemistry (≈1 kcal/mol, "chemical accuracy") to the properties that decide real materials choices, at a cost low enough to apply to thousands of systems of thousands of atoms. Density functional theory is exact in principle for the ground-state energy and density [1], but all unknown physics is carried by the exchange–correlation (XC) functional, which must be approximated; current approximations are off by 2–10 kcal/mol and fail qualitatively for strong correlation and charge transfer [2]. Band gaps lie outside ground-state DFT altogether: the Kohn–Sham gap differs from the fundamental gap by the derivative discontinuity even for the exact functional [3], so the theory itself must be extended. And standard Kohn–Sham DFT scales as O(N³), capping routine calculations near 10³ atoms although locality arguments say linear scaling is possible [4].

Given a set of molecules or materials, predict their energetics, electronic structure and response properties to within stated tolerances of CCSD(T), quantum Monte Carlo or experiment, at a cost that scales linearly with the number of atoms, with the accuracy traceable to physics (exact constraints, derived forms, or machine-learned constructions that satisfy exact conditions by construction) rather than to fitting on the test data.

The method need not be a single functional or algorithm. New exact conditions, new functional forms, machine-learned functionals, linear-scaling solvers, embedding schemes and extensions of the theory (generalized Kohn–Sham, time-dependent, relativistic, superconducting DFT) are all acceptable. The goal is a shared, physically grounded description whose accuracy transfers to systems and property classes it was not developed on, not a collection of unrelated fits.

Calibration Systems (public development set)

| System | Required predictions |
|---|---|
| First-row small molecules (H₂, LiH, CH₄, NH₃, H₂O, HF, N₂, CO, F₂, C₂H₂, HCN, H₂CO, OH, ³CH₂) and atoms | Atomization energies to ≤ 1 kcal/mol against CCSD(T)/CBS; vertical ionization potentials of B, F and HF. Baselines (atomization / IP MAE, kcal/mol): PBE 7.3 / 5.2, r²SCAN 2.6 / 6.3, B3LYP 1.7 / 6.4, ωB97M-V 1.3 / 1.3. |
| Exact and near-exact densities (H, He, Li, Be, hydrogenic ions, uniform gas) | Satisfaction of the known exact conditions [5, 6]: uniform-gas limit, scaling relations, gradient expansion, Lieb–Oxford bound, one-electron self-interaction freedom, fractional-charge linearity. Each violated condition must be declared. |
| Bulk silicon and a periodic water box, 64 to 10⁴ atoms | Energy and forces within ≤ 1 meV/atom and ≤ 10 meV/Å of the conventional O(N³) result for the same functional; measured scaling exponent. |

What makes the method general?

Use shared physical descriptions: exact constraints and appropriate norms for the functional, locality for the solver, and a declared set of ingredients. Parameters are allowed only when derived from exact conditions or from data the submitter generates independently of the benchmark, and their origin must be stated before prediction. Do not tune any parameter to the development or hidden sets; selection among a few physically motivated candidate forms is permitted if every candidate and its score is disclosed.

Recover established limits where their assumptions hold (uniform gas, slowly varying gas, one-electron limit), but permit different mechanisms where the physics differs: charge transfer and stretched bonds are governed by fractional-charge behaviour rather than semilocal corrections [3]; band gaps need a generalized Kohn–Sham treatment [7]; metals defeat locality-based truncation at zero temperature [4].

Evidence and validation

The method applied to the calibration systems must transfer, without re-parametrization, to the systems below. Scoring uses hidden sets; a submission is ranked on the property classes absent from calibration and must not be worse than the best existing functional on any class by more than a stated tolerance.

| Validation system | Required predictions |
|---|---|
| Unseen first-row molecules (C₂H₄, CO₂, H₂O₂, O₂, NO, CH₃, NH₂, CH₃OH, HNO, LiF, BH₃, BF) and ionization potentials | Atomization energies and IPs to ≤ 1 kcal/mol, showing calibration accuracy transfers to new molecules of the same classes. |
| Second-row molecules (SiH₄, PH₃, H₂S, HCl, Cl₂, SO₂, CS), the H + H₂ barrier, the water dimer, proton affinities of NH₃ and H₂O | The same tolerances on classes and elements absent from calibration. Existing reference on these four classes: ωB97M-V 2.7 / 0.7 / 0.2 / 2.1 kcal/mol; B3LYP 1.9 / 5.3 / 0.6 / 0.6. |
| Main-group and transition-metal thermochemistry databases (GMTKN55 [8] and transition-metal sets) | Weighted error below the best functional of the same rung, with no subset worse than it by more than the tolerance. |
| Solids: lattice constants, cohesive and formation energies, and fundamental band gaps of semiconductors, insulators, metals and transition-metal oxides | Energies within ≤ 1 kcal/mol per formula unit; band gaps within 0.2 eV through a generalized Kohn–Sham scheme [7]. |
| Large, disordered systems (amorphous SiO₂, a metal slab, a solvated protein), 10⁴–10⁵ atoms | Calibration-level accuracy where the O(N³) reference is computable; measured O(N) scaling over two decades of size; stable molecular dynamics. |
| Extensions: TDDFT excitation energies, relativistic spin–orbit splittings, superconducting T_c | The analogous accuracy or cost metrics against established benchmarks, derived from the same exact conditions. |

Theory track

The challenge also accepts results on the exact functional itself: tighter bounds or higher-order expansions of known properties, and new exact conditions that the exact functional satisfies and existing approximations violate. These are scored by proof or expert review, numerical verification on exact densities, and the improvement obtained on the validation systems when the condition is enforced.

Source notes

[1] P. Hohenberg and W. Kohn, Physical Review 136, B864 (1964), DOI: 10.1103/PhysRev.136.B864; W. Kohn and L. J. Sham, Physical Review 140, A1133 (1965), DOI: 10.1103/PhysRev.140.A1133.

[2] A. J. Cohen, P. Mori-Sánchez and W. Yang, Insights into current limitations of density functional theory, Science 321, 792 (2008). DOI: 10.1126/science.1158722.

[3] J. P. Perdew, R. G. Parr, M. Levy and J. L. Balduz, Density-functional theory for fractional particle number: derivative discontinuities of the energy, Physical Review Letters 49, 1691 (1982). DOI: 10.1103/PhysRevLett.49.1691.

[4] S. Goedecker, Linear scaling electronic structure methods, Reviews of Modern Physics 71, 1085 (1999). DOI: 10.1103/RevModPhys.71.1085.

[5] J. P. Perdew, K. Burke and M. Ernzerhof, Generalized gradient approximation made simple, Physical Review Letters 77, 3865 (1996). DOI: 10.1103/PhysRevLett.77.3865.

[6] J. Sun, A. Ruzsinszky and J. P. Perdew, Strongly constrained and appropriately normed semilocal density functional, Physical Review Letters 115, 036402 (2015). DOI: 10.1103/PhysRevLett.115.036402.

[7] J. P. Perdew, W. Yang, K. Burke et al., Understanding band gaps of solids in generalized Kohn–Sham theory, Proceedings of the National Academy of Sciences 114, 2801 (2017). DOI: 10.1073/pnas.1621352114.

[8] L. Goerigk et al., A look at the density functional theory zoo with the advanced GMTKN55 database, Physical Chemistry Chemical Physics 19, 32184 (2017). DOI: 10.1039/C7CP04913G.
