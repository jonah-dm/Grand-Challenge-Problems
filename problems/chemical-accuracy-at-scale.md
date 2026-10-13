Problem Statement: Chemical Accuracy at Scale for Density Functional Theory

Prompt: Develop methods that bring the accuracy of high-level quantum chemistry (≈1 kcal/mol, "chemical accuracy") to the properties that decide real materials choices, at a cost low enough to apply to thousands of systems of thousands of atoms. Density functional theory (DFT) is exact in principle for the ground-state energy and density [1, 2], but all unknown physics is carried by the exchange–correlation (XC) functional, which must be approximated; today's approximations are off by 2–10 kcal/mol and fail qualitatively for strong correlation and charge transfer [3, 4]. For quantities that are not ground-state properties, above all the band gap, the Kohn–Sham eigenvalue gap differs from the fundamental gap by the derivative discontinuity even for the exact functional [5], so the theory itself must be extended [6, 7]. And standard Kohn–Sham DFT scales as O(N³) with system size, hybrids worse, capping routine calculations near 10³ atoms although locality arguments say linear scaling is possible [8, 9].

Given a set of molecules or materials (composition, geometry, charge and spin state), predict their ground-state energetics, electronic structure and response properties to within stated tolerances of CCSD(T), quantum Monte Carlo or experimental references, with a computational cost that scales linearly with the number of atoms, and with the accuracy traceable to physics (exact constraints, derived forms, or machine-learned constructions that satisfy exact conditions by construction) rather than to fitting on the test data.

The method need not be a single functional or a single algorithm. New exact conditions on the functional, new semilocal, nonlocal or orbital-dependent forms, machine-learned functionals, linear-scaling solvers, embedding schemes and extensions of the theory (generalized Kohn–Sham, time-dependent, relativistic, superconducting DFT) are all acceptable. The goal is a shared, physically grounded description whose accuracy transfers to systems and property classes it was not developed on, not a collection of unrelated fits.

The method should predict the following systems. The outputs below are proposed challenge requirements. References provide the theoretical basis and existing benchmarks for each class; they do not establish that the requested accuracy at the requested cost has already been demonstrated.

Calibration Systems (public development set)

| System | Required predictions |
|---|---|
| First-row small molecules (H₂, LiH, CH₄, NH₃, H₂O, HF, N₂, CO, F₂, C₂H₂, HCN, H₂CO, OH, ³CH₂) and atoms | Atomization energies to ≤ 1 kcal/mol against CCSD(T)/CBS at fixed geometries; vertical ionization potentials of B, F and HF. Baselines on this set: PBE 7.3 / 5.2, r²SCAN 2.6 / 6.3, B3LYP 1.7 / 6.4, ωB97M-V 1.3 / 1.3 kcal/mol (atomization / IP mean absolute error). |
| Exact and near-exact densities (H, He, Li, Be, hydrogenic ions, Hooke's atom, uniform gas) | Satisfaction of known exact conditions [10–13]: uniform-gas limit, spin and coordinate scaling, gradient expansion, Lieb–Oxford bound, one-electron self-interaction freedom, piecewise linearity in fractional charge and the flat-plane condition. Each violated condition must be declared. |
| Bulk crystalline silicon and a periodic water box, 64 to 10⁴ atoms | Total energy and forces within ≤ 1 meV/atom and ≤ 10 meV/Å of the conventional O(N³) calculation with the same functional; measured wall-time and memory scaling exponent. |

What makes the method general?

Use shared physical descriptions: exact constraints and appropriate norms for the functional, locality (nearsightedness) for the solver, and a declared set of ingredients (density, gradient, kinetic-energy density, exact exchange, nonlocal correlation). Material-specific or class-specific parameters are allowed only when they are derived from exact conditions or from data the submitter generates independently of the benchmark; their origin must be stated before prediction. Do not fit a new parameter set to each property class, and do not tune any parameter to the development or hidden sets. Selection among a small number of physically motivated candidate forms using development-set scores is permitted if every candidate and its score is disclosed.

Recover established limits where their assumptions hold (the uniform gas, the slowly varying gas, the one-electron limit), but permit different mechanisms where the physics differs: charge-transfer and stretched-bond systems are governed by the derivative discontinuity and fractional-charge behaviour rather than by semilocal corrections [5, 14]; band gaps require a generalized Kohn–Sham or many-body treatment [6]; metals and small-gap systems defeat locality-based truncation at zero temperature and need finite-temperature or Fermi-operator formulations [9].

Evidence and validation

The method applied to the calibration systems must be applicable, without re-parametrization, to the systems below. Final scoring uses hidden sets of different molecules and property classes; a submission is ranked on the hidden property classes absent from calibration and must not be worse than the best existing functional on any class by more than a stated tolerance.

| Validation system | Required predictions |
|---|---|
| Unseen first-row molecules (e.g. C₂H₄, CO₂, H₂O₂, O₂, NO, CH₃, NH₂, CH₃OH, HNO, LiF, BH₃, BF) and atomic/molecular ionization potentials | Atomization energies and IPs to ≤ 1 kcal/mol against CCSD(T)/CBS, demonstrating that calibration-set accuracy transfers to new molecules of the same classes. |
| Second-row molecules (SiH₄, PH₃, H₂S, HCl, Cl₂, SO₂, CS), reaction barrier heights (H + H₂), noncovalent complexes (water dimer), proton affinities (NH₃, H₂O) | The same tolerances on property classes and elements absent from calibration; self-interaction-sensitive barriers and dispersion-bound complexes must both be captured. Existing reference: ωB97M-V achieves 2.7 / 0.7 / 0.2 / 2.1 kcal/mol on these four classes; B3LYP 1.9 / 5.3 / 0.6 / 0.6. |
| Main-group and transition-metal thermochemistry databases (GMTKN55 [15], transition-metal reaction sets) | Weighted mean absolute error below the best functional of the same rung of Jacob's ladder, with no subset worse than that functional by more than the tolerance. |
| Solids: lattice constants, cohesive energies, formation enthalpies and fundamental band gaps of a fixed set of semiconductors, insulators, metals and transition-metal oxides | Energies within ≤ 1 kcal/mol per formula unit of QMC or curated experiment; band gaps within 0.2 eV of experiment through a generalized Kohn–Sham scheme whose eigenvalue gap is the fundamental gap [6]. |
| Large and disordered systems: amorphous SiO₂, a metal slab, a solvated protein or an electrochemical interface, 10⁴–10⁵ atoms | Energies and forces within the calibration tolerances of the O(N³) reference where that reference is computable, O(N) or O(N log N) measured scaling over at least two decades of system size, and stable energy-conserving molecular dynamics. |
| Extensions: excitation energies (linear-response and real-time TDDFT [7]), spin–orbit splittings in heavy-element systems (relativistic DFT), superconducting transition temperatures (SCDFT [16]) | The analogous accuracy or cost metrics for the chosen generalization against established benchmark sets, with the new functional or algorithm derived from the same exact conditions. |

Theory track

Alongside the predictive method, the challenge accepts results on the exact functional itself: tighter bounds or higher-order expansions of known properties (the Lieb–Oxford constant [12], the gradient expansion, scaling relations), and new exact sum rules or constraints that the exact XC functional satisfies and existing approximations violate. A theory result is scored by proof (formal or expert review), by numerical verification on exact or near-exact densities, and by the improvement obtained on the validation systems when the condition is enforced in a functional.

Source notes

[1] P. Hohenberg and W. Kohn, Inhomogeneous electron gas, Physical Review 136, B864–B871 (1964). DOI: 10.1103/PhysRev.136.B864.

[2] W. Kohn and L. J. Sham, Self-consistent equations including exchange and correlation effects, Physical Review 140, A1133–A1138 (1965). DOI: 10.1103/PhysRev.140.A1133.

[3] A. J. Cohen, P. Mori-Sánchez and W. Yang, Insights into current limitations of density functional theory, Science 321, 792–794 (2008). DOI: 10.1126/science.1158722.

[4] N. Mardirossian and M. Head-Gordon, Thirty years of density functional theory in computational chemistry: an overview and extensive assessment of 200 density functionals, Molecular Physics 115, 2315–2372 (2017). DOI: 10.1080/00268976.2017.1333644.

[5] J. P. Perdew, R. G. Parr, M. Levy and J. L. Balduz, Density-functional theory for fractional particle number: derivative discontinuities of the energy, Physical Review Letters 49, 1691–1694 (1982). DOI: 10.1103/PhysRevLett.49.1691.

[6] J. P. Perdew, W. Yang, K. Burke et al., Understanding band gaps of solids in generalized Kohn–Sham theory, Proceedings of the National Academy of Sciences 114, 2801–2806 (2017). DOI: 10.1073/pnas.1621352114.

[7] E. Runge and E. K. U. Gross, Density-functional theory for time-dependent systems, Physical Review Letters 52, 997–1000 (1984). DOI: 10.1103/PhysRevLett.52.997.

[8] E. Prodan and W. Kohn, Nearsightedness of electronic matter, Proceedings of the National Academy of Sciences 102, 11635–11638 (2005). DOI: 10.1073/pnas.0505436102.

[9] S. Goedecker, Linear scaling electronic structure methods, Reviews of Modern Physics 71, 1085–1123 (1999). DOI: 10.1103/RevModPhys.71.1085.

[10] J. P. Perdew, K. Burke and M. Ernzerhof, Generalized gradient approximation made simple, Physical Review Letters 77, 3865–3868 (1996). DOI: 10.1103/PhysRevLett.77.3865.

[11] J. Sun, A. Ruzsinszky and J. P. Perdew, Strongly constrained and appropriately normed semilocal density functional, Physical Review Letters 115, 036402 (2015). DOI: 10.1103/PhysRevLett.115.036402.

[12] E. H. Lieb and S. Oxford, Improved lower bound on the indirect Coulomb energy, International Journal of Quantum Chemistry 19, 427–439 (1981). DOI: 10.1002/qua.560190306.

[13] P. Mori-Sánchez, A. J. Cohen and W. Yang, Discontinuous nature of the exchange-correlation functional in strongly correlated systems, Physical Review Letters 102, 066403 (2009). DOI: 10.1103/PhysRevLett.102.066403.

[14] P. Mori-Sánchez, A. J. Cohen and W. Yang, Localization and delocalization errors in density functional theory and implications for band-gap prediction, Physical Review Letters 100, 146401 (2008). DOI: 10.1103/PhysRevLett.100.146401.

[15] L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi and S. Grimme, A look at the density functional theory zoo with the advanced GMTKN55 database for general main group thermochemistry, kinetics and noncovalent interactions, Physical Chemistry Chemical Physics 19, 32184–32215 (2017). DOI: 10.1039/C7CP04913G.

[16] M. Lüders, M. A. L. Marques, N. N. Lathiotakis et al., Ab initio theory of superconductivity. I. Density functional formalism and approximate functionals, Physical Review B 72, 024545 (2005). DOI: 10.1103/PhysRevB.72.024545.
