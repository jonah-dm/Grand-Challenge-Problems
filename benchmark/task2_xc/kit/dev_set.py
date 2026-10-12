"""Public development set. Geometries are fixed inputs (approx. experimental, Angstrom).

Each species: name -> (atoms, charge, spin) where spin = N_alpha - N_beta.
Each item: (class, expression) where expression is {species: coefficient};
the target value is sum(coef * E[species]) in kcal/mol.
"""
import numpy as np
from geom import linear, diatomic, bent, pyramid, tetra, d2r

ATOM_SPIN = {"H": 1, "Li": 1, "B": 1, "C": 2, "N": 3, "O": 2, "F": 1}

SPECIES = {f"{a}_atom": ([(a, (0, 0, 0))], 0, s) for a, s in ATOM_SPIN.items()}
SPECIES.update({
    "H2": (diatomic("H", "H", 0.741), 0, 0),
    "LiH": (diatomic("Li", "H", 1.595), 0, 0),
    "CH4": (tetra("C", "H", 1.087), 0, 0),
    "NH3": (pyramid("N", "H", 1.012, 106.7), 0, 0),
    "H2O": (bent("O", "H", 0.958, 104.5), 0, 0),
    "HF": (diatomic("F", "H", 0.917), 0, 0),
    "N2": (diatomic("N", "N", 1.098), 0, 0),
    "CO": (diatomic("C", "O", 1.128), 0, 0),
    "F2": (diatomic("F", "F", 1.412), 0, 0),
    "C2H2": (linear(("H", 0), ("C", 1.063), ("C", 1.203), ("H", 1.063)), 0, 0),
    "HCN": (linear(("H", 0), ("C", 1.066), ("N", 1.153)), 0, 0),
    "H2CO": ([("C", (0, 0, 0)), ("O", (0, 0, 1.205)),
              ("H", (1.111 * np.sin(58.05 * d2r), 0, -1.111 * np.cos(58.05 * d2r))),
              ("H", (-1.111 * np.sin(58.05 * d2r), 0, -1.111 * np.cos(58.05 * d2r)))], 0, 0),
    "OH": (diatomic("O", "H", 0.970), 0, 1),
    "CH2_triplet": (bent("C", "H", 1.078, 133.9), 0, 2),
    # cations for vertical ionization potentials
    "B_cation": ([("B", (0, 0, 0))], 1, 0),
    "F_cation": ([("F", (0, 0, 0))], 1, 2),
    "HF_cation": (diatomic("F", "H", 0.917), 1, 1),
})


def atomization(name):
    from collections import Counter
    atoms = SPECIES[name][0]
    expr = {f"{s}_atom": c for s, c in Counter(a for a, _ in atoms).items()}
    expr[name] = -1
    return expr


AE = ["H2", "LiH", "CH4", "NH3", "H2O", "HF", "N2", "CO", "F2", "C2H2", "HCN", "H2CO", "OH", "CH2_triplet"]
ITEMS = {f"AE:{m}": ("AE", atomization(m)) for m in AE}
ITEMS.update({
    "IP:B": ("IP", {"B_cation": 1, "B_atom": -1}),
    "IP:F": ("IP", {"F_cation": 1, "F_atom": -1}),
    "IP:HF": ("IP", {"HF_cation": 1, "HF": -1}),
})
