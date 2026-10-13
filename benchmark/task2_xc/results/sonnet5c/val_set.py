"""Own validation set (NOT the dev set, NOT hidden): classes the dev set lacks. References: CCSD(T)/CBS via
reference.py machinery (cc-pVTZ/QZ). Used only for the disclosed candidate comparison."""
import numpy as np, sys
sys.path.insert(0, '/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit')
from geom import linear, diatomic, bent, pyramid, tetra, d2r
SPECIES = {}
for a, s in {"H": 1, "P": 3, "S": 2, "Cl": 1}.items():
    SPECIES[f"{a}_atom"] = ([(a, (0, 0, 0))], 0, s)
SPECIES.update({
    "H2": (diatomic("H", "H", 0.741), 0, 0),
    "HCl": (diatomic("Cl", "H", 1.275), 0, 0),
    "H2S": (bent("S", "H", 1.336, 92.1), 0, 0),
    "PH3": (pyramid("P", "H", 1.42, 93.3), 0, 0),
    "H3_TS": (linear(("H", 0), ("H", 0.930), ("H", 0.930)), 0, 1),
    "NH3": (pyramid("N", "H", 1.012, 106.7), 0, 0),
    "NH4+": (tetra("N", "H", 1.022), 1, 0),
    "H2O": (bent("O", "H", 0.958, 104.5), 0, 0),
    "H3O+": (pyramid("O", "H", 0.976, 111.3), 1, 0),
})
def ae(name):
    from collections import Counter
    atoms = SPECIES[name][0]
    expr = {f"{s}_atom": c for s, c in Counter(a for a, _ in atoms).items()}
    expr[name] = -1
    return expr
ITEMS = {f"AE2:{m}": ("AE2", ae(m)) for m in ["HCl", "H2S", "PH3"]}
ITEMS.update({
    "BAR:H+H2": ("BAR", {"H3_TS": 1, "H_atom": -1, "H2": -1}),
    "PA:NH3": ("PA", {"NH3": 1, "NH4+": -1}),
    "PA:H2O": ("PA", {"H2O": 1, "H3O+": -1}),
})
