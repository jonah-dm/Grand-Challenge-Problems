# sonnet notes (INCOMPLETE - no valid score)
Idea: linear blend of libxc hybrid meta-GGAs (MN15 0.4, M06-2X 0.4, PW6B95 0.2; HYB = weighted sum 0.448).
Screening with built-in functionals at def2-TZVP (energies mixed linearly, non-self-consistent): WTMAD m062x 0.84, mn15 0.77,
pw6b95 0.89, blend 0.4/0.4/0.2 ~0.67 (def2-tzvp, not the qzvp scoring basis). Would likely not beat wB97M-V (0.56).
BUG (unresolved): submission.py wrapping libxc.eval_xc(name,...) inside define_xc_ gives wrong energies
(H atom -0.36 vs -0.497 built-in), even for a single functional with its own HYB. exc matches the built-in on random rho,
so suspect vxc/vlapl layout or hyb handling in the custom path. No valid score.json was produced.
Next: debug the vxc layout (compare with ni.eval_xc_eff output), then use a range-separated/NLC-free form, add own parameters.
