## Literature report: observable-aware compression of Wannier functions into GTOs

1. Wannier-to-GTO fitting
Bakhta, Cancès, Cazeaux, Fang, Kaxiras (2018), "Compression of Wannier functions into Gaussian-type orbitals", Comput. Phys. Commun. (arXiv:1712.02996). The direct ancestor of this project: a greedy algorithm compresses MLWFs (graphene, hBN, FeSe, Si) into symmetry-adapted Gaussian-polynomial orbitals by minimizing a Sobolev-norm approximation error, adding one GTO at a time. The stated motivation is exactly yours — fast evaluation of TB matrix elements — but the loss is purely function-norm; TB observables are not in the objective, and accuracy of t(R) is not controlled.
Zhou, Gull, Zgid (2021), "Material-specific optimization of Gaussian basis sets against plane-wave data", J. Chem. Theory Comput. (arXiv:2105.00589). Optimizes GTO exponents/contractions for a specific solid by simultaneously minimizing total energy and band-energy deviations from a plane-wave reference, with an explicit penalty on the overlap-matrix condition number. This is the closest published analogue of a "mixed loss": function/energy term + derived-observable term + conditioning term.
pySCF↔Wannier90 interoperability (pyWannier90, Pham 2019+) goes the other way (build MLWFs from GTO calculations); nobody there fits GTOs to plane-wave Wannier functions with observable targets.

Verdict: fitting itself is established (your baseline paper); an observable-augmented loss for the Wannier→GTO direction appears not to exist.

2. Basis optimization on derived observables (quantum chemistry / NAO codes)
Dunning (1989) correlation-consistent sets: exponents optimized to maximize contribution to correlation energy (an observable), not wavefunction norm — the canonical precedent for observable-targeted basis design. Jensen's polarization-consistent and later polarizability-optimized cc sets (e.g. JCTC 2024 multiresolution benchmark) show that energy-optimized bases converge irregularly for other properties, motivating property-targeted reoptimization — direct conceptual support for your finding that H¹-optimal ≠ observable-optimal.
SIESTA (Junquera, Artacho et al., early 2000s): PAO shapes/cutoffs optimized by simplex minimization of the total energy — variational, observable-targeted, with known transferability bias.
ABACUS / spillage (Chen, Guo, He 2010; Lin, Ren, He 2021): NAOs minimize the "spillage" — norm lost when projecting plane-wave Kohn-Sham states onto the AO basis — later a "generalized spillage" including a gradient term (effectively an H¹-seminorm on the occupied subspace). Notably they too found plain L² spillage gives poor derived quantities and added kinetic/gradient terms — a step in your direction, but still norm-based rather than S(R)/T(R)-based.
BasisOpt (Shaw & Hill 2023, J. Chem. Phys.): general Python framework where the objective is user-definable (energies, properties) — "property-targeted basis optimization" as a tool exists.

Verdict: optimizing bases on energies/properties rather than function norms is standard practice in quantum chemistry; your specific targets (real-space hopping/overlap tables) are a known variant of this philosophy, not previously instantiated.

3. TB downfolding for TBG
Fang & Kaxiras (2016), PRB 93, 235153 (arXiv:1604.05371): DFT + Wannier transformation of bilayer graphene at sampled stackings; interlayer coupling fitted to a smooth functional form t(r, θ) with angular harmonics (beyond two-center Slater-Koster) reflecting higher angular-momentum content of the Wannier p_z. Reference data = Wannier matrix elements themselves; they fit hopping tables, not orbitals.
Moon & Koshino (2012/2013), PRB: two-center Slater-Koster exponential t(d) = V_ppπ, V_ppσ parametrization with decay length 0.319a₀ — fitted to reproduce bands, universally used for TBG.
Koshino, Yuan, Koretsune, Ochi, Kuroki, Fu (2018), PRX 8, 031087: MLWFs of the continuum model (three-peak orbitals) → extended Hubbard model; again the orbitals are an output, hoppings tabulated numerically, no analytic orbital fit.
Carr, Fang, Po, Vishwanath, Kaxiras (2019) (arXiv:1907.06282): multistep Wannier projection from k·p to minimal TB models with relaxation; hoppings tabulated.
Li & Chan (2026), "Intrinsic Wannier functions for Hamiltonian downfolding" (arXiv:2608.15557): non-iterative downfolding judged on downfolded-Hamiltonian quality — sign that the community increasingly evaluates Wannierization by the Hamiltonian it produces, not orbital localization.

Verdict: TBG downfolding fits functional forms for t(R) against Wannier/DFT matrix elements or bands; no one fits analytic localized orbitals whose induced S(R), T(R) match reference — that is your gap.

4. Observable-targeted compression
White (2017); White & Stoudenmire (2019) gausslets / multisliced gausslet bases: system-adapted Gaussian-wavelet bases engineered so that Hamiltonian matrix elements (notably the two-electron tensor) acquire diagonal, compressed structure — compression judged by operator representation, not state fidelity. Philosophically aligned but a different mechanism (basis construction, not least-squares fitting).
SCDM (Damle, Lin, Ying 2015, arXiv:1408.4926): compressed representation of Kohn-Sham orbitals via selected columns of the density matrix — localization by interpolative decomposition; again subspace-fidelity, not observable-targeted; no follow-up found that fits GTOs to SCDM orbitals with matrix-element losses.
The Zhou–Gull–Zgid paper (area 1) is the best example of compression optimized against band energies (spectra of the downfolded operator).

Verdict: "compress so the operator's matrix elements are right" exists as a design philosophy (gausslets, material-specific GTOs), but a least-squares loss mixing H¹ error with penalties on ⟨w_0|A|w_R⟩ appears novel in this literature.

5. Orthonormality-constrained fitting
Löwdin (1950) symmetric orthogonalization, with the classical result (Carlson–Keller; Mayer 2002, Int. J. Quantum Chem.) that S^(-1/2) is the orthonormal set closest in L² to the original — the standard post-hoc route rather than an in-loop penalty.
SIESTA/ABACUS/OpenMX NAOs are per-atom orthonormal by construction; inter-site non-orthogonality is kept (generalized eigenproblem) rather than penalized. Direct minimization on the Stiefel manifold (e.g. arXiv:2412.18807 and the orbital-minimization-method literature, Corsetti 2014) is the standard way to impose orthonormality as a hard constraint during optimization; soft orthonormality penalties are common in machine-learning contexts but rare in orbital fitting.

Verdict: orthonormality handling is standard (Löwdin post-processing or manifold constraints); a soft penalty is a known, unremarkable variant. Consider Löwdin-then-fit or Stiefel-constrained optimization as alternatives to a penalty weight you must tune.

Overall conclusion

Your planned loss — H¹ approximation error + penalties on {S(R), T(R)} + orthonormality — is a novel combination of standard ingredients. Each ingredient has precedent: observable-targeted basis optimization (Dunning; Zhou–Gull–Zgid), norm-based Wannier→GTO compression (Bakhta et al.), t(R) fitting against Wannier references (Fang–Kaxiras). What appears genuinely absent from the literature is fitting localized analytic orbitals under a loss that directly penalizes errors in the TB observables they induce. The ABACUS experience (L² spillage → add gradient term) and the polarizability-basis literature independently confirm your diagnosis that norm-optimality underweights what observables probe.

Top 5 most relevant papers:

Bakhta, Cancès, Cazeaux, Fang, Kaxiras 2018 (arXiv:1712.02996) — your baseline, pure-norm loss.
Zhou, Gull, Zgid 2021 (arXiv:2105.00589) — closest existing mixed loss (energy + bands + conditioning).
Fang & Kaxiras 2016, PRB 93, 235153 — reference t(R) functional forms for graphene bilayers.
Koshino et al. 2018, PRX 8, 031087 — TBG Wannier/Hubbard downfolding standard.
Lin, Ren, He 2021 (ABACUS generalized spillage) — norm-based orbital fitting augmented with derivative terms, the nearest NAO analogue.

Sources: 
arXiv:1712.02996,
ScienceDirect CPC (https://www.sciencedirect.com/science/article/abs/pii/S0010465518301231),
arXiv:2105.00589, 
arXiv:1604.05371, 
PRX 8, 031087,
arXiv:1907.06282,
arXiv:2608.15557,
SIESTA basis optimization, 
OpenMX NAO note, 
BasisOpt JCP (https://pubs.aip.org/aip/jcp/article/159/4/044802/2904098/BasisOpt-A-Python-package-for-quantum-chemistry), 
JCTC polarizability benchmark, 
Multisliced gausslets, 
SCDM arXiv:1408.4926,
Löwdin minimum-perturbation JCP,
Stiefel-manifold DFT arXiv:2412.18807, 
ABACUS arXiv:2501.08697, 
pyWannier90
