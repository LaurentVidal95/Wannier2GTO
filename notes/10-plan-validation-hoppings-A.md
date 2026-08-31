# Plan — validation hoppings, approche A

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Mesurer l'erreur du Wannier gaussien-compressé sur les hoppings
S(R) = ⟨w₀|w_R⟩ et T(R) = ⟨w₀|−½∇²|w_R⟩, R≠0, contre la référence plane-wave
(design : [09-design-validation-hoppings.md](09-design-validation-hoppings.md)).

**Architecture:** Une fonction `compare_hopping` dans
`src/tight_binding/hamiltonian_scalar_product.jl` (jumelle de
`compare_onsite_kinetic` : référence par phase de Fourier e^{−iG·R}, gaussien
par `translate` + `integral`), testée sur formes closes gaussiennes et
identités plane-wave, puis un script `workflow/validate_tb_hoppings.jl` qui
boucle sur les 6 R du design + R=0 (contrôle) pour Ecut ∈ {15, 50}.

**Tech Stack:** Julia, DFTK (basis supercellule), intégrales gaussiennes
natives (`julia_integrals.jl`), données existantes
`data/compressed_wannier_55_H1*.json` et
`workflow/wannier_functions/wannier_pz_Ecut-*.json`.

**Conventions du repo** : tests lancés fichier par fichier
(`julia --project=. test/test_hoppings.jl`), pas de runtests.jl. Commits
proposés à l'utilisateur, **jamais créés sans validation** (règle
research-workflow) — les étapes "Commit" ci-dessous signifient *proposer* le
commit.

**Fait vérifié sur le code** (à ne pas redécouvrir) :
- `translate(Wc, R)` ([CompressedWannier.jl:33](../src/CompressedWannier.jl))
  translate les centres SAGTO de `R` **cartésien en Bohr** ;
- `integral(Wc₁, Wc₂; type=:overlap/:kinetic)`
  ([integrals.jl:63](../src/integrals/integrals.jl)) est natif Julia ;
- `GaussianPolynomial([exps], [1.0], center, ζ)` **normalise** la gaussienne
  (cf. [test_integrals_laplacian.jl](../test/test_integrals_laplacian.jl)) ;
- `BasisFunction(coeffs, SAGTOs)` est un struct nu, pas de normalisation ;
- le constructeur bas niveau
  `CompressedWannier(center, basis_functions, coefficients, nothing, ComplexF64[], ComplexF64[], 0.0, 0)`
  permet de fabriquer un Wc synthétique pour les tests ;
- piège périodique : la référence supercellule contient les images périodiques
  de w_R (5×5), le GTO est en espace libre. Pour |R| ≪ 5a c'est le même
  plancher d'erreur que la compression elle-même — on le note dans le bilan,
  on ne le "corrige" pas.

---

### Task 1: `compare_hopping` (src + tests)

**Files:**
- Test: `test/test_hoppings.jl` (create)
- Modify: `src/tight_binding/hamiltonian_scalar_product.jl` (append after `compare_onsite_kinetic`)
- Modify: `src/Wannier2GTO.jl:67` (add export)

- [ ] **Step 1: Write the failing tests**

`test/test_hoppings.jl` :

```julia
using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra

# Synthetic compressed "Wannier": one normalized s-Gaussian at the origin.
# Closed forms for two normalized s-Gaussians with equal spread ζ, separation R
# (Szabo–Ostrund A.9/A.11 with α=β=ζ):
#   S(R) = exp(-ζ|R|²/2)
#   T(R) = ⟨g|-½∇²|g_R⟩ = (ζ/2)(3 - ζ|R|²) exp(-ζ|R|²/2)
function single_gaussian_Wc(ζ)
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)  # normalized
    Φ = W2G.BasisFunction([1.0], [g])
    W2G.CompressedWannier(zeros(3), [Φ], [1.0],
                          nothing, ComplexF64[], ComplexF64[], 0.0, 0)
end

# Minimal Γ-only plane-wave basis (kinetic-only model, no atoms needed).
function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr, anisotropic on purpose
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "compare_hopping: GTO side against closed forms" begin
    ζ = 0.7
    Wc = single_gaussian_Wc(ζ)
    basis = tiny_basis()
    # Any admissible w_fourier does for this testset: take the Gaussian itself
    # expanded on the basis, so the reference side runs without erroring.
    w = Wc(basis)
    w = w / norm(w)
    for R in ([1.0, 0.0, 0.0], [0.5, -0.3, 0.2], [0.0, 0.0, 1.3])
        out = W2G.compare_hopping(Wc, w, basis, R)
        R² = sum(abs2, R)
        @test out.S_gto ≈ exp(-ζ * R² / 2)                      rtol = 1e-10
        @test out.T_gto ≈ (ζ/2) * (3 - ζ*R²) * exp(-ζ*R²/2)    rtol = 1e-10
    end
end

@testset "compare_hopping: reference side identities" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    n_G = length(G_vectors(basis, kpt))
    # Real random real-space field -> Fourier coefficients with w_{-G} = conj(w_G),
    # i.e. a real "Wannier", as in production.
    w_real = randn(basis.fft_size...)
    w = DFTK.fft(basis, kpt, complex.(w_real))
    w = w / norm(w)
    @assert length(w) == n_G
    Wc = single_gaussian_Wc(0.7)

    out0 = W2G.compare_onsite_kinetic(Wc, w, basis)
    # R = 0 must reproduce the on-site diagnostics exactly.
    outR0 = W2G.compare_hopping(Wc, w, basis, zeros(3))
    @test outR0.S_ref ≈ 1.0            rtol = 1e-12
    @test outR0.T_ref ≈ out0.T_ref     rtol = 1e-12
    @test outR0.T_gto ≈ out0.T_gto     rtol = 1e-12

    # Translation by a lattice vector is the identity on the periodic reference.
    L = basis.model.lattice[:, 1]
    outL = W2G.compare_hopping(Wc, w, basis, L)
    @test outL.S_ref ≈ outR0.S_ref     rtol = 1e-10
    @test outL.T_ref ≈ outR0.T_ref     rtol = 1e-10

    # Hermiticity: S(R) = S(-R), T(R) = T(-R) for a real w.
    R = [1.1, -0.4, 0.6]
    outp = W2G.compare_hopping(Wc, w, basis, R)
    outm = W2G.compare_hopping(Wc, w, basis, -R)
    @test outp.S_ref ≈ outm.S_ref      rtol = 1e-9
    @test outp.T_ref ≈ outm.T_ref      rtol = 1e-9
end
```

Note pour l'ingénieur : si `DFTK.fft(basis, kpt, ...)` n'a pas cette signature
dans la version DFTK du projet, remplacer par la construction directe d'un
vecteur hermitien : tirer `w` complexe aléatoire puis symétriser sur les paires
(G, −G) via `G_vectors(basis, kpt)` — le point du test est w réel en espace
physique, pas la méthode pour l'obtenir.

- [ ] **Step 2: Run tests, verify they fail**

Run: `julia --project=. test/test_hoppings.jl`
Expected: FAIL — `UndefVarError: compare_hopping not defined` (ou
`compare_hopping` pas une propriété de W2G).

- [ ] **Step 3: Implement `compare_hopping`**

Append to `src/tight_binding/hamiltonian_scalar_product.jl` :

```julia
@doc raw"""
Hopping comparison, second link of the tight-binding validation chain:
``S(R) = \langle w_0, w_R\rangle`` and
``T(R) = \langle w_0, -\tfrac12\nabla^2 w_R\rangle`` from the true plane-wave
Wannier vs its Gaussian-compressed counterpart. `R_cart` is Cartesian, Bohr.

- reference: exact periodic translation ``w(\cdot-R) \leftrightarrow
  w_G e^{-iG\cdot R}``, then grid dot products (no new machinery);
- Gaussian: `translate(Wc, R_cart)` then native `integral`.

The hoppings of a real Wannier are real; the imaginary parts of the reference
values are numerical leakage and asserted small, not returned.

Beware small denominators on distant hoppings: `rel_err_*` divides by the
reference value, `abs_err_*_per_T0` divides by the on-site kinetic ``T(0)`` —
report both (design note 09).
"""
function compare_hopping(Wc::CompressedWannier, w_fourier::AbstractVector,
                         basis_supercell::PlaneWaveBasis, R_cart::AbstractVector)
    @assert length(R_cart) == 3 "R_cart must be a 3-vector, got $(length(R_cart))"
    @assert all(isfinite, w_fourier) "non-finite coefficients in the reference Wannier"
    kpt = only(basis_supercell.kpoints)
    Gs = G_vectors_cart(basis_supercell, kpt)
    @assert length(w_fourier) == length(Gs) "Wannier/basis mismatch"

    phase = [cis(-dot(G, R_cart)) for G in Gs]
    wR = w_fourier .* phase
    G² = [sum(abs2, G) for G in Gs]

    S_ref_c = dot(w_fourier, wR)                 # ⟨w₀, w_R⟩
    T_ref_c = 0.5 * dot(w_fourier, G² .* wR)     # ⟨w₀, -½∇² w_R⟩
    norm²_ref = real(dot(w_fourier, w_fourier))
    T0_ref = 0.5 * real(dot(w_fourier, G² .* w_fourier))  # on-site scale

    IMAG_TOL = 1e-6  # real-Wannier sanity: imaginary leakage bound (relative)
    @assert abs(imag(S_ref_c)) ≤ IMAG_TOL * norm²_ref "Im S(R) leakage: $(imag(S_ref_c))"
    @assert abs(imag(T_ref_c)) ≤ IMAG_TOL * abs(T0_ref) "Im T(R) leakage: $(imag(T_ref_c))"
    S_ref, T_ref = real(S_ref_c), real(T_ref_c)

    WcR = translate(Wc, R_cart)
    S_gto = integral(Wc, WcR; type=:overlap)
    T_gto = integral(Wc, WcR; type=:kinetic)
    norm²_gto = real(integral(Wc, Wc; type=:overlap))

    # Norm-corrected values (translation is unitary: ‖w_R‖ = ‖w₀‖).
    S_ref_n, T_ref_n = S_ref / norm²_ref, T_ref / norm²_ref
    S_gto_n, T_gto_n = S_gto / norm²_gto, T_gto / norm²_gto

    (; S_ref = S_ref_n, S_gto = S_gto_n,
       T_ref = T_ref_n, T_gto = T_gto_n,
       T0_ref = T0_ref / norm²_ref,
       rel_err_S = abs(S_gto_n - S_ref_n) / abs(S_ref_n),
       rel_err_T = abs(T_gto_n - T_ref_n) / abs(T_ref_n),
       abs_err_S_per_T0 = abs(S_gto_n - S_ref_n) / abs(T0_ref / norm²_ref),
       abs_err_T_per_T0 = abs(T_gto_n - T_ref_n) / abs(T0_ref / norm²_ref),
       norm_ref = sqrt(norm²_ref), norm_gto = sqrt(norm²_gto))
end
```

Attention cohérence avec les tests du Step 1 : les tests GTO comparent
`out.S_gto` aux formes closes de gaussiennes **normalisées** — comme
`norm²_gto = 1` pour la gaussienne synthétique normalisée, les valeurs
norm-corrected coïncident avec les valeurs brutes et les tests passent tels
quels. Même chose côté référence (`w` normalisé ⇒ `norm²_ref = 1`).

- [ ] **Step 4: Export**

Dans `src/Wannier2GTO.jl`, à côté de l'export existant ligne 67 :

```julia
export compare_onsite_kinetic
export compare_hopping
```

- [ ] **Step 5: Run tests, verify they pass**

Run: `julia --project=. test/test_hoppings.jl`
Expected: PASS (deux testsets verts). Si le testset "reference side" casse sur
la signature `DFTK.fft`, appliquer la note du Step 1 et relancer.

- [ ] **Step 6: Vérifier que l'existant ne casse pas**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: PASS.

- [ ] **Step 7: Commit (proposer)**

```bash
git add test/test_hoppings.jl src/tight_binding/hamiltonian_scalar_product.jl src/Wannier2GTO.jl
git commit -m "TB validation: compare_hopping (S,T at R≠0), closed-form and plane-wave identity tests"
```

---

### Task 2: script workflow `validate_tb_hoppings.jl`

**Files:**
- Create: `workflow/validate_tb_hoppings.jl`

- [ ] **Step 1: Write the script**

Même squelette que [validate_tb_onsite.jl](../workflow/validate_tb_onsite.jl)
(mêmes fichiers de données, mêmes gardes), en gardant le basis unit-cell pour
extraire a₁, a₂ :

```julia
#
# Hopping validation of the tight-binding chain (approach A, design note 09).
#
# Compares S(R) = ⟨w₀|w_R⟩ and T(R) = ⟨w₀|−½∇²|w_R⟩ from the true plane-wave
# Wannier (Fourier phase translation + grid) against the Gaussian-compressed
# counterpart (rigid translation + analytic integrals).
#
# R sampling (design note 09): intralayer a₁, a₁+a₂, 2a₁; pseudo-interlayer
# z-shift of 3.35 Å with lateral offset δ ∈ {0, (a₁+a₂)/3, a₁/2}. R = 0 is
# included as a consistency check against note 08.
#
# Usage:
#   julia --project=. workflow/validate_tb_hoppings.jl          # Ecut = 50
#   W2G_ECUT=15 julia --project=. workflow/validate_tb_hoppings.jl   # local smoke
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using LinearAlgebra
using Printf

const ECUT = parse(Int, get(ENV, "W2G_ECUT", "50"))
const D_SLAB = 10.0        # Å, interlayer vacuum used for the wannierization
const KGRID = [5, 5, 1]
const D_INTERLAYER = austrip(3.35u"Å")   # Bernal-like interlayer distance

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
include(joinpath(ROOT, "workflow", "monolayer_graphene.jl"))

compressed_file(Ecut) = joinpath(ROOT, "data",
    Ecut == 15 ? "compressed_wannier_55_H1.json" :
                 "compressed_wannier_55_H1_Ecut$(Ecut).json")
reference_file(Ecut) = joinpath(ROOT, "workflow", "wannier_functions",
                                "wannier_pz_Ecut-$(Ecut).json")

for f in (compressed_file(ECUT), reference_file(ECUT))
    isfile(f) || error("missing data file: $f")
end

@info "Building bases" ECUT D_SLAB KGRID
basis_uc = Graphene(; d=D_SLAB * u"Å", kgrid=KGRID, Ecut=ECUT).basis()
basis_sc = DFTK.cell_to_supercell(basis_uc)

# Unit-cell lattice vectors, Cartesian Bohr (DFTK stores them as columns).
a₁ = basis_uc.model.lattice[:, 1]
a₂ = basis_uc.model.lattice[:, 2]
ẑd = [0.0, 0.0, D_INTERLAYER]

hoppings = [
    ("R=0 (check)",        zeros(3)),
    ("intra a1",           a₁),
    ("intra a1+a2",        a₁ + a₂),
    ("intra 2a1",          2a₁),
    ("inter AA (δ=0)",     ẑd),
    ("inter AB-like",      ẑd + (a₁ + a₂) / 3),
    ("inter mid (δ=a1/2)", ẑd + a₁ / 2),
]

ref = W2G.read_wannier_function(reference_file(ECUT))
w_fourier = ref.wannier
n_G = length(DFTK.G_vectors(basis_sc, only(basis_sc.kpoints)))
@assert length(w_fourier) == n_G "Wannier/basis mismatch: len(w)=$(length(w_fourier)) vs n_G=$n_G — check (d, Ecut, kgrid)"
@assert norm(w_fourier) ≈ 1 "reference Wannier is not normalized (‖w‖ = $(norm(w_fourier)))"

Wc = W2G.CompressedWannier(compressed_file(ECUT))
@assert !isempty(Wc.basis_functions) "compressed Wannier carries no basis function"

println()
println("Hoppings S(R), T(R)   (Ecut = $ECUT, d = $D_SLAB Å, kgrid = $KGRID, norm-corrected)")
println("  n basis functions (GTO) = ", length(Wc.basis_functions))
println("  compression error (H¹)  = ", Wc.error)
println()
@printf("  %-20s %8s │ %12s %12s %9s │ %12s %12s %9s %9s\n",
        "R", "|R| (Å)", "S_ref", "S_gto", "relerr_S",
        "T_ref", "T_gto", "relerr_T", "err/T(0)")
for (label, R) in hoppings
    out = W2G.compare_hopping(Wc, w_fourier, basis_sc, R)
    @printf("  %-20s %8.3f │ %12.3e %12.3e %9.2e │ %12.3e %12.3e %9.2e %9.2e\n",
            label, ustrip(auconvert(u"Å", norm(R))),
            out.S_ref, out.S_gto, out.rel_err_S,
            out.T_ref, out.T_gto, out.rel_err_T, out.abs_err_T_per_T0)
end
println()
println("  T(0) reference (norm-corrected) = ",
        @sprintf("%.10f", W2G.compare_hopping(Wc, w_fourier, basis_sc, zeros(3)).T0_ref))
println("  On-site reference values (note 08): rel err T(0) = 0.41 % (Ecut 15), 1.89 % (Ecut 50)")
println("  H¹ compression errors           : 5.59 % (Ecut 15), 12.38 % (Ecut 50)")
println()
```

Si `auconvert`/`ustrip` posent problème avec les imports du projet, remplacer
la colonne |R| par `norm(R) / austrip(1u"Å")` — cosmétique, ne pas y perdre du
temps.

- [ ] **Step 2: Sanity de la ligne R=0**

Run: `W2G_ECUT=15 julia --project=. workflow/validate_tb_hoppings.jl`
Expected: la ligne "R=0 (check)" affiche relerr_T ≈ 4.1e-3 (valeur note 08,
Ecut 15) et S_gto ≈ 1 à l'erreur de norme près. Toute autre valeur = bug de
câblage, **stop et debug avant d'interpréter le reste**.

- [ ] **Step 3: Commit (proposer)**

```bash
git add workflow/validate_tb_hoppings.jl
git commit -m "TB validation: hopping workflow script (approach A, design note 09)"
```

---

### Task 3: runs et bilan

- [ ] **Step 1: Run Ecut = 15 (local)**

Run: `W2G_ECUT=15 julia --project=. workflow/validate_tb_hoppings.jl`
Sauvegarder la sortie (rediriger vers
`workflow/diag_outputs/hoppings_Ecut15.log`).

- [ ] **Step 2: Run Ecut = 50**

Run: `julia --project=. workflow/validate_tb_hoppings.jl > workflow/diag_outputs/hoppings_Ecut50.log`
Le coût dominant est la construction du basis supercellule 5×5 à Ecut 50 et le
`n_G`-produit scalaire — même gabarit que le run on-site Ecut 50. Si la RAM
locale ne suffit pas, même procédure cluster que pour
`validate_tb_onsite.jl` (note 07).

- [ ] **Step 3: Bilan**

Écrire `notes/11-bilan-validation-hoppings.md` : tableau erreurs vs distance,
confrontation aux deux issues du design (§Critères de succès de la note 09) :
erreur hopping ~H¹ (premier ordre, prédiction confirmée → réorientation loss
phase B) ou ~on-site (superconvergence robuste → chaîne TB inchangée).
Mentionner le plancher d'images périodiques (préambule de ce plan) dans
l'interprétation des hoppings les plus lointains (2a₁, |R| ≈ 5.3 Å vs
supercellule 13.2 Å).

- [ ] **Step 4: Commit (proposer)**

```bash
git add workflow/diag_outputs/hoppings_Ecut*.log notes/11-bilan-validation-hoppings.md
git commit -m "TB validation: hopping results (Ecut 15/50) and bilan"
```
