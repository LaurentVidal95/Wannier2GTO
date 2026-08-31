# Plan — loss hopping-ciblée (implémentation)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implémenter la loss du design [12](12-design-loss-hoppings.md)
(H¹ + μ·pénalités {S(R),T(R)} + ν·orthonormalité) dans le joint-optim de la
phase B, avec cibles précalculées, K-restarts, évaluation par critère, et
paires contrôlées μ=0/μ>0.

**Architecture:** (1) un module `hopping_targets` (jeux de R, références
plane-wave, JSON, critères) ; (2) extension de `joint_loss` avec pénalité
hoppings analytique différentiable ; (3) `run_joint_optim` (Task 10 phase B,
note 04) adapté avec passthrough μ/ν ; (4) scripts workflow précalcul + scan μ
apparié. Init greedy non touchée (véhicule (c) : le raffinement part
d'`init_params` comme prévu phase B ; la comparaison contrôlée est le jumeau
μ=0).

**Tech Stack:** Julia, DFTK, Optim/L-BFGS + ForwardDiff (backend validé Task
8v2), intégrales gaussiennes natives, JSON3.

**Écart assumé au design 12 §3 (init)** : le design évoquait une init par la
sortie greedy phase A ; or la machinerie phase B initialise par `init_params`
heuristique + K restarts (note 04, Task 5), et encoder les SAGTOs greedy dans
la paramétrisation contrainte `JointLayout` est un composant qui n'existe pas.
v1 = K-restarts `init_params` (les paires μ=0/μ>0 restent contrôlées — mêmes
seeds). Si le smoke (Task 8) montre que les restarts n'atteignent pas le
niveau H¹ du greedy (~5.6 % à Ecut 15) même à μ=0, l'encodage greedy→layout
devient un préalable et fera l'objet d'un mini-design dédié.

**Conventions repo** : tests fichier par fichier, commits **proposés** à
l'utilisateur (règle research-workflow) — chaque étape "Commit" = proposer.
Julia startup ~minutes : timeouts Bash généreux (600000 ms).

**Faits vérifiés sur le code** (ne pas redécouvrir) :
- `joint_loss` (`src/joint_optim/loss.jl:41`) : variable projection — Φs via
  `params_to_basis_functions`, coefficients c via solve Tikhonov, retour ‖résidu‖²_H¹.
- `run_lbfgs_once` (`src/joint_optim/runner.jl:31`) : L-BFGS Optim +
  `ForwardDiff.gradient` sur closure de `joint_loss`.
- `run_joint_optim` n'existe **pas encore** (Task 10 phase B en pause, code
  prêt dans `notes/04-plan-phase-B.md:952-1074` — le reprendre).
- `translate` n'existe que pour `CompressedWannier`
  (`src/CompressedWannier.jl:33`), pas pour `BasisFunction`.
- `integral(Φ₁::BasisFunction, Φ₂; type=:overlap/:kinetic)` natif
  (`src/integrals/integrals.jl:51`).
- `compare_hopping` (`src/tight_binding/hamiltonian_scalar_product.jl`, testé
  par `test/test_hoppings.jl`, 13 tests verts) contient déjà le calcul de
  référence par phases — Task 1 le factorise.
- Helpers de test réutilisables dans `test/test_hoppings.jl` :
  `single_gaussian_Wc(ζ)`, `tiny_basis()` (recopier, ne pas importer).
- Échelles (bilan 11, Ecut 15) : T(0)≈1.47 Ha, ΔT(a₁)≈7.5e-3 Ha,
  terme H¹ ≈ (5.6e-2)² ≈ 3.1e-3 ; terme hopping ω-pondéré ≈ 1e-4 → μ utile
  ~1–1e3.

---

### Task 1: `translate(::BasisFunction)` + `reference_hopping` (refactor DRY)

**Files:**
- Modify: `src/CompressedWannier.jl:33-41`
- Modify: `src/tight_binding/hamiltonian_scalar_product.jl` (compare_hopping)
- Modify: `src/Wannier2GTO.jl:67` (exports)
- Test: `test/test_hoppings.jl` (append)

- [ ] **Step 1: Write the failing tests** (append à `test/test_hoppings.jl`)

```julia
@testset "translate(BasisFunction) and reference_hopping" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φ = W2G.BasisFunction([1.0], [g])
    R = [0.4, -0.2, 0.9]
    ΦR = W2G.translate(Φ, R)
    # Same closed form as the CompressedWannier path.
    @test W2G.integral(Φ, ΦR; type=:overlap) ≈ exp(-ζ * sum(abs2, R) / 2)  rtol = 1e-10

    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w_real = randn(basis.fft_size...)
    w = DFTK.fft(basis, kpt, complex.(w_real))
    w = w / norm(w)
    out0 = W2G.reference_hopping(w, basis, zeros(3))
    @test out0.S_ref ≈ 1.0  rtol = 1e-12
    @test out0.T_ref ≈ out0.T0_ref  rtol = 1e-12
    # Must agree with compare_hopping's reference side.
    Wc = single_gaussian_Wc(ζ)
    Rb = [1.1, -0.4, 0.6]
    @test W2G.reference_hopping(w, basis, Rb).S_ref ≈ W2G.compare_hopping(Wc, w, basis, Rb).S_ref  rtol = 1e-12
    @test W2G.reference_hopping(w, basis, Rb).T_ref ≈ W2G.compare_hopping(Wc, w, basis, Rb).T_ref  rtol = 1e-12
end
```

- [ ] **Step 2: Run, verify FAIL**

Run: `julia --project=. test/test_hoppings.jl`
Expected: FAIL — `translate` sans méthode `BasisFunction` /
`reference_hopping not defined`. Les 13 tests existants restent verts.

- [ ] **Step 3: Implement**

Dans `src/CompressedWannier.jl`, remplacer `translate(Wc, R)` par :

```julia
function translate(Φ::BasisFunction, R::AbstractVector)
    BasisFunction(Φ.coeffs,
                  [GaussianPolynomial(X, X.center + R) for X in Φ.SAGTOs])
end

function translate(Wc::CompressedWannier, R::AbstractVector{T}) where T
    basis_functions = Vector{BasisFunction}(map(Φ -> translate(Φ, R), Wc.basis_functions))
    CompressedWannier(Wc.center + R, basis_functions, Wc.coefficients,
                      Wc.basis_supercell, Wc.wannier, Wc.residual, Wc.error, Wc.error_norm)
end
```

Dans `src/tight_binding/hamiltonian_scalar_product.jl`, ajouter avant
`compare_hopping` :

```julia
@doc raw"""
Reference (plane-wave) hopping values, norm-corrected: ``S(R), T(R)`` of the
periodic Wannier under exact phase translation ``w_G e^{-iG\cdot R}``, plus the
on-site kinetic scale ``T(0)``. Extracted from `compare_hopping` so target
precomputation does not need a Gaussian side.
"""
function reference_hopping(w_fourier::AbstractVector,
                           basis_supercell::PlaneWaveBasis, R_cart::AbstractVector)
    @assert length(R_cart) == 3 "R_cart must be a 3-vector, got $(length(R_cart))"
    @assert all(isfinite, w_fourier) "non-finite coefficients in the reference Wannier"
    kpt = only(basis_supercell.kpoints)
    Gs = G_vectors_cart(basis_supercell, kpt)
    @assert length(w_fourier) == length(Gs) "Wannier/basis mismatch"

    phase = [cis(-dot(G, R_cart)) for G in Gs]
    wR = w_fourier .* phase
    G² = [sum(abs2, G) for G in Gs]

    S_c = dot(w_fourier, wR)
    T_c = 0.5 * dot(w_fourier, G² .* wR)
    norm² = real(dot(w_fourier, w_fourier))
    T0 = 0.5 * real(dot(w_fourier, G² .* w_fourier)) / norm²

    IMAG_TOL = 1e-6  # real-Wannier sanity: imaginary leakage bound (relative)
    @assert abs(imag(S_c)) ≤ IMAG_TOL * norm² "Im S(R) leakage: $(imag(S_c))"
    @assert abs(imag(T_c)) ≤ IMAG_TOL * abs(T0) * norm² "Im T(R) leakage: $(imag(T_c))"

    (; S_ref = real(S_c) / norm², T_ref = real(T_c) / norm², T0_ref = T0)
end
```

Puis dans `compare_hopping`, remplacer tout le bloc référence (du calcul de
`kpt`/`Gs`/`phase` jusqu'à `S_ref, T_ref = ...` et les normes réf) par :

```julia
    ref = reference_hopping(w_fourier, basis_supercell, R_cart)
    norm²_ref = real(dot(w_fourier, w_fourier))
```

et utiliser ensuite `ref.S_ref`, `ref.T_ref`, `ref.T0_ref` (déjà
norm-corrected — supprimer les divisions par `norm²_ref` côté référence ;
conserver `T0_ref = ref.T0_ref` dans le tuple de retour, et
`norm_ref = sqrt(norm²_ref)`). Les erreurs `rel_err_*` / `abs_err_*_per_T0`
gardent les mêmes formules avec ces valeurs.

Dans `src/Wannier2GTO.jl`, sous `export compare_hopping` :

```julia
export reference_hopping
```

- [ ] **Step 4: Run tests, verify PASS**

Run: `julia --project=. test/test_hoppings.jl`
Expected: PASS, 13 anciens + nouveaux tests verts (le refactor de
`compare_hopping` est protégé par les 13 anciens).

- [ ] **Step 5: Commit (proposer)**

```bash
git add src/CompressedWannier.jl src/tight_binding/hamiltonian_scalar_product.jl src/Wannier2GTO.jl test/test_hoppings.jl
git commit -m "Hopping loss: translate(BasisFunction), reference_hopping factored out of compare_hopping"
```

---

### Task 2: `HoppingTargets` — jeux de R, références, JSON, critères

**Files:**
- Create: `src/tight_binding/hopping_targets.jl`
- Modify: `src/Wannier2GTO.jl` (include + exports)
- Test: `test/test_hopping_targets.jl` (create)

- [ ] **Step 1: Write the failing tests**

`test/test_hopping_targets.jl` :

```julia
using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra

function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "hopping_R_sets: determinism and geometry" begin
    a₁ = [4.99, 0.0, 0.0]; a₂ = [-2.49, 4.32, 0.0]
    kw = (d_inter=6.33, z_val=(5.67, 6.99), seed=1234, n_random=5,
          r_min=3.78, r_max=7.56)   # Bohr ≈ (3.35, 3.0, 3.7, 2.0, 4.0) Å
    rs1 = W2G.hopping_R_sets(a₁, a₂; kw...)
    rs2 = W2G.hopping_R_sets(a₁, a₂; kw...)
    @test rs1.labels == rs2.labels
    @test all(rs1.Rs .≈ rs2.Rs)                      # same seed → identical
    @test count(==(:training), rs1.sets) == 5
    @test count(==(:validation), rs1.sets) == 4 + 5   # deterministic core + random
    @test rs1.ortho_idx == findall(l -> startswith(l, "intra"), rs1.labels)
    for i in eachindex(rs1.labels)
        startswith(rs1.labels[i], "rand") || continue
        @test kw.r_min ≤ norm(rs1.Rs[i]) ≤ kw.r_max
    end
    rs3 = W2G.hopping_R_sets(a₁, a₂; kw..., seed=4321)
    @test !all(rs1.Rs .≈ rs3.Rs)                     # seed changes randoms
end

@testset "build_hopping_targets + JSON round-trip" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...)))
    w = w / norm(w)
    a₁ = basis.model.lattice[:, 1]; a₂ = basis.model.lattice[:, 2]
    rsets = W2G.hopping_R_sets(a₁, a₂; d_inter=6.33, z_val=(5.67, 6.99),
                               seed=7, n_random=2, r_min=3.78, r_max=7.56)
    targets = W2G.build_hopping_targets(w, basis, rsets)
    @test targets.T0_ref ≈ W2G.reference_hopping(w, basis, zeros(3)).T0_ref
    i1 = findfirst(==("intra a1"), targets.labels)
    @test targets.S_ref[i1] ≈ W2G.reference_hopping(w, basis, a₁).S_ref  rtol = 1e-12

    file = joinpath(mktempdir(), "targets.json")
    W2G.store(targets; file)
    loaded = W2G.HoppingTargets(file)
    @test loaded.labels == targets.labels
    @test loaded.sets == targets.sets
    @test all(loaded.Rs .≈ targets.Rs)
    @test loaded.S_ref ≈ targets.S_ref
    @test loaded.T_ref ≈ targets.T_ref
    @test loaded.T0_ref ≈ targets.T0_ref
    @test loaded.ortho_idx == targets.ortho_idx
end
```

- [ ] **Step 2: Run, verify FAIL**

Run: `julia --project=. test/test_hopping_targets.jl`
Expected: FAIL — `hopping_R_sets not defined`.

- [ ] **Step 3: Implement `src/tight_binding/hopping_targets.jl`**

```julia
using Random

@doc raw"""
Hopping-loss targets (design note 12 §4): labeled displacement set split into
training/validation, with reference values ``S_{\mathrm{ref}}(R),
T_{\mathrm{ref}}(R)`` (norm-corrected) and the on-site scale ``T(0)``.
`ortho_idx` indexes the training intralayer entries whose exact reference is
``S = 0`` (MLWF orthonormality) — they receive the extra ν penalty.
All lengths in Bohr, all Rs Cartesian.
"""
struct HoppingTargets{T<:Real}
    labels    :: Vector{String}
    sets      :: Vector{Symbol}       # :training | :validation
    Rs        :: Vector{Vector{T}}
    S_ref     :: Vector{T}
    T_ref     :: Vector{T}
    T0_ref    :: T
    ortho_idx :: Vector{Int}
end

@doc raw"""
Displacement sets of design note 12 §4. Training: note-09 set minus 2a₁
(periodic-image contamination). Validation: deterministic core (δ=a₁/4 shift,
two alternative z distances, a₂ for the C₃ check) plus `n_random` seeded
random draws in the overlap shell `r_min ≤ |R| ≤ r_max`.
`d_inter` and `z_val` are z distances in Bohr; the caller converts from Å.
"""
function hopping_R_sets(a₁::AbstractVector, a₂::AbstractVector;
                        d_inter::Real, z_val::Tuple, seed::Int,
                        n_random::Int, r_min::Real, r_max::Real)
    @assert 0 < r_min < r_max "invalid shell: [$r_min, $r_max]"
    ẑ(d) = [0.0, 0.0, float(d)]
    entries = Tuple{String, Vector{Float64}, Symbol}[
        ("intra a1",        Vector{Float64}(a₁),                :training),
        ("intra a1+a2",     Vector{Float64}(a₁ + a₂),           :training),
        ("inter AA",        ẑ(d_inter),                          :training),
        ("inter AB-like",   ẑ(d_inter) + (a₁ + a₂) / 3,          :training),
        ("inter mid",       ẑ(d_inter) + a₁ / 2,                 :training),
        ("val inter d=a1/4", ẑ(d_inter) + a₁ / 4,                :validation),
        ("val inter z1",     ẑ(z_val[1]),                        :validation),
        ("val inter z2",     ẑ(z_val[2]),                        :validation),
        ("val intra a2",     Vector{Float64}(a₂),                :validation),
    ]
    rng = MersenneTwister(seed)
    for k in 1:n_random
        u = randn(rng, 3); u /= norm(u)
        r = r_min + (r_max - r_min) * rand(rng)
        push!(entries, ("rand $k", r * u, :validation))
    end
    labels = [e[1] for e in entries]
    Rs     = [e[2] for e in entries]
    sets   = [e[3] for e in entries]
    ortho_idx = findall(l -> startswith(l, "intra"), labels)
    (; labels, Rs, sets, ortho_idx)
end

function build_hopping_targets(w_fourier::AbstractVector,
                               basis_supercell::PlaneWaveBasis, rsets)
    refs = [reference_hopping(w_fourier, basis_supercell, R) for R in rsets.Rs]
    T0 = reference_hopping(w_fourier, basis_supercell, zeros(3)).T0_ref
    HoppingTargets(rsets.labels, rsets.sets, rsets.Rs,
                   [r.S_ref for r in refs], [r.T_ref for r in refs],
                   T0, rsets.ortho_idx)
end

function store(targets::HoppingTargets; file::String)
    data = Dict("labels" => targets.labels,
                "sets" => String.(targets.sets),
                "Rs" => targets.Rs,
                "S_ref" => targets.S_ref, "T_ref" => targets.T_ref,
                "T0_ref" => targets.T0_ref, "ortho_idx" => targets.ortho_idx)
    open(io -> JSON3.write(io, data), file, "w")
    nothing
end

function HoppingTargets(file::String)
    d = open(JSON3.read, file)
    HoppingTargets(String.(d.labels), Symbol.(d.sets),
                   [Vector{Float64}(R) for R in d.Rs],
                   Vector{Float64}(d.S_ref), Vector{Float64}(d.T_ref),
                   Float64(d.T0_ref), Vector{Int}(d.ortho_idx))
end
```

Dans `src/Wannier2GTO.jl` : ajouter
`include("tight_binding/hopping_targets.jl")` **après**
`include("tight_binding/hamiltonian_scalar_product.jl")` (dépend de
`reference_hopping`), et exporter :

```julia
export HoppingTargets, hopping_R_sets, build_hopping_targets
```

(`store` est déjà défini pour `CompressedWannier` dans le module — la nouvelle
méthode s'ajoute par dispatch, ne pas ré-exporter.)

- [ ] **Step 4: Run tests, verify PASS**

Run: `julia --project=. test/test_hopping_targets.jl` — Expected: PASS.
Run: `julia --project=. test/test_hoppings.jl` — Expected: toujours PASS.

- [ ] **Step 5: Commit (proposer)**

```bash
git add src/tight_binding/hopping_targets.jl src/Wannier2GTO.jl test/test_hopping_targets.jl
git commit -m "Hopping loss: HoppingTargets (R sets, references, JSON round-trip)"
```

---

### Task 3: script de précalcul des cibles + données Ecut 15/50

**Files:**
- Create: `workflow/precompute_hopping_targets.jl`
- Create (données): `data/hopping_targets_Ecut15.json`, `data/hopping_targets_Ecut50.json`

- [ ] **Step 1: Write the script**

```julia
#
# Precompute hopping-loss targets (design note 12 §4): reference S(R), T(R)
# from the true plane-wave Wannier, for the training/validation R sets.
# Stored once in data/, so the optimization loop never touches plane waves
# for the hopping terms.
#
# Usage:
#   julia --project=. workflow/precompute_hopping_targets.jl          # Ecut = 50
#   W2G_ECUT=15 julia --project=. workflow/precompute_hopping_targets.jl
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using LinearAlgebra
using Printf

const ECUT = parse(Int, get(ENV, "W2G_ECUT", "50"))
const D_SLAB = 10.0      # Å, interlayer vacuum used for the wannierization
const KGRID = [5, 5, 1]
const SEED = 1234        # validation random-shell draws (reproducibility)
const N_RANDOM = 5

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
include(joinpath(ROOT, "workflow", "monolayer_graphene.jl"))

reference_file(Ecut) = joinpath(ROOT, "workflow", "wannier_functions",
                                "wannier_pz_Ecut-$(Ecut).json")
target_file(Ecut) = joinpath(ROOT, "data", "hopping_targets_Ecut$(Ecut).json")
isfile(reference_file(ECUT)) || error("missing data file: $(reference_file(ECUT))")

@info "Building bases" ECUT D_SLAB KGRID
basis_uc = Graphene(; d=D_SLAB * u"Å", kgrid=KGRID, Ecut=ECUT).basis()
basis_sc = DFTK.cell_to_supercell(basis_uc)
a₁ = basis_uc.model.lattice[:, 1]
a₂ = basis_uc.model.lattice[:, 2]

ref = W2G.read_wannier_function(reference_file(ECUT))
w_fourier = ref.wannier
n_G = length(DFTK.G_vectors(basis_sc, only(basis_sc.kpoints)))
@assert length(w_fourier) == n_G "Wannier/basis mismatch: len(w)=$(length(w_fourier)) vs n_G=$n_G"
@assert norm(w_fourier) ≈ 1 "reference Wannier is not normalized"

rsets = W2G.hopping_R_sets(a₁, a₂;
    d_inter=austrip(3.35u"Å"), z_val=(austrip(3.0u"Å"), austrip(3.7u"Å")),
    seed=SEED, n_random=N_RANDOM,
    r_min=austrip(2.0u"Å"), r_max=austrip(4.0u"Å"))
targets = W2G.build_hopping_targets(w_fourier, basis_sc, rsets)
W2G.store(targets; file=target_file(ECUT))

println()
println("Hopping targets (Ecut = $ECUT, seed = $SEED)  →  $(target_file(ECUT))")
println("  T(0) = ", @sprintf("%.10f", targets.T0_ref))
@printf("  %-18s %-11s %8s │ %12s %12s\n", "label", "set", "|R| (Å)", "S_ref", "T_ref")
for i in eachindex(targets.labels)
    @printf("  %-18s %-11s %8.3f │ %12.4e %12.4e\n",
            targets.labels[i], targets.sets[i],
            norm(targets.Rs[i]) / austrip(1u"Å"),
            targets.S_ref[i], targets.T_ref[i])
end
```

- [ ] **Step 2: Run both Ecuts**

```bash
W2G_ECUT=15 julia --project=. workflow/precompute_hopping_targets.jl
```
```bash
julia --project=. workflow/precompute_hopping_targets.jl
```
Expected : tables affichées ; sanity : "intra a1" doit donner
S_ref ~1e-16 et T_ref ≈ 1.30e-2 (Ecut 15) / 3.89e-2 (Ecut 50), "inter AA"
S_ref ≈ −4.95e-2 / −5.52e-2 (valeurs bilan 11). Écart = bug, stop.

- [ ] **Step 3: Commit (proposer)** — les JSON de cibles sont de petites
données de référence versionnées (comme `data/compressed_wannier_*.json`) :

```bash
git add workflow/precompute_hopping_targets.jl data/hopping_targets_Ecut15.json data/hopping_targets_Ecut50.json
git commit -m "Hopping loss: precomputed S/T targets (Ecut 15/50, seeded validation shell)"
```

---

### Task 4: pénalité hoppings dans `joint_loss` (différentiable)

**Files:**
- Modify: `src/joint_optim/loss.jl`
- Modify: `src/Wannier2GTO.jl` (exports)
- Test: `test/test_loss_hoppings.jl` (create)

- [ ] **Step 1: Write the failing tests**

`test/test_loss_hoppings.jl` :

```julia
using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra
using ForwardDiff

function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "gto_hoppings: closed forms" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    for R in ([1.0, 0.0, 0.0], [0.5, -0.3, 0.2])
        R² = sum(abs2, R)
        h = W2G.gto_hoppings(Φs, c, R)
        @test h.S ≈ exp(-ζ * R² / 2)                     rtol = 1e-10
        @test h.T ≈ (ζ/2) * (3 - ζ*R²) * exp(-ζ*R²/2)   rtol = 1e-10
    end
end

# Synthetic targets around a single training R + one ortho index.
function synthetic_targets(; S1, T1, T0)
    W2G.HoppingTargets(["intra a1"], [:training], [[1.0, 0.0, 0.0]],
                       [S1], [T1], T0, [1])
end

@testset "hopping_penalty: zero at target, positive off target" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    h = W2G.gto_hoppings(Φs, c, [1.0, 0.0, 0.0])   # norm² = 1 here
    t_exact = synthetic_targets(S1=h.S, T1=h.T, T0=1.5)
    @test W2G.hopping_penalty(Φs, c, t_exact; μ=1.0, ν=0.0) ≈ 0.0  atol = 1e-14
    t_off = synthetic_targets(S1=h.S + 0.1, T1=h.T, T0=1.5)
    @test W2G.hopping_penalty(Φs, c, t_off; μ=1.0, ν=0.0) ≈ 0.1^2 / 1.5^2  rtol = 1e-10
    # ν penalty adds S² for ortho-indexed entries.
    @test W2G.hopping_penalty(Φs, c, t_exact; μ=0.0, ν=1.0) ≈ h.S^2  rtol = 1e-10
end

@testset "joint_loss: targets=nothing ≡ μ=0, gradient finite" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...)))
    w = w / norm(w)
    layout = W2G.JointLayout(N_centered=1, N_pibond=1)
    flat = W2G.init_params(W2G.MersenneTwister(11), layout)
    kw = (log_ζ_min=log(1e-2), log_ζ_max=log(4.0), π_bond_unit=[1.0, 0.0, 0.0])

    L0 = W2G.joint_loss(flat, layout, w, basis; kw...)
    targets = synthetic_targets(S1=0.02, T1=0.01, T0=1.5)
    @test W2G.joint_loss(flat, layout, w, basis; kw..., targets, μ=0.0, ν=0.0) ≈ L0  rtol = 1e-12
    Lμ = W2G.joint_loss(flat, layout, w, basis; kw..., targets, μ=10.0, ν=10.0)
    @test Lμ > L0

    f = θ -> W2G.joint_loss(θ, layout, w, basis; kw..., targets, μ=10.0, ν=10.0)
    G = ForwardDiff.gradient(f, flat)
    @test all(isfinite, G)
    G0 = ForwardDiff.gradient(θ -> W2G.joint_loss(θ, layout, w, basis; kw...), flat)
    @test !(G ≈ G0)   # the penalty actually contributes to the gradient
end
```

Note : si `W2G.MersenneTwister` n'est pas accessible via le module, faire
`using Random` et `MersenneTwister(11)` directement — détail d'import, pas de
design.

- [ ] **Step 2: Run, verify FAIL**

Run: `julia --project=. test/test_loss_hoppings.jl`
Expected: FAIL — `gto_hoppings not defined`.

- [ ] **Step 3: Implement** (append à `src/joint_optim/loss.jl`, et modifier
`joint_loss`)

```julia
@doc raw"""
Raw (NOT norm-corrected) Gaussian-side hoppings of the fitted function
``w_{\mathrm{GTO}} = \sum_i c_i \Phi_i``:
``S = \langle w, w(\cdot - R)\rangle``, ``T = \langle w, -\tfrac12\nabla^2
w(\cdot - R)\rangle``, via analytic integrals. ForwardDiff-generic in the
parameters carried by `Φs` and `c`.
"""
function gto_hoppings(Φs::AbstractVector{<:BasisFunction}, c::AbstractVector,
                      R::AbstractVector)
    ΦsR = [translate(Φ, R) for Φ in Φs]
    S = zero(eltype(c)); T = zero(eltype(c))
    for i in eachindex(Φs), j in eachindex(ΦsR)
        S += c[i] * c[j] * integral(Φs[i], ΦsR[j]; type=:overlap)
        T += c[i] * c[j] * integral(Φs[i], ΦsR[j]; type=:kinetic)
    end
    (; S, T)
end

@doc raw"""
Hopping penalty of design note 12 §2:
``\mu \sum_{R\in\mathcal{T}} [\Delta T(R)^2 + \Delta S(R)^2]/T(0)^2
 + \nu \sum_{R_{\mathrm{intra}}} S_{\mathrm{GTO}}(R)^2``,
training entries only, norm-corrected GTO values.
"""
function hopping_penalty(Φs::AbstractVector{<:BasisFunction}, c::AbstractVector,
                         targets::HoppingTargets; μ::Real, ν::Real)
    norm² = gto_hoppings(Φs, c, zeros(3)).S
    pen = zero(norm²); ortho = zero(norm²)
    for i in eachindex(targets.labels)
        targets.sets[i] == :training || continue
        h = gto_hoppings(Φs, c, targets.Rs[i])
        Sn, Tn = h.S / norm², h.T / norm²
        pen += ((Tn - targets.T_ref[i])^2 + (Sn - targets.S_ref[i])^2) / targets.T0_ref^2
        (i in targets.ortho_idx) && (ortho += Sn^2)
    end
    μ * pen + ν * ortho
end
```

Dans `joint_loss`, étendre la signature avec
`targets::Union{Nothing, HoppingTargets} = nothing, μ::Real = 0.0,
ν::Real = 0.0` et remplacer la dernière ligne
`Hˢ_norm(basis_supercell, residual; s=s)^2` par :

```julia
    L = Hˢ_norm(basis_supercell, residual; s=s)^2
    if !isnothing(targets) && (μ > 0 || ν > 0)
        L += hopping_penalty(Φs, c, targets; μ, ν)
    end
    L
```

Dans `src/Wannier2GTO.jl` : `export gto_hoppings, hopping_penalty`.

- [ ] **Step 4: Run tests, verify PASS**

Run: `julia --project=. test/test_loss_hoppings.jl` — Expected: PASS.
Si le test ForwardDiff échoue en `MethodError` sur des `Dual` dans
`integral`/`translate` : le point de friction attendu est une annotation de
type trop stricte (`Vector{Float64}`, `Float64`) dans `GaussianPolynomial`,
`translate` ou `overlap_julia`/`kinetic_julia` — relâcher vers des types
paramétriques `<:Real` au point exact de l'erreur, sans toucher aux
algorithmes, puis relancer aussi `test/test_integrals_laplacian.jl` et
`test/test_loss.jl`.
Run: `julia --project=. test/test_loss.jl` — Expected: PASS (non-régression
de la loss H¹, `targets=nothing` par défaut).

- [ ] **Step 5: Commit (proposer)**

```bash
git add src/joint_optim/loss.jl src/Wannier2GTO.jl test/test_loss_hoppings.jl
git commit -m "Hopping loss: gto_hoppings + hopping_penalty wired into joint_loss (ForwardDiff-checked)"
```

---

### Task 5: `run_joint_optim` (Task 10 phase B) avec passthrough μ/ν

**Files:**
- Modify: `src/joint_optim/runner.jl`

- [ ] **Step 1: Étendre `run_lbfgs_once`**

Ajouter à sa signature les kwargs
`targets::Union{Nothing, HoppingTargets} = nothing, μ::Real = 0.0,
ν::Real = 0.0`, et les passer dans la closure :

```julia
    f = θ -> joint_loss(θ, layout, w_z_fourier, basis_supercell;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit=π_bond_unit, s=s, ε=ε,
                        targets, μ, ν)
```

- [ ] **Step 2: Implémenter `run_joint_optim`**

Reprendre **intégralement** le code de `notes/04-plan-phase-B.md:962-1073`
(struct `JointOptimResult` + fonction), avec quatre modifications :

1. signature : ajouter `targets::Union{Nothing, HoppingTargets} = nothing,
   μ::Real = 0.0, ν::Real = 0.0` ;
2. l'appel `run_lbfgs_once(...)` : ajouter `targets, μ, ν` ;
3. `config.json` : ajouter au Dict `"μ" => μ, "ν" => ν,
   "targets_file" => "see caller"` (le chemin exact est loggé par le script
   appelant) ;
4. en tête de fichier : `using JSON3, Dates, Random` (MersenneTwister).

- [ ] **Step 3: Smoke K=2 minuscule**

```bash
julia --project=. -e '
using Random, LinearAlgebra, DFTK, Wannier2GTO; import Wannier2GTO as W2G
lattice = diagm([8.0, 9.0, 10.0])
model = Model(lattice; terms=[Kinetic()], n_electrons=1,
              spin_polarization=:spinless, symmetries=false)
basis = PlaneWaveBasis(model; Ecut=8, kgrid=(1,1,1))
kpt = only(basis.kpoints)
w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...))); w /= norm(w)
layout = W2G.JointLayout(N_centered=2, N_pibond=1)
targets = W2G.HoppingTargets(["intra a1"], [:training], [[1.0,0.0,0.0]],
                             [0.02], [0.01], 1.5, [1])
result = W2G.run_joint_optim(layout, w, basis;
    K=2, master_seed=100, log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
    π_bond_unit=[1.0,0.0,0.0], max_iter=5, targets, μ=10.0, ν=10.0)
println("OK  dir=", result.output_dir, "  best=", result.best_idx)
run(`ls $(result.output_dir)`)'
```

Expected: termine, dossier `workflow/joint_outputs/run_*/` avec `config.json`
(contenant μ=10), `restart_01/02.json`, `best.json`, `summary.txt`.

- [ ] **Step 4: Commit (proposer)**

```bash
git add src/joint_optim/runner.jl
git commit -m "Phase B Task 10 + hopping loss: K-restart run_joint_optim with mu/nu passthrough"
```

---

### Task 6: évaluation par critère

**Files:**
- Modify: `src/tight_binding/hopping_targets.jl` (append)
- Modify: `src/Wannier2GTO.jl` (exports)
- Test: `test/test_hopping_targets.jl` (append)

- [ ] **Step 1: Write the failing tests** (append)

```julia
@testset "evaluate_hoppings + hopping_criteria" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    R1 = [1.0, 0.0, 0.0]; R2 = [0.0, 0.0, 1.5]
    h1 = W2G.gto_hoppings(Φs, c, R1); h2 = W2G.gto_hoppings(Φs, c, R2)
    # Targets: entry 1 = training "intra a1" met exactly; entry 2 = validation
    # with a deliberately wrong reference sign.
    targets = W2G.HoppingTargets(["intra a1", "val x"], [:training, :validation],
                                 [R1, R2], [h1.S, h2.S], [h1.T, -h2.T],
                                 1.5, [1])
    rows = W2G.evaluate_hoppings(Φs, c, targets)
    @test length(rows) == 2
    @test rows[1].rel_err_T ≈ 0.0  atol = 1e-12
    @test rows[1].sign_ok
    @test !rows[2].sign_ok                      # wrong sign detected

    crit = W2G.hopping_criteria(rows, targets)
    @test length(crit) == 3
    c1 = crit[findfirst(x -> x.name == "relerr_T(a1) ≤ 5%", crit)]
    @test c1.pass && c1.value ≈ 0.0  atol = 1e-12
    c2 = crit[findfirst(x -> x.name == "validation signs", crit)]
    @test !c2.pass                              # |T_ref| = h2.T > floor, sign wrong
    c3 = crit[findfirst(x -> x.name == "|S(a1)| ≤ 1e-3", crit)]
    @test c3.value ≈ abs(h1.S)  rtol = 1e-10    # here S(a1) is large → FAIL
    @test !c3.pass
end
```

- [ ] **Step 2: Run, verify FAIL**, puis **Step 3: Implement** (append à
`src/tight_binding/hopping_targets.jl`) :

```julia
@doc raw"""
Norm-corrected GTO hoppings vs targets, one row per entry: values, errors
(relative and in units of ``T(0)``), sign check. Pure data — printing is the
caller's concern.
"""
function evaluate_hoppings(Φs::AbstractVector{<:BasisFunction}, c::AbstractVector,
                           targets::HoppingTargets)
    norm² = gto_hoppings(Φs, c, zeros(3)).S
    map(eachindex(targets.labels)) do i
        h = gto_hoppings(Φs, c, targets.Rs[i])
        Sn, Tn = h.S / norm², h.T / norm²
        (; label = targets.labels[i], set = targets.sets[i],
           S_ref = targets.S_ref[i], S_gto = Sn,
           T_ref = targets.T_ref[i], T_gto = Tn,
           rel_err_T = abs(Tn - targets.T_ref[i]) / abs(targets.T_ref[i]),
           abs_err_T_per_T0 = abs(Tn - targets.T_ref[i]) / abs(targets.T0_ref),
           sign_ok = sign(Tn) == sign(targets.T_ref[i]))
    end
end

# Design note 12 §5 thresholds. Criteria 1 and 3 are measured on "intra a1"
# (training: they test achievability of the loss's own target); criterion 2 on
# the validation set (generalization). Criterion 4 (H¹ ≤ 1.5× the μ=0 twin) is
# computed by the workflow script, which owns both runs of the pair.
function hopping_criteria(rows, targets::HoppingTargets;
                          tol_relT::Real = 0.05, tol_S::Real = 1e-3,
                          T_sign_floor::Real = 5e-4)
    i1 = findfirst(r -> r.label == "intra a1", rows)
    @assert !isnothing(i1) "targets must contain an 'intra a1' entry"
    val = [r for r in rows if r.set == :validation && abs(r.T_ref) > T_sign_floor]
    [(name = "relerr_T(a1) ≤ 5%", value = rows[i1].rel_err_T,
      pass = rows[i1].rel_err_T ≤ tol_relT),
     (name = "validation signs", value = count(r -> !r.sign_ok, val),
      pass = all(r -> r.sign_ok, val)),
     (name = "|S(a1)| ≤ 1e-3", value = abs(rows[i1].S_gto),
      pass = abs(rows[i1].S_gto) ≤ tol_S)]
end
```

Exports : `export evaluate_hoppings, hopping_criteria`.

- [ ] **Step 4: Run tests, verify PASS**

Run: `julia --project=. test/test_hopping_targets.jl` — Expected: PASS.

- [ ] **Step 5: Commit (proposer)**

```bash
git add src/tight_binding/hopping_targets.jl src/Wannier2GTO.jl test/test_hopping_targets.jl
git commit -m "Hopping loss: per-criterion evaluation (independent PASS/FAIL logging)"
```

---

### Task 7: script workflow scan μ apparié

**Files:**
- Create: `workflow/joint_compression_hoppings.jl`

- [ ] **Step 1: Write the script**

Basé sur le squelette de `notes/04-plan-phase-B.md:1113-1196`
(joint_compression.jl, jamais créé — ce script le remplace) :

```julia
#
# Joint compression with the hopping-targeted loss (design note 12).
#
# Runs a paired μ-scan: the μ = 0 twin (pure H¹) and each μ > 0 share the same
# master seed, hence identical per-restart inits — the controlled comparison
# required by design §5. After each run, rebuilds the best basis and logs the
# per-criterion evaluation (independent PASS/FAIL) on the validation targets.
#
# Usage:
#   W2G_ECUT=15 julia --project=. workflow/joint_compression_hoppings.jl   # smoke
#   julia --project=. workflow/joint_compression_hoppings.jl              # Ecut = 50
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using Random
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "monolayer_graphene.jl"))

# --- Configuration --------------------------------------------------------
const ECUT        = parse(Int, get(ENV, "W2G_ECUT", "50"))
const KGRID       = [5, 5, 1]
const D           = 10.0u"Å"
const N_CENTERED  = 5
const N_PIBOND    = 5
const K_RESTARTS  = 3
const MASTER_SEED = 42          # shared across the whole scan (paired runs)
const MAX_ITER    = parse(Int, get(ENV, "W2G_MAX_ITER", "100"))
const ζ_MIN       = 1e-2
const ζ_MAX       = ECUT / 4.0
const ε_TIKHONOV  = 1e-8
const MUS         = [0.0, 1e0, 1e1, 1e2, 1e3]   # μ scan; ν = μ (design §2)
const H1_REGRESSION_FACTOR = 1.5                 # criterion 4 threshold

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
wannier_file = joinpath(ROOT, "workflow", "wannier_functions",
                        "wannier_pz_Ecut-$(ECUT).json")
targets_file = joinpath(ROOT, "data", "hopping_targets_Ecut$(ECUT).json")
for f in (wannier_file, targets_file)
    isfile(f) || error("missing data file: $f — run precompute_hopping_targets.jl first")
end

# --- Setup ----------------------------------------------------------------
println("Building plane-wave basis (Ecut=$ECUT)...")
basis = Graphene(; d=D, kgrid=KGRID, Ecut=ECUT).basis()
basis_supercell = DFTK.cell_to_supercell(basis)

data = W2G.read_wannier_function(wannier_file)
w_z_fourier = normalize(data.wannier)
π_bond_unit = let θ = data.π_bond.θ; [cos(θ), sin(θ), 0.0] end
targets = W2G.HoppingTargets(targets_file)
layout = W2G.JointLayout(N_centered=N_CENTERED, N_pibond=N_PIBOND)
kw = (log_ζ_min=log(ζ_MIN), log_ζ_max=log(ζ_MAX), π_bond_unit=π_bond_unit)

rebuild(flat) = begin
    Φs = W2G.params_to_basis_functions(flat, layout; kw...)
    c = W2G.joint_inner_solve(Φs, w_z_fourier, basis_supercell; ε=ε_TIKHONOV)
    (Φs, c)
end

# --- Paired scan ----------------------------------------------------------
results = Dict{Float64, Any}()
for μ in MUS
    ν = μ
    println("\n", "=" ^ 70)
    println("RUN μ = $μ  (ν = $ν, master_seed = $MASTER_SEED)")
    println("=" ^ 70)
    res = W2G.run_joint_optim(layout, w_z_fourier, basis_supercell;
                              K=K_RESTARTS, master_seed=MASTER_SEED,
                              max_iter=MAX_ITER, ε=ε_TIKHONOV, kw...,
                              targets=(μ > 0 ? targets : nothing), μ, ν)
    best = res.restart_results[res.best_idx]
    Φs, c = rebuild(best.flat_final)
    # H¹ part alone (comparable across μ): recompute loss without penalty.
    h1² = W2G.joint_loss(best.flat_final, layout, w_z_fourier, basis_supercell;
                         kw..., ε=ε_TIKHONOV)
    rows = W2G.evaluate_hoppings(Φs, c, targets)
    crit = W2G.hopping_criteria(rows, targets)
    results[μ] = (; res, h1², rows, crit)

    @printf("H¹ relative error: %.4f%%   (loss with penalty: %.6e)\n",
            100 * sqrt(h1²), best.loss_final)
    @printf("%-18s %-11s │ %11s %11s │ %11s %11s %9s %9s %5s\n",
            "label", "set", "S_ref", "S_gto", "T_ref", "T_gto",
            "relerr_T", "err/T(0)", "sign")
    for r in rows
        @printf("%-18s %-11s │ %11.3e %11.3e │ %11.3e %11.3e %9.2e %9.2e %5s\n",
                r.label, r.set, r.S_ref, r.S_gto, r.T_ref, r.T_gto,
                r.rel_err_T, r.abs_err_T_per_T0, r.sign_ok ? "ok" : "BAD")
    end
    println("Criteria (independent):")
    for cr in crit
        @printf("  [%s] %-22s value = %.3e\n", cr.pass ? "PASS" : "FAIL",
                cr.name, float(cr.value))
    end
end

# --- Criterion 4: H¹ anti-regression vs the μ=0 twin ----------------------
println("\n", "=" ^ 70)
println("PAIRED SUMMARY (criterion 4: H¹ ≤ $(H1_REGRESSION_FACTOR)× the μ=0 twin)")
println("=" ^ 70)
h1_twin = results[0.0].h1²
@printf("%-10s %14s %10s %6s │ criteria 1-3\n", "μ", "H¹ rel err", "vs μ=0", "C4")
for μ in MUS
    r = results[μ]
    ratio = sqrt(r.h1² / h1_twin)
    c4 = ratio ≤ H1_REGRESSION_FACTOR
    c123 = join([cr.pass ? "PASS" : "FAIL" for cr in r.crit], " ")
    @printf("%-10.1e %13.4f%% %9.3fx %6s │ %s\n",
            μ, 100 * sqrt(r.h1²), ratio, c4 ? "PASS" : "FAIL", c123)
end
```

- [ ] **Step 2: Vérifier la cohérence apparillage**

Relire le script : `MASTER_SEED` identique pour tous les μ ⇒
`init_params(MersenneTwister(seed))` identiques par restart ⇒ comparaison
contrôlée. C'est **le** point du design §5 — ne pas "améliorer" en variant les
seeds.

- [ ] **Step 3: Commit (proposer)**

```bash
git add workflow/joint_compression_hoppings.jl
git commit -m "Hopping loss: paired mu-scan driver (shared seeds, per-criterion logging)"
```

---

### Task 8: smoke run Ecut 15 et mini-bilan

- [ ] **Step 1: Smoke réduit** (vérifie le câblage de bout en bout, pas les
critères — N et max_iter trop petits pour conclure) :

```bash
W2G_ECUT=15 W2G_MAX_ITER=30 julia --project=. -t auto workflow/joint_compression_hoppings.jl 2>&1 | tee workflow/diag_outputs/joint_hoppings_smoke15.log
```

Expected (~1-2 h laptop, dominé par ForwardDiff ; cold start ~10 min de
compilation) : le scan des 5 μ termine ; pour μ=0 la table hoppings ressemble
au bilan 11 (ordre de grandeur) ; quand μ croît, `relerr_T(a1)` et `|S(a1)|`
**décroissent de façon monotone ou quasi** tandis que l'erreur H¹ croît
doucement. Toute non-monotonie brutale ou NaN = bug (suspects : normalisation
norm², signe de la pénalité, Duals dans `integral`).

- [ ] **Step 2: Mini-bilan**

Écrire `notes/14-bilan-loss-hoppings-smoke.md` (court) : tableau
μ ↦ (H¹, critères 1-4), lecture du trade-off, et le go/no-go pour les runs
sérieux Ecut 50 (N=15, K=5, max_iter=200 — budgets note 04 §Task 12 :
~12 h laptop ou cluster). Les runs Ecut 50 et le bilan final (verdict design
12 §5, quel critère mord) sont **hors de ce plan** — session dédiée, décision
utilisateur sur le budget calcul.

- [ ] **Step 3: Commit (proposer)**

```bash
git add notes/14-bilan-loss-hoppings-smoke.md
git commit -m "Hopping loss: Ecut 15 smoke scan — mu trade-off and go/no-go for Ecut 50"
```
