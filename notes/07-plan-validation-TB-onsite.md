# Plan — validation on-site de la chaîne TB un corps

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute the native Julia kinetic/Laplacian integral for Gaussian-polynomial primitives and compare the on-site one-body integral ⟨w|−½∇²|w⟩ between the true plane-wave Wannier and the Gaussian-compressed Wannier.

**Architecture:** `laplacian_julia` expands ∇²X₂ (a Gaussian-polynomial) into a finite sum of monomials of the same spread/center, each integrated against X₁ via the existing exact `overlap_julia`; sign convention makes it the positive form ⟨∇X₁,∇X₂⟩ = −⟨X₁,∇²X₂⟩, with `kinetic_julia = ½ laplacian_julia`. Both are wired into `symb_to_integral`, replacing the dead GaIn calls. A thin harness then contrasts the analytic Gaussian kinetic energy against the grid kinetic energy ½⟨w,|G|²w⟩ of the reference Wannier.

**Tech Stack:** Julia 1.10+, DFTK (PlaneWaveBasis / G-vectors), ForwardDiff (AD-compatibility of the new integral), JSON3 (data loading).

**Spec:** `notes/06-design-validation-TB-onsite.md`. Scope this plan = L-1 + L-2a (kinetic). V_KS (L-2b) is a follow-up plan.

---

## File Structure

- `src/integrals/julia_integrals.jl` — **Modify.** Add `_axis_laplacian_terms`, `laplacian_julia`, `kinetic_julia` (same scalar signature as `overlap_julia`).
- `src/integrals/integrals.jl` — **Modify** (lines 9-16). Route `:laplacian` and `:kinetic` to the native functions instead of GaIn.
- `src/tight_binding/hamiltonian_scalar_product.jl` — **Modify.** Add `compare_onsite_kinetic(Wc, w_fourier, basis_supercell)`.
- `test/test_integrals_laplacian.jl` — **Create.** Closed-form, grid-coherence, symmetry, H¹-overlap, and harness integration tests.
- `workflow/validate_tb_onsite.jl` — **Create.** Loads the Ecut=50 data and prints the on-site kinetic comparison (run on cluster by the user).

Convention reminder (julia_integrals.jl): `X(r) = x^{nx} y^{ny} z^{nz} exp(−ζ‖r−R‖²)`, monomials in **absolute** coordinates, primitives passed unnormalized at the scalar level (normalization handled by `GaussianPolynomial`).

---

## Task 1: `laplacian_julia` / `kinetic_julia` + closed-form unit test

**Files:**
- Create: `test/test_integrals_laplacian.jl`
- Modify: `src/integrals/julia_integrals.jl` (append), `src/integrals/integrals.jl:9-16`

- [ ] **Step 1: Write the failing closed-form test**

Closed form: for a normalized s-Gaussian `g = exp(−ζr²)/‖·‖`, `⟨g,−∇²g⟩ = 3ζ`, so the positive Laplacian form equals `3ζ` and the kinetic equals `3ζ/2`.

Create `test/test_integrals_laplacian.jl`:

```julia
using Test
using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra

@testset "laplacian_julia: closed form (normalized s-Gaussian)" begin
    # g(r) = exp(-ζ r²), normalized to L² norm 1 by the GaussianPolynomial ctor.
    # ⟨∇g, ∇g⟩ = -⟨g, ∇²g⟩ = 3ζ ;  ⟨g, -½∇² g⟩ = 3ζ/2.
    for ζ in (0.3, 1.0, 2.5)
        g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)  # normalized
        @test W2G.integral(g, g; type=:laplacian) ≈ 3ζ      rtol = 1e-10
        @test W2G.integral(g, g; type=:kinetic)   ≈ 3ζ / 2  rtol = 1e-10
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: ERROR/FAIL — `symb_to_integral[:laplacian]` is `GaIn.laplacian`, whose backend is unavailable on this machine (dead hard-coded path).

- [ ] **Step 3: Implement the native Laplacian and wire it**

Append to `src/integrals/julia_integrals.jl`:

```julia
@doc raw"""
Per-axis monomial terms of ``\partial^2/\partial u^2`` acting on
``u^a\, e^{-\zeta (u-R)^2}``, written in the absolute coordinate ``u`` (the
Gaussian factor is unchanged). Returns `(coeff, new_exponent)` tuples:

```math
\partial_u^2 = a(a-1)\,u^{a-2}
             + \big(-2\zeta(1+2a) + 4\zeta^2 R^2\big)\,u^{a}
             + 4\zeta a R\, u^{a-1}
             - 8\zeta^2 R\, u^{a+1}
             + 4\zeta^2\, u^{a+2}.
```
"""
function _axis_laplacian_terms(a::Int, ζ::T, R) where {T}
    C = promote_type(T, typeof(R))
    terms = Tuple{C, Int}[]
    a ≥ 2 && push!(terms, (C(a * (a - 1)),                a - 2))
    push!(terms,           (-2ζ * (1 + 2a) + 4ζ^2 * R^2,   a    ))
    a ≥ 1 && push!(terms, (4ζ * a * R,                     a - 1))
    push!(terms,           (-8ζ^2 * R,                      a + 1))
    push!(terms,           (4ζ^2,                           a + 2))
    terms
end

@doc raw"""
Pure-Julia positive Laplacian form of two unnormalized Gaussian-polynomial
primitives, matching the GaIn convention used by `Hˢ_overlap(Ms; s=1)`:

```math
\texttt{laplacian\_julia}(X_1, X_2) = \langle \nabla X_1, \nabla X_2\rangle
    = -\langle X_1, \nabla^2 X_2\rangle .
```

``\nabla^2 X_2`` is a finite sum of monomials of the same spread/center as
``X_2``; each is integrated against ``X_1`` via the exact `overlap_julia`.
Same scalar signature as `overlap_julia` (drop-in for `symb_to_integral`),
ForwardDiff-compatible.
"""
function laplacian_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                         ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    n2 = (nx2, ny2, nz2)
    acc = zero(promote_type(typeof(ζ1), typeof(ζ2), eltype(R1), eltype(R2)))
    for axis in 1:3                      # ∇²X₂ = Σ_axis ∂²_axis X₂
        for (c, e) in _axis_laplacian_terms(n2[axis], ζ2, R2[axis])
            m = ntuple(i -> i == axis ? e : n2[i], 3)
            acc += c * overlap_julia(ζ1, R1, nx1, ny1, nz1,
                                     ζ2, R2, m[1], m[2], m[3])
        end
    end
    -acc    # ⟨∇X₁,∇X₂⟩ = -⟨X₁,∇²X₂⟩
end

"Kinetic-energy form ``\\langle X_1, -\\tfrac12\\nabla^2 X_2\\rangle = \\tfrac12\\,\\texttt{laplacian\\_julia}``."
function kinetic_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                       ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    laplacian_julia(ζ1, R1, nx1, ny1, nz1, ζ2, R2, nx2, ny2, nz2) / 2
end
```

Edit `src/integrals/integrals.jl` lines 9-16, replacing the GaIn routing for these two keys:

```julia
symb_to_integral = Dict([:overlap => overlap_julia,
                         :overlap_upper_bound => GaIn.overlap_upper_bound,
                         :laplacian => laplacian_julia,
                         :kinetic => kinetic_julia,
                         :coulomb => GaIn.coulomb,
                         :atomic => GaIn.atomic,
                         ]
                        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: PASS (closed-form testset green).

- [ ] **Step 5: Commit**

```bash
git add src/integrals/julia_integrals.jl src/integrals/integrals.jl test/test_integrals_laplacian.jl
git commit -m "TB integrals: native laplacian_julia/kinetic_julia (GaIn-free)

Analytic Laplacian of a Gaussian-polynomial via overlap_julia on the expanded
∇²X₂; positive form ⟨∇X₁,∇X₂⟩=-⟨X₁,∇²X₂⟩ (matches Hˢ_overlap), kinetic=½ that.
Closed-form ⟨g,-½∇²g⟩=3ζ/2 for normalized s-Gaussian verified.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: grid auto-validation, symmetry, and H¹-overlap tests

**Files:**
- Modify: `test/test_integrals_laplacian.jl` (append testsets + heavy fixture)

- [ ] **Step 1: Write the failing/anchored tests**

Append to `test/test_integrals_laplacian.jl`:

```julia
using DFTK
using DFTK.Unitful

# Heavy fixture: a Γ-point supercell plane-wave basis (built once).
include(joinpath(@__DIR__, "..", "workflow", "monolayer_graphene.jl"))
const _BASIS_SC = let
    basis = Graphene(; d=10u"Å", kgrid=[5, 5, 1], Ecut=15).basis()
    DFTK.cell_to_supercell(basis)
end

# ⟨∇X₁,∇X₂⟩ on the grid: ∇² ↦ -|G|² in Fourier, so the positive form is
# dot(X₁_four, |G|² .* X₂_four).
function _grid_laplacian(X1::W2G.GaussianPolynomial, X2::W2G.GaussianPolynomial)
    kpt = only(_BASIS_SC.kpoints)
    G2 = [sum(abs2, q) for q in DFTK.G_vectors_cart(_BASIS_SC, kpt)]
    real(dot(X1(_BASIS_SC), G2 .* X2(_BASIS_SC)))
end

@testset "laplacian_julia: analytic vs grid" begin
    # Moderate ζ (≲ Ecut/4) to keep the Gaussian resolved on the grid
    # (avoids the phase-A numerical-escape pathology).
    pairs = [
        (W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), 1.0),
         W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), 0.5)),
        (W2G.GaussianPolynomial([(1, 0, 1)], [1.0], zeros(3), 1.5),
         W2G.GaussianPolynomial([(0, 1, 1)], [1.0], [0.3, 0.0, 0.0], 0.8)),
    ]
    for (X1, X2) in pairs
        @test W2G.integral(X1, X2; type=:laplacian) ≈ _grid_laplacian(X1, X2) rtol = 1e-3
    end
end

@testset "laplacian_julia: symmetry" begin
    X1 = W2G.GaussianPolynomial([(1, 0, 1)], [1.0], zeros(3), 1.5)
    X2 = W2G.GaussianPolynomial([(0, 1, 1)], [1.0], [0.3, 0.0, 0.0], 0.8)
    @test W2G.integral(X1, X2; type=:laplacian) ≈ W2G.integral(X2, X1; type=:laplacian) rtol = 1e-12
end

@testset "Hˢ_overlap(Ms; s=1) is GaIn-free and SPD" begin
    Ms = [W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), ζ) for ζ in (0.5, 1.0, 2.0)]
    S = W2G.Hˢ_overlap(Ms; s=1)
    @test issymmetric(round.(S; digits=10))
    @test isposdef(Symmetric(S))
end
```

- [ ] **Step 2: Run to verify**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: all testsets PASS. (The fixture build takes ~tens of seconds.) If the `analytic vs grid` testset fails by a sign, the convention in `laplacian_julia` is inverted — revisit the trailing `-acc`.

- [ ] **Step 3: Commit**

```bash
git add test/test_integrals_laplacian.jl
git commit -m "TB integrals: grid-coherence, symmetry, H¹-overlap tests for laplacian_julia

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: `compare_onsite_kinetic` harness + integration test

**Files:**
- Modify: `src/tight_binding/hamiltonian_scalar_product.jl` (append)
- Modify: `test/test_integrals_laplacian.jl` (append integration testset)

- [ ] **Step 1: Write the failing integration test**

Single normalized s-Gaussian wrapped as a one-function `CompressedWannier`; its own grid transform is the "reference" Wannier. Reference grid kinetic, analytic Gaussian kinetic, and the closed form `3ζ/2` must all agree.

Append to `test/test_integrals_laplacian.jl`:

```julia
@testset "compare_onsite_kinetic: single-Gaussian triple agreement" begin
    ζ = 1.0
    X = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)   # normalized
    Φ = W2G.BasisFunction([1.0], [X])
    Wc = W2G.CompressedWannier(zeros(3), [Φ], [1.0],
                               _BASIS_SC, ComplexF64[], ComplexF64[], 0.0, 0.0)
    w_fourier = normalize(X(_BASIS_SC))   # reference == the same s-Gaussian
    out = W2G.compare_onsite_kinetic(Wc, w_fourier, _BASIS_SC)
    @test out.T_gto ≈ 3ζ / 2  rtol = 1e-10
    @test out.rel_err < 1e-3            # grid vs analytic, same function
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: FAIL — `compare_onsite_kinetic` is undefined.

- [ ] **Step 3: Implement the harness**

Append to `src/tight_binding/hamiltonian_scalar_product.jl`:

```julia
"""
On-site one-body kinetic comparison for the TB validation chain.

Contrasts the analytic kinetic energy of the Gaussian-compressed Wannier `Wc`
against the grid kinetic energy of the reference plane-wave Wannier `w_fourier`
(both on `basis_supercell`):

  - reference (grid):   ⟨w|−½∇²|w⟩ = ½ ⟨w, |G|² w⟩    (∇² ↦ −|G|² in Fourier)
  - Gaussian (analytic): integral(Wc, Wc; type=:kinetic)

Returns a NamedTuple with both values, the relative error, and the L² norms
(for context — the Gaussian Wannier need not be exactly normalized).
"""
function compare_onsite_kinetic(Wc::CompressedWannier, w_fourier::AbstractVector,
                                basis_supercell::PlaneWaveBasis)
    kpt = only(basis_supercell.kpoints)
    G2 = [sum(abs2, q) for q in G_vectors_cart(basis_supercell, kpt)]
    T_ref = 0.5 * real(dot(w_fourier, G2 .* w_fourier))
    T_gto = integral(Wc, Wc; type=:kinetic)
    rel_err = abs(T_gto - T_ref) / abs(T_ref)
    (; T_ref, T_gto, rel_err,
       norm_ref = norm(w_fourier),
       norm_gto = sqrt(real(integral(Wc, Wc; type=:overlap))))
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=. test/test_integrals_laplacian.jl`
Expected: all testsets PASS.

- [ ] **Step 5: Commit**

```bash
git add src/tight_binding/hamiltonian_scalar_product.jl test/test_integrals_laplacian.jl
git commit -m "TB validation: compare_onsite_kinetic harness (Wannier grid vs GTO analytic)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: `workflow/validate_tb_onsite.jl` (cluster run)

**Files:**
- Create: `workflow/validate_tb_onsite.jl`

No automated test — this is a data-loading script the user runs on the cluster (Ecut=50). It must load the **two** files (the compressed Wannier carries an empty `wannier` field after JSON round-trip; the reference comes from the `wannier_pz` file) on a **shared** supercell basis.

- [ ] **Step 1: Write the script**

Create `workflow/validate_tb_onsite.jl`:

```julia
# On-site one-body kinetic validation of the TB chain (Ecut=50, cluster run).
# Compares ⟨w|−½∇²|w⟩ for the true plane-wave Wannier vs the Gaussian-
# compressed Wannier. See notes/06-design-validation-TB-onsite.md.
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using LinearAlgebra

const ECUT = 50
dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]..., "workflow")
include(joinpath(dir, "monolayer_graphene.jl"))   # brings `Graphene` into scope

# Shared Γ-point supercell basis (must match the wannierization grid).
basis_sc = let
    basis = Graphene(; d=10u"Å", kgrid=[5, 5, 1], Ecut=ECUT).basis()
    DFTK.cell_to_supercell(basis)
end

# Reference Wannier (true plane-wave, Fourier) and Gaussian-compressed Wannier.
ref = W2G.read_wannier_function(joinpath(dir, "wannier_functions/wannier_pz_Ecut-$ECUT.json"))
w_fourier = normalize(ref.wannier)
Wc = W2G.CompressedWannier(joinpath(splitpath(dir)[1:end-1]..., "data",
                                     "compressed_wannier_55_H1_Ecut$ECUT.json"))

out = W2G.compare_onsite_kinetic(Wc, w_fourier, basis_sc)

println("On-site kinetic ⟨w|-½∇²|w⟩ @ Ecut=$ECUT")
println("  reference (grid)   T_ref   = ", out.T_ref)
println("  Gaussian (analytic) T_gto  = ", out.T_gto)
println("  relative error             = ", out.rel_err)
println("  ‖w‖ ref / gto              = ", out.norm_ref, " / ", out.norm_gto)
```

- [ ] **Step 2: Sanity check the script loads (local, low Ecut optional)**

The user runs this on the cluster at Ecut=50. Optionally, a local smoke pass at `ECUT = 15` (data files present) confirms paths and wiring before the cluster job. No assertion here — visual inspection of the printed relative error.

- [ ] **Step 3: Commit**

```bash
git add workflow/validate_tb_onsite.jl
git commit -m "TB validation: on-site kinetic workflow script (Ecut=50, cluster)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Out of scope (follow-up plans)

- **L-2b — V_KS term**: add the potential contribution to the on-site comparison (reference via real-space ∫w²V_KS; Gaussian via `potential_scalar_prod(Wc,Wc)`), resolving the cell↔supercell wiring of the KS potential. Needs `scfres` from `Graphene().scf()`.
- Hoppings t(R≠0), band reconstruction (T-band), bilayer/TBG, two-body terms.
