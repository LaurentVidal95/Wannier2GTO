# Phase B Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement joint optimization of SAGTO basis functions for monolayer graphene pz wannier compression, replacing the greedy procedure to break through the 8.6% H¹ ceiling identified in phase A.

**Architecture:** Variable projection (closed-form Tikhonov-regularized solve for linear coefficients $c_i$), L-BFGS with K independent restarts on shape parameters. Two function types: centered (D3-invariant) and π-bond (rotation-equivalent triplets). Box-via-sigmoid bounded log spreads.

**Tech Stack:** Julia 1.11, DFTK (planewave basis), Optim.jl (L-BFGS), Zygote (AD), ChainRulesCore. New module organization under `src/joint_optim/`.

**Reference spec:** `notes/03-design-phase-B.md`. All design decisions are locked there; this plan is the build sequence.

---

## File Structure

**To create:**
- `src/joint_optim/parametrization.jl` — `JointLayout` struct, flat ↔ structured params, sigmoid encoding, conversion to `BasisFunction`s
- `src/joint_optim/init.jl` — heuristic init for restarts (RNG-driven)
- `src/joint_optim/loss.jl` — variable projection inner solve + H¹ residual loss
- `src/joint_optim/runner.jl` — single L-BFGS run, K-restarts orchestration, JSON persistence
- `workflow/joint_compression.jl` — top-level driver script with smoke + run configs
- `test/test_parametrization.jl` — sigmoid round-trip, basis construction
- `test/test_loss.jl` — manual loss check on toy example, gradient finiteness
- `test/test_init.jl` — bounds & determinism

**To modify:**
- `src/Wannier2GTO.jl` — `include` new files, export top-level entry points
- `Project.toml` — add Zygote, ChainRulesCore, Random (via Pkg, not by hand)
- `.gitignore` — add `workflow/joint_outputs/`

---

## Task 1: Dependencies and module skeleton

**Files:**
- Modify: `Project.toml` (via Pkg)
- Modify: `src/Wannier2GTO.jl`
- Create: `src/joint_optim/parametrization.jl` (stub)
- Create: `src/joint_optim/init.jl` (stub)
- Create: `src/joint_optim/loss.jl` (stub)
- Create: `src/joint_optim/runner.jl` (stub)
- Modify: `.gitignore`

- [ ] **Step 1: Add Zygote and ChainRulesCore via Pkg**

```bash
julia --project=. -e 'using Pkg; Pkg.add(["Zygote", "ChainRulesCore"])'
```
Expected: both packages added to `[deps]` in `Project.toml`, no errors.

- [ ] **Step 2: Create stub files for the four new modules**

```bash
mkdir -p src/joint_optim
```

Create `src/joint_optim/parametrization.jl` with:
```julia
# JointLayout, flatten/unflatten, sigmoid encoding, params → BasisFunctions
# (filled in Task 2)
```

Same one-line stub for `init.jl`, `loss.jl`, `runner.jl`.

- [ ] **Step 3: Wire the new files into the main module**

Modify `src/Wannier2GTO.jl`. After `include("integrals/integrals.jl")` and before the tight-binding includes, add:
```julia
# Phase B: joint optimization of SAGTOs
export JointLayout, init_params, joint_loss, run_joint_optim
include("joint_optim/parametrization.jl")
include("joint_optim/init.jl")
include("joint_optim/loss.jl")
include("joint_optim/runner.jl")
```

Also add `using Random` and `using Zygote` near the top of the file (below `using ForwardDiff`).

- [ ] **Step 4: Add `joint_outputs` to `.gitignore`**

Append to `.gitignore`:
```
workflow/joint_outputs
```

- [ ] **Step 5: Verify the package precompiles**

```bash
julia --project=. -e 'using Wannier2GTO; println("OK")'
```
Expected: prints `OK` after precompilation, no errors.

- [ ] **Step 6: Commit**

```bash
git add Project.toml src/Wannier2GTO.jl src/joint_optim .gitignore
git commit -m "Phase B setup: deps and module skeleton"
```

---

## Task 2: `JointLayout` and flat ↔ structured parameter conversion

**Files:**
- Modify: `src/joint_optim/parametrization.jl`
- Create: `test/test_parametrization.jl`

**Math:** A `JointLayout` describes the structure (counts), and `JointParams` (or just a flat vector) holds values. For each *centered* function: latent $u_i$ (1) + polynomial coefs $\lambda$ (6) = 7 floats. For each *π-bond* function: $r_i$ (1) + $u_i$ (1) + polynomial coefs $\lambda$ (4) = 6 floats. The flat vector concatenates all centered then all π-bond entries.

- [ ] **Step 1: Write failing test for layout & flat dim**

Create `test/test_parametrization.jl`:
```julia
using Test
using Wannier2GTO

@testset "JointLayout" begin
    layout = JointLayout(N_centered=3, N_pibond=2,
                         poly_dim_centered=6, poly_dim_pibond=4)
    # 3*(1+6) + 2*(1+1+4) = 21 + 12 = 33
    @test flat_dim(layout) == 33
end
```

- [ ] **Step 2: Run the test, expect failure**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: errors with `JointLayout` undefined.

- [ ] **Step 3: Implement `JointLayout` and `flat_dim`**

Add to `src/joint_optim/parametrization.jl`:
```julia
"""
Describes the structure of a joint parameterization. Used to encode/decode a
flat parameter vector into structured per-function parameters.

`poly_dim_centered`/`poly_dim_pibond` are the number of polynomial groups
returned by `select_orders` for the D3-symmetric / non-D3-symmetric cases
(typically 6 and 4 with `max_xy_order=3, max_z_order=3`).
"""
Base.@kwdef struct JointLayout
    N_centered::Int
    N_pibond::Int
    poly_dim_centered::Int = 6
    poly_dim_pibond::Int = 4
end

"Number of floats in the flat parameter vector for this layout."
flat_dim(layout::JointLayout) =
    layout.N_centered * (1 + layout.poly_dim_centered) +
    layout.N_pibond  * (1 + 1 + layout.poly_dim_pibond)
```

- [ ] **Step 4: Run the test, expect pass**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: 1 test passed.

- [ ] **Step 5: Add structured access functions with tests**

Append to `test/test_parametrization.jl`:
```julia
@testset "flat ↔ structured" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = collect(1.0:flat_dim(layout))   # 2*(1+6) + 1*(1+1+4) = 20

    # Centered Φ_1: u = 1.0, λ = [2, 3, 4, 5, 6, 7]
    u_c, λ_c = get_centered(flat, layout, 1)
    @test u_c == 1.0
    @test λ_c == [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    # Centered Φ_2: u = 8.0, λ = [9, 10, 11, 12, 13, 14]
    u_c2, _ = get_centered(flat, layout, 2)
    @test u_c2 == 8.0

    # π-bond Φ_1: r = 15.0, u = 16.0, λ = [17, 18, 19, 20]
    r_p, u_p, λ_p = get_pibond(flat, layout, 1)
    @test r_p == 15.0 && u_p == 16.0
    @test λ_p == [17.0, 18.0, 19.0, 20.0]
end
```

Implement in `parametrization.jl`:
```julia
"Slice the flat vector to extract the i-th centered function (1-indexed)."
function get_centered(flat::AbstractVector, layout::JointLayout, i::Int)
    @assert 1 ≤ i ≤ layout.N_centered
    block = 1 + layout.poly_dim_centered
    offset = (i - 1) * block
    u = flat[offset + 1]
    λ = flat[offset + 2 : offset + block]
    (u, λ)
end

"Slice the flat vector to extract the j-th π-bond function (1-indexed)."
function get_pibond(flat::AbstractVector, layout::JointLayout, j::Int)
    @assert 1 ≤ j ≤ layout.N_pibond
    centered_total = layout.N_centered * (1 + layout.poly_dim_centered)
    block = 2 + layout.poly_dim_pibond
    offset = centered_total + (j - 1) * block
    r = flat[offset + 1]
    u = flat[offset + 2]
    λ = flat[offset + 3 : offset + block]
    (r, u, λ)
end
```

- [ ] **Step 6: Run the tests, expect pass**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: 2 testsets, all green.

- [ ] **Step 7: Commit**

```bash
git add src/joint_optim/parametrization.jl test/test_parametrization.jl
git commit -m "Phase B: JointLayout struct and flat ↔ structured conversion"
```

---

## Task 3: Box-via-sigmoid encoding for log ζ

**Files:**
- Modify: `src/joint_optim/parametrization.jl`
- Modify: `test/test_parametrization.jl`

**Math:** $\log\zeta = \log\zeta_\min + (\log\zeta_\max - \log\zeta_\min)\cdot\sigma(u)$ with $\sigma(u) = 1/(1+e^{-u})$. Inverse: $u = \mathrm{logit}((\log\zeta - \log\zeta_\min)/(\log\zeta_\max - \log\zeta_\min))$.

- [ ] **Step 1: Write failing test for sigmoid round-trip**

Append to `test/test_parametrization.jl`:
```julia
@testset "log ζ encoding round-trip" begin
    log_ζ_min = log(1e-2)
    log_ζ_max = log(4.0)
    for ζ in (0.05, 0.5, 1.0, 3.0, 3.99)
        log_ζ = log(ζ)
        u = log_ζ_to_u(log_ζ, log_ζ_min, log_ζ_max)
        log_ζ_back = log_ζ_from_u(u, log_ζ_min, log_ζ_max)
        @test log_ζ_back ≈ log_ζ rtol=1e-12
    end
    # Out-of-bounds u stays in the box
    log_ζ_huge = log_ζ_from_u(100.0, log_ζ_min, log_ζ_max)
    @test log_ζ_huge ≈ log_ζ_max rtol=1e-10
    log_ζ_neg = log_ζ_from_u(-100.0, log_ζ_min, log_ζ_max)
    @test log_ζ_neg ≈ log_ζ_min rtol=1e-10
end
```

- [ ] **Step 2: Run test, expect failure**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: `log_ζ_from_u` undefined.

- [ ] **Step 3: Implement encoding/decoding**

Add to `src/joint_optim/parametrization.jl`:
```julia
"Standard logistic sigmoid."
@inline _sigmoid(u) = inv(one(u) + exp(-u))

"Logit (inverse sigmoid), valid for x ∈ (0, 1)."
@inline _logit(x) = log(x / (one(x) - x))

"""
Decode a latent `u ∈ ℝ` into `log ζ ∈ [log_ζ_min, log_ζ_max]` via
`log ζ = log_ζ_min + (log_ζ_max - log_ζ_min) · σ(u)`.
"""
@inline function log_ζ_from_u(u, log_ζ_min, log_ζ_max)
    log_ζ_min + (log_ζ_max - log_ζ_min) * _sigmoid(u)
end

"""
Encode a target `log ζ ∈ (log_ζ_min, log_ζ_max)` to its latent `u`.
"""
@inline function log_ζ_to_u(log_ζ, log_ζ_min, log_ζ_max)
    x = (log_ζ - log_ζ_min) / (log_ζ_max - log_ζ_min)
    _logit(x)
end
```

- [ ] **Step 4: Run test, expect pass**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: 3 testsets, all green.

- [ ] **Step 5: Commit**

```bash
git add src/joint_optim/parametrization.jl test/test_parametrization.jl
git commit -m "Phase B: box-via-sigmoid encoding for log ζ"
```

---

## Task 4: Flat params → Vector{BasisFunction} conversion

**Files:**
- Modify: `src/joint_optim/parametrization.jl`
- Modify: `test/test_parametrization.jl`

**Math:** Given `JointLayout`, the flat vector, the π-bond unit direction $\hat{u}_\pi$, and the bounds, build the `Vector{BasisFunction}` of size $N = N_\text{centered} + N_\text{pibond}$. Each centered function: one polynomial group per allowed order with coef $\lambda$, all sharing center 0 and spread $\zeta$. Each π-bond function: same with center $r\,\hat{u}_\pi$, plus D3 enforcement is applied automatically in the conversion (so that downstream `BasisFunction`s already have all 3 rotational copies).

- [ ] **Step 1: Write failing test for shape and counts**

Append to `test/test_parametrization.jl`:
```julia
@testset "params_to_basis_functions: shape" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = zeros(flat_dim(layout))    # all zeros
    log_ζ_min, log_ζ_max = log(1e-2), log(4.0)
    π_bond_dir = [cos(-π/6), sin(-π/6), 0.0]

    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_dir,
                                    max_xy_order=3, max_z_order=3)
    @test length(Φs) == 3
    # Centered Φ_1, Φ_2: 6 SAGTOs each
    @test length(Φs[1].SAGTOs) == 6
    @test length(Φs[2].SAGTOs) == 6
    # π-bond Φ_3: 4 SAGTOs × 3 D3 rotations = 12 SAGTOs
    @test length(Φs[3].SAGTOs) == 12
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: `params_to_basis_functions` undefined.

- [ ] **Step 3: Implement the conversion**

Add to `src/joint_optim/parametrization.jl`:
```julia
"""
Convert a flat parameter vector to a `Vector{BasisFunction}` of length
`N = N_centered + N_pibond`. The first `N_centered` are D3-invariant centered
functions; the remaining `N_pibond` are π-bond functions whose D3 symmetry
has already been enforced (each contains 3 rotational copies of its SAGTOs).

Parameters
- `flat`: parameter vector of length `flat_dim(layout)`
- `layout`: structure descriptor
- `log_ζ_min`, `log_ζ_max`: spread box bounds (in log-scale)
- `π_bond_unit`: unit vector along the π-bond direction in 3D (z-component zero)
- `max_xy_order`, `max_z_order`: passed to `select_orders` for SAGTO basis
"""
function params_to_basis_functions(flat::AbstractVector, layout::JointLayout;
                                   log_ζ_min, log_ζ_max,
                                   π_bond_unit::AbstractVector,
                                   max_xy_order=3, max_z_order=3)
    Φs = BasisFunction[]

    # Centered functions
    xy_c, z_c = select_orders(max_xy_order, max_z_order, true)
    for i in 1:layout.N_centered
        u, λ = get_centered(flat, layout, i)
        ζ = exp(log_ζ_from_u(u, log_ζ_min, log_ζ_max))
        SAGTOs = SAGTO_basis(zeros(eltype(flat), 3), ζ, xy_c, z_c)
        # Centered functions are inherently D3-invariant — no enforcement
        push!(Φs, BasisFunction(collect(λ), SAGTOs))
    end

    # π-bond functions, D3 enforced by construction (3 rotational copies)
    xy_p, z_p = select_orders(max_xy_order, max_z_order, false)
    for j in 1:layout.N_pibond
        r, u, λ = get_pibond(flat, layout, j)
        ζ = exp(log_ζ_from_u(u, log_ζ_min, log_ζ_max))
        center = r .* π_bond_unit
        SAGTOs = SAGTO_basis(center, ζ, xy_p, z_p)
        Φ = BasisFunction(collect(λ), SAGTOs)
        push!(Φs, enforce_D3_symmetry(Φ))
    end

    Φs
end
```

- [ ] **Step 4: Run test, expect pass**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: 4 testsets, all green.

- [ ] **Step 5: Add a regression test for centers and spreads**

Append to `test/test_parametrization.jl`:
```julia
@testset "params_to_basis_functions: centers and spreads" begin
    layout = JointLayout(N_centered=1, N_pibond=1)
    flat = zeros(flat_dim(layout))
    flat[1] = 0.0                  # centered u = 0 -> log ζ at midpoint
    flat[8] = 2.0                  # π-bond r = 2.0
    flat[9] = 1.0                  # π-bond u = 1.0
    log_ζ_min, log_ζ_max = log(1e-2), log(4.0)
    π_bond_dir = [cos(-π/6), sin(-π/6), 0.0]

    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_dir,
                                    max_xy_order=3, max_z_order=3)
    # Φ_1 centered at origin
    @test Φs[1].SAGTOs[1].center ≈ zeros(3)
    # Φ_2 first SAGTO at r * π_bond_dir, the others at the rotated positions
    @test Φs[2].SAGTOs[1].center ≈ 2.0 .* π_bond_dir rtol=1e-12
    # Spread reflects sigmoid(u=1.0) ≈ 0.731
    expected_log_ζ = log(1e-2) + (log(4.0) - log(1e-2)) * (1 / (1 + exp(-1.0)))
    @test log(Φs[2].SAGTOs[1].spread) ≈ expected_log_ζ rtol=1e-10
end
```

- [ ] **Step 6: Run, expect pass**

```bash
julia --project=. test/test_parametrization.jl
```
Expected: 5 testsets, all green.

- [ ] **Step 7: Commit**

```bash
git add src/joint_optim/parametrization.jl test/test_parametrization.jl
git commit -m "Phase B: params → Vector{BasisFunction} conversion with D3 enforcement"
```

---

## Task 5: Init heuristic per restart

**Files:**
- Modify: `src/joint_optim/init.jl`
- Create: `test/test_init.jl`

**Math:** Sample $u \sim \mathcal{N}(0, 1.5)$ for each spread, $r \sim \mathcal{U}([0.5, 5.0])$ for each π-bond, $\lambda$ initialized at table-5.2 values + $\mathcal{N}(0, 0.05)$.

- [ ] **Step 1: Write failing test for determinism and bounds**

Create `test/test_init.jl`:
```julia
using Test
using Random
using Wannier2GTO

@testset "init_params: determinism" begin
    layout = JointLayout(N_centered=3, N_pibond=2)
    flat1 = init_params(MersenneTwister(42), layout)
    flat2 = init_params(MersenneTwister(42), layout)
    @test flat1 == flat2
    flat3 = init_params(MersenneTwister(43), layout)
    @test flat1 != flat3
end

@testset "init_params: dimension and finiteness" begin
    layout = JointLayout(N_centered=4, N_pibond=3)
    flat = init_params(MersenneTwister(1), layout)
    @test length(flat) == flat_dim(layout)
    @test all(isfinite, flat)
end

@testset "init_params: r bounds" begin
    layout = JointLayout(N_centered=2, N_pibond=4)
    flat = init_params(MersenneTwister(7), layout; r_min=0.5, r_max=5.0)
    for j in 1:layout.N_pibond
        r, _, _ = get_pibond(flat, layout, j)
        @test 0.5 ≤ r ≤ 5.0
    end
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. test/test_init.jl
```
Expected: `init_params` undefined.

- [ ] **Step 3: Implement init_params**

Add to `src/joint_optim/init.jl`:
```julia
using Random

"""
Default polynomial coefficients (matching the symmetry-adapted table 5.2 of
the manuscript) for centered (D3) and non-centered functions, in the order
returned by `select_orders(max_xy_order=3, max_z_order=3, ·)`.

These provide a physically-motivated starting point that is then perturbed
per restart. For now we use a simple all-ones init; if a more carefully
chosen baseline proves useful, edit here.
"""
const _DEFAULT_λ_CENTERED = ones(6)
const _DEFAULT_λ_PIBOND   = ones(4)

"""
Sample a flat parameter vector for a single L-BFGS restart. Uses `rng` for
all random draws — pass a fresh `MersenneTwister(seed)` per restart for
diversity.

Defaults match `notes/03-design-phase-B.md` §5.
"""
function init_params(rng::AbstractRNG, layout::JointLayout;
                     u_std::Real = 1.5,
                     r_min::Real = 0.5,
                     r_max::Real = 5.0,
                     λ_perturb_std::Real = 0.05)
    flat = zeros(flat_dim(layout))

    block_c = 1 + layout.poly_dim_centered
    for i in 1:layout.N_centered
        offset = (i - 1) * block_c
        flat[offset + 1] = u_std * randn(rng)
        flat[offset + 2 : offset + block_c] .=
            _DEFAULT_λ_CENTERED .+ λ_perturb_std .* randn(rng, layout.poly_dim_centered)
    end

    centered_total = layout.N_centered * block_c
    block_p = 2 + layout.poly_dim_pibond
    for j in 1:layout.N_pibond
        offset = centered_total + (j - 1) * block_p
        flat[offset + 1] = r_min + (r_max - r_min) * rand(rng)
        flat[offset + 2] = u_std * randn(rng)
        flat[offset + 3 : offset + block_p] .=
            _DEFAULT_λ_PIBOND .+ λ_perturb_std .* randn(rng, layout.poly_dim_pibond)
    end

    flat
end
```

- [ ] **Step 4: Run tests, expect pass**

```bash
julia --project=. test/test_init.jl
```
Expected: 3 testsets, all green.

- [ ] **Step 5: Commit**

```bash
git add src/joint_optim/init.jl test/test_init.jl
git commit -m "Phase B: deterministic per-restart init heuristic"
```

---

## Task 6: Variable projection inner solve

**Files:**
- Modify: `src/joint_optim/loss.jl`

**Math:** Given Φs, $w_z$ Fourier, basis_supercell, and Tikhonov $\varepsilon$: compute $S = \langle\Phi_i, \Phi_j\rangle_{H^1}$, $\Gamma_i = \langle w_z, \Phi_i\rangle_{H^1}$, then $\mathbf{c}^\star = (S + \varepsilon I)^{-1}\Gamma$.

We can reuse `_make_inner_solver(:tikhonov, ε)` from BasisFunctions.jl which is already AD-friendly.

- [ ] **Step 1: Implement `joint_inner_solve`**

Add to `src/joint_optim/loss.jl`:
```julia
"""
Compute the optimal linear coefficients `c_i` via Tikhonov-regularized
projection onto the basis `Φs` in the H^s-inner-product space.

`s = 1` corresponds to H¹ (the regime used in compression).

This wraps `_make_inner_solver(:tikhonov, ε)` from BasisFunctions.jl. The
matrix `S` and the right-hand side `Γ` are computed in plane-wave Fourier.
Returned `c` has the same length as `Φs`.
"""
function joint_inner_solve(Φs::Vector{BasisFunction},
                           w_z_fourier::AbstractVector,
                           basis_supercell::PlaneWaveBasis;
                           s::Int = 1,
                           ε::Real = 1e-8)
    Φs_Four = [Φ(basis_supercell) for Φ in Φs]
    Γ = [Hˢ_dot(basis_supercell, w_z_fourier, Φ; s=s) for Φ in Φs_Four]
    S = Hˢ_overlap(basis_supercell, Φs_Four; s=s)
    solver = _make_inner_solver(:tikhonov, ε)
    c = solver(S, Γ)
    real(c)   # eigenvalues of S are real, c should be real up to noise
end
```

- [ ] **Step 2: Manual sanity check (no automated test yet — defer to Task 7 integration)**

```bash
julia --project=. -e '
using Wannier2GTO, DFTK, DFTK.Unitful
import Wannier2GTO as W2G
include("workflow/monolayer_graphene.jl")

basis = Graphene(; d=10u"Å", kgrid=[5,5,1], Ecut=15).basis()
sc = DFTK.cell_to_supercell(basis)
data = W2G.read_wannier_function("workflow/wannier_functions/wannier_pz_Ecut-15.json")

layout = JointLayout(N_centered=2, N_pibond=1)
flat = init_params(Random.MersenneTwister(1), layout)
log_ζ_min, log_ζ_max = log(1e-2), log(4.0)
π_dir = [cos(-1.31), sin(-1.31), 0.0]; π_dir ./= sqrt(sum(abs2, π_dir))

Φs = params_to_basis_functions(flat, layout;
        log_ζ_min, log_ζ_max, π_bond_unit=π_dir)
c = W2G.joint_inner_solve(Φs, normalize(data.wannier), sc)
println("c = ", c)
println("|c| = ", abs.(c))
'
```
Expected: prints a 3-element vector of floats, all finite, magnitude < 100.

- [ ] **Step 3: Commit**

```bash
git add src/joint_optim/loss.jl
git commit -m "Phase B: Tikhonov-regularized variable projection solve"
```

---

## Task 7: H¹ residual loss (joint_loss)

**Files:**
- Modify: `src/joint_optim/loss.jl`
- Create: `test/test_loss.jl`

**Math:** $\mathcal{L}(\theta) = \|w_z - \sum_i c_i^\star(\theta)\,\Phi_i(\theta)\|^2_{H^1}$. Computed in plane-wave space via Plancherel. The full pipeline takes a flat $\theta$, calls `params_to_basis_functions`, runs `joint_inner_solve`, then computes the residual norm.

- [ ] **Step 1: Write a baseline correctness test (loss > 0, finite, depends on params)**

Create `test/test_loss.jl`:
```julia
using Test
using Random
using DFTK
using DFTK.Unitful
using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra

# Heavy fixture: bring `Graphene` into scope and build the supercell basis
# once for all testsets in this file.
include(joinpath(@__DIR__, "..", "workflow", "monolayer_graphene.jl"))

const _BASIS_SC = let
    basis = Graphene(; d=10u"Å", kgrid=[5,5,1], Ecut=15).basis()
    DFTK.cell_to_supercell(basis)
end

const _DATA = W2G.read_wannier_function(joinpath(@__DIR__, "..",
                            "workflow/wannier_functions/wannier_pz_Ecut-15.json"))
const _W_FOURIER = normalize(_DATA.wannier)
const _π_BOND_DIR = let
    θ = _DATA.π_bond.θ
    [cos(θ), sin(θ), 0.0]
end

@testset "joint_loss: positive and finite" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = init_params(MersenneTwister(0), layout)
    L = joint_loss(flat, layout, _W_FOURIER, _BASIS_SC;
                   log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                   π_bond_unit=_π_BOND_DIR)
    @test isfinite(L)
    @test L > 0
end

@testset "joint_loss: different params → different losses" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat1 = init_params(MersenneTwister(0), layout)
    flat2 = init_params(MersenneTwister(1), layout)
    L1 = joint_loss(flat1, layout, _W_FOURIER, _BASIS_SC;
                    log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                    π_bond_unit=_π_BOND_DIR)
    L2 = joint_loss(flat2, layout, _W_FOURIER, _BASIS_SC;
                    log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                    π_bond_unit=_π_BOND_DIR)
    @test L1 != L2
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. test/test_loss.jl
```
Expected: `joint_loss` undefined.

- [ ] **Step 3: Implement joint_loss**

Add to `src/joint_optim/loss.jl`:
```julia
"""
Joint-fit loss: H¹ squared residual after Tikhonov-regularized
projection of `w_z` onto the basis described by `flat`.

Parameters
- `flat`: free parameters (length `flat_dim(layout)`)
- `layout`: parameter layout descriptor
- `w_z_fourier`: target wannier in plane-wave Fourier (supercell Γ point)
- `basis_supercell`: discretization basis
- `log_ζ_min`, `log_ζ_max`: spread box bounds (passed to encoding)
- `π_bond_unit`: unit vector along π-bond axis
- `s`: Sobolev exponent (default 1 for H¹)
- `ε`: Tikhonov regularization
- `max_xy_order`, `max_z_order`: SAGTO orders

Returns a positive real scalar.
"""
function joint_loss(flat::AbstractVector, layout::JointLayout,
                    w_z_fourier::AbstractVector,
                    basis_supercell::PlaneWaveBasis;
                    log_ζ_min, log_ζ_max,
                    π_bond_unit::AbstractVector,
                    s::Int = 1,
                    ε::Real = 1e-8,
                    max_xy_order::Int = 3,
                    max_z_order::Int = 3)
    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_unit,
                                    max_xy_order, max_z_order)
    Φs_Four = [Φ(basis_supercell) for Φ in Φs]
    Γ = [Hˢ_dot(basis_supercell, w_z_fourier, Φ; s=s) for Φ in Φs_Four]
    S = Hˢ_overlap(basis_supercell, Φs_Four; s=s)
    solver = _make_inner_solver(:tikhonov, ε)
    c = real(solver(S, Γ))
    residual = w_z_fourier - sum(c[i] .* Φs_Four[i] for i in eachindex(c))
    Hˢ_norm(basis_supercell, residual; s=s)^2
end
```

- [ ] **Step 4: Run tests, expect pass**

```bash
julia --project=. test/test_loss.jl
```
Expected: 2 testsets pass.

- [ ] **Step 5: Sanity-check the loss scale on a small case**

```bash
julia --project=. -e '
using Random; using Wannier2GTO; import Wannier2GTO as W2G
include("test/test_loss.jl")
layout = JointLayout(N_centered=8, N_pibond=7)
flat = init_params(MersenneTwister(0), layout)
L = joint_loss(flat, layout, _W_FOURIER, _BASIS_SC;
               log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
               π_bond_unit=_π_BOND_DIR)
println("Loss at random init (N=15): ", L)
println("  sqrt(L) = ", sqrt(L), "  (should be O(1) since w_z is L²-normalized)")
'
```
Expected: Loss in range [0.05, 1.0] at random init (much higher than greedy's 0.07² ≈ 5e-3, but finite and informative).

- [ ] **Step 6: Commit**

```bash
git add src/joint_optim/loss.jl test/test_loss.jl
git commit -m "Phase B: joint_loss with variable projection"
```

---

## Task 8: Zygote gradient verification

**Files:**
- Modify: `test/test_loss.jl`

**Goal:** Confirm Zygote can differentiate through the entire pipeline (params → SAGTOs → S\Γ → residual → norm). If this fails, we'll need to switch to Enzyme.

- [ ] **Step 1: Add a gradient finiteness test**

Append to `test/test_loss.jl`:
```julia
using Zygote

@testset "joint_loss: Zygote gradient is finite" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = init_params(MersenneTwister(0), layout)
    f = θ -> joint_loss(θ, layout, _W_FOURIER, _BASIS_SC;
                        log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                        π_bond_unit=_π_BOND_DIR)
    g = first(Zygote.gradient(f, flat))
    @test length(g) == flat_dim(layout)
    @test all(isfinite, g)
    @test norm(g) > 0
end
```

- [ ] **Step 2: Run the test**

```bash
julia --project=. test/test_loss.jl
```
Expected: gradient test passes. If it fails with a Zygote-related error (typically "Mutating arrays is not supported" or similar), follow Step 3 troubleshooting; otherwise skip to Step 4.

- [ ] **Step 3: Troubleshooting (only if Step 2 fails)**

If Zygote fails on a mutation, the most likely culprits in our codebase are:
- `Hˢ_overlap` with `for μ … for ν …` loop building a mutable matrix
- `enforce_D3_symmetry` with `vcat` and explicit array building

Workarounds (in priority):
1. Wrap problematic operations with `Zygote.Buffer` or convert loops to comprehensions
2. Replace `Hˢ_overlap` numerical loop with the comprehension form already at `Hs_scalar_prods.jl:36-40`
3. As a last resort, switch backend to Enzyme: replace `Zygote.gradient` with `Enzyme.gradient` and add Enzyme to deps

The key file likely to need adaptation is `src/common/Hs_scalar_prods.jl` line 24-34. Make `Hˢ_overlap` AD-friendly by replacing the mutating pattern with a comprehension:
```julia
function Hˢ_overlap(basis_supercell::PlaneWaveBasis, Xs_fourier; s=0)
    n = length(Xs_fourier)
    [Hˢ_dot(basis_supercell, Xs_fourier[i], Xs_fourier[j]; s) for i in 1:n, j in 1:n]
end
```
(Drops the `Symmetric` wrapper, but the matrix is symmetric by construction.)

Re-run the test, and if it still fails, document the actual Zygote error in `notes/04-plan-phase-B.md` and decide whether to switch to Enzyme before continuing.

- [ ] **Step 4: Commit (success path) or escalate (failure path)**

If gradient test passes:
```bash
git add test/test_loss.jl src/common/Hs_scalar_prods.jl  # second only if modified
git commit -m "Phase B: Zygote gradient verified end-to-end"
```

If gradient test still fails after Step 3 troubleshooting, **stop and report the issue** — do not proceed to Task 9 with broken AD.

---

## Task 9: Single L-BFGS run

**Files:**
- Modify: `src/joint_optim/runner.jl`

**Goal:** Wire up Optim.jl's L-BFGS with Zygote gradients on `joint_loss`. Returns the final flat params, final loss, iteration trace.

- [ ] **Step 1: Implement `run_lbfgs_once`**

Add to `src/joint_optim/runner.jl`:
```julia
using Optim
using LineSearches

"""
Result of a single L-BFGS run.

- `flat_final`: parameter vector at convergence
- `loss_final`: loss value at convergence
- `loss_history`: vector of loss values per iteration
- `iterations`: number of iterations actually run
- `converged`: did Optim report convergence?
- `wallclock_seconds`: wall-clock time
"""
struct LBFGSResult
    flat_final::Vector{Float64}
    loss_final::Float64
    loss_history::Vector{Float64}
    iterations::Int
    converged::Bool
    wallclock_seconds::Float64
end

"""
Run L-BFGS on `joint_loss` from `flat_init`, with Zygote-provided gradients.

Hyperparameters follow `notes/03-design-phase-B.md` §3 defaults.
"""
function run_lbfgs_once(flat_init::AbstractVector, layout::JointLayout,
                        w_z_fourier::AbstractVector,
                        basis_supercell::PlaneWaveBasis;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit::AbstractVector,
                        max_iter::Int = 200,
                        g_abstol::Real = 1e-5,
                        f_reltol::Real = 1e-8,
                        s::Int = 1,
                        ε::Real = 1e-8)
    f = θ -> joint_loss(θ, layout, w_z_fourier, basis_supercell;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit=π_bond_unit, s=s, ε=ε)
    g! = (G, θ) -> (G .= first(Zygote.gradient(f, θ)); nothing)

    history = Float64[]
    cb = function (state)
        push!(history, state.value)
        false   # keep going
    end

    options = Optim.Options(g_abstol=g_abstol, f_reltol=f_reltol,
                            iterations=max_iter, callback=cb,
                            show_trace=false)

    t0 = time()
    res = optimize(f, g!, collect(flat_init),
                   LBFGS(; linesearch=LineSearches.HagerZhang()),
                   options)
    elapsed = time() - t0

    LBFGSResult(Optim.minimizer(res), Optim.minimum(res),
                history, Optim.iterations(res),
                Optim.converged(res), elapsed)
end
```

- [ ] **Step 2: Smoke test the runner manually**

```bash
julia --project=. -e '
using Random, Wannier2GTO; import Wannier2GTO as W2G
include("test/test_loss.jl")
layout = JointLayout(N_centered=2, N_pibond=1)
flat0 = init_params(MersenneTwister(0), layout)
res = W2G.run_lbfgs_once(flat0, layout, _W_FOURIER, _BASIS_SC;
                         log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                         π_bond_unit=_π_BOND_DIR, max_iter=20)
println("Initial loss: ", first(res.loss_history))
println("Final loss:   ", res.loss_final)
println("Decreased?    ", res.loss_final < first(res.loss_history))
println("Iters:        ", res.iterations)
println("Wallclock:    ", round(res.wallclock_seconds; digits=2), " s")
'
```
Expected: loss decreases, iterations between 5 and 20, wallclock ~10–60 s.

- [ ] **Step 3: Commit**

```bash
git add src/joint_optim/runner.jl
git commit -m "Phase B: single-run L-BFGS with Zygote gradients"
```

---

## Task 10: K-restarts orchestration and persistence

**Files:**
- Modify: `src/joint_optim/runner.jl`

**Goal:** Loop K independent restarts (different RNG seeds), keep best, save per-restart and aggregate results to JSON in `workflow/joint_outputs/run_<timestamp>/`.

- [ ] **Step 1: Implement `run_joint_optim`**

Add to `src/joint_optim/runner.jl`:
```julia
using JSON3
using Dates

"""
Aggregate result of a K-restart joint optimization run.
"""
struct JointOptimResult
    layout::JointLayout
    K::Int
    seeds::Vector{Int}
    restart_results::Vector{LBFGSResult}
    best_idx::Int
    output_dir::String
end

"""
Run K independent L-BFGS restarts, save artifacts, return the best result.

Output layout under `output_root/run_<timestamp>/`:
  - config.json     : hyperparameters and seeds
  - restart_<k>.json: per-restart trace (loss curve + final flat)
  - best.json       : best restart payload
  - summary.txt     : human-readable comparison
"""
function run_joint_optim(layout::JointLayout,
                         w_z_fourier::AbstractVector,
                         basis_supercell::PlaneWaveBasis;
                         K::Int = 5,
                         master_seed::Int = 42,
                         log_ζ_min, log_ζ_max,
                         π_bond_unit::AbstractVector,
                         output_root::String = joinpath(@__DIR__,
                                                        "..", "..",
                                                        "workflow", "joint_outputs"),
                         max_iter::Int = 200,
                         g_abstol::Real = 1e-5,
                         f_reltol::Real = 1e-8,
                         s::Int = 1,
                         ε::Real = 1e-8)
    timestamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    output_dir = joinpath(output_root, "run_" * timestamp)
    mkpath(output_dir)

    seeds = [master_seed + k for k in 0:(K - 1)]
    results = LBFGSResult[]

    for (k, seed) in enumerate(seeds)
        println("--- Restart $k/$K (seed=$seed) ---"); flush(stdout)
        flat0 = init_params(MersenneTwister(seed), layout)
        res = run_lbfgs_once(flat0, layout, w_z_fourier, basis_supercell;
                             log_ζ_min, log_ζ_max,
                             π_bond_unit=π_bond_unit,
                             max_iter, g_abstol, f_reltol, s, ε)
        push!(results, res)
        @info "  loss_init=$(first(res.loss_history))  loss_final=$(res.loss_final)" *
              "  iters=$(res.iterations)  wallclock=$(round(res.wallclock_seconds; digits=1))s"
        # Per-restart persistence
        open(joinpath(output_dir, "restart_$(lpad(k, 2, '0')).json"), "w") do io
            JSON3.write(io,
                Dict("seed" => seed,
                     "loss_history" => res.loss_history,
                     "loss_final" => res.loss_final,
                     "iterations" => res.iterations,
                     "converged" => res.converged,
                     "wallclock_seconds" => res.wallclock_seconds,
                     "flat_final" => res.flat_final))
        end
    end

    best_idx = argmin([r.loss_final for r in results])

    # Aggregate config and best
    open(joinpath(output_dir, "config.json"), "w") do io
        JSON3.write(io,
            Dict("K" => K, "master_seed" => master_seed,
                 "seeds" => seeds, "max_iter" => max_iter,
                 "log_ζ_min" => log_ζ_min, "log_ζ_max" => log_ζ_max,
                 "ε" => ε, "s" => s,
                 "N_centered" => layout.N_centered,
                 "N_pibond" => layout.N_pibond,
                 "poly_dim_centered" => layout.poly_dim_centered,
                 "poly_dim_pibond" => layout.poly_dim_pibond))
    end

    open(joinpath(output_dir, "best.json"), "w") do io
        JSON3.write(io,
            Dict("best_restart" => best_idx,
                 "best_seed" => seeds[best_idx],
                 "best_loss" => results[best_idx].loss_final,
                 "best_flat" => results[best_idx].flat_final))
    end

    open(joinpath(output_dir, "summary.txt"), "w") do io
        println(io, "Joint optim summary — $(timestamp)")
        println(io, "K = $K, layout = (N_centered=$(layout.N_centered), N_pibond=$(layout.N_pibond))")
        println(io, "")
        println(io, "Restart  Seed     Loss_init      Loss_final     Iters  Wall(s)")
        for (k, r) in enumerate(results)
            println(io, lpad(k, 6), "  ",
                    lpad(seeds[k], 5), "  ",
                    lpad(round(first(r.loss_history); digits=4), 12), "  ",
                    lpad(round(r.loss_final; digits=6), 14), "  ",
                    lpad(r.iterations, 5), "  ",
                    lpad(round(r.wallclock_seconds; digits=1), 7))
        end
        println(io, "")
        println(io, "Best: restart $best_idx with loss = $(results[best_idx].loss_final)")
    end

    JointOptimResult(layout, K, seeds, results, best_idx, output_dir)
end
```

- [ ] **Step 2: Smoke test K=2 on a tiny configuration**

```bash
julia --project=. -e '
using Random, Wannier2GTO; import Wannier2GTO as W2G
include("test/test_loss.jl")
layout = JointLayout(N_centered=2, N_pibond=1)
result = W2G.run_joint_optim(layout, _W_FOURIER, _BASIS_SC;
                             K=2, master_seed=100,
                             log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                             π_bond_unit=_π_BOND_DIR, max_iter=10)
println("Output dir: ", result.output_dir)
println("Best restart: ", result.best_idx)
println("Best loss: ", result.restart_results[result.best_idx].loss_final)
run(`ls $(result.output_dir)`)
'
```
Expected: output dir contains `config.json`, `restart_01.json`, `restart_02.json`, `best.json`, `summary.txt`.

- [ ] **Step 3: Commit**

```bash
git add src/joint_optim/runner.jl
git commit -m "Phase B: K-restarts loop with JSON persistence"
```

---

## Task 11: End-to-end smoke test workflow

**Files:**
- Create: `workflow/joint_compression.jl`

**Goal:** A single CLI script that takes Ecut and N as inputs, runs the joint optim, prints summary. Smoke test with Ecut=15, N=10, K=3.

- [ ] **Step 1: Write the script**

Create `workflow/joint_compression.jl`:
```julia
#
# Phase B end-to-end joint optimization driver.
#
# Reads a precomputed wannier from JSON, reconstructs the matching plane-wave
# basis, runs K independent L-BFGS restarts on the joint loss, saves artifacts.
#
# Usage:
#   julia --project=. -t auto workflow/joint_compression.jl
#
# Tweak the configuration block below.
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
const ECUT          = 15
const KGRID         = [5, 5, 1]
const D             = 10.0u"Å"
const WANNIER_JSON  = joinpath(@__DIR__,
                                "wannier_functions/wannier_pz_Ecut-15.json")
const N_CENTERED    = 5    # smoke: small total
const N_PIBOND      = 5
const K_RESTARTS    = 3
const MASTER_SEED   = 42
const MAX_ITER      = 100
const ζ_MIN         = 1e-2
const ζ_MAX         = ECUT / 4.0
const ε_TIKHONOV    = 1e-8

# --- Setup ----------------------------------------------------------------
println("=" ^ 70)
println("Phase B joint compression — smoke run")
println("=" ^ 70)
println("Ecut=$(ECUT)  kgrid=$(KGRID)  d=$(D)")
println("N_centered=$(N_CENTERED) N_pibond=$(N_PIBOND) K=$(K_RESTARTS)")
println("ζ box: [$(ζ_MIN), $(ζ_MAX)]   ε=$(ε_TIKHONOV)")
println()

println("Building plane-wave basis...")
basis = Graphene(; d=D, kgrid=KGRID, Ecut=ECUT).basis()
basis_supercell = DFTK.cell_to_supercell(basis)
println("  Γ G-vectors: ",
        length(DFTK.G_vectors(basis_supercell, basis_supercell.kpoints[1])))

println("Loading wannier from $(WANNIER_JSON)")
data = W2G.read_wannier_function(WANNIER_JSON)
w_z_fourier = normalize(data.wannier)
π_bond_unit = let θ = data.π_bond.θ
    [cos(θ), sin(θ), 0.0]
end

# --- Run ------------------------------------------------------------------
layout = JointLayout(N_centered=N_CENTERED, N_pibond=N_PIBOND)
println("Layout: ", layout)
println("flat_dim = ", W2G.flat_dim(layout))
println()

result = run_joint_optim(layout, w_z_fourier, basis_supercell;
                         K=K_RESTARTS, master_seed=MASTER_SEED,
                         log_ζ_min=log(ζ_MIN), log_ζ_max=log(ζ_MAX),
                         π_bond_unit=π_bond_unit,
                         max_iter=MAX_ITER, ε=ε_TIKHONOV)

# --- Report ---------------------------------------------------------------
println()
println("=" ^ 70)
println("RESULT")
println("=" ^ 70)
best = result.restart_results[result.best_idx]
println("Best restart index: ", result.best_idx, " (seed=$(result.seeds[result.best_idx]))")
@printf("Best loss (squared H¹): %.6e\n", best.loss_final)
@printf("Best H¹ relative error: %.4f%%\n", 100 * sqrt(best.loss_final))
println("All artifacts: ", result.output_dir)
```

- [ ] **Step 2: Run the smoke test**

```bash
julia --project=. -t auto workflow/joint_compression.jl 2>&1 | tee /tmp/joint_smoke.log
```
Expected: completes within ~10 min on laptop, prints final relative error. Target for smoke: H¹ relative error < 30%, just to validate end-to-end. If error is < 10%, that's excellent for a smoke run with N=10 / K=3.

- [ ] **Step 3: Inspect artifacts**

```bash
ls workflow/joint_outputs/run_*/
cat workflow/joint_outputs/run_*/summary.txt
```
Expected: directory contains 4 files (`config.json`, `best.json`, `summary.txt`, plus K=3 restart_*.json), summary shows monotone loss decrease across all restarts.

- [ ] **Step 4: Commit**

```bash
git add workflow/joint_compression.jl
git commit -m "Phase B: joint compression driver script with smoke run"
```

---

## Task 12: Compare against phase A baseline

**Files:**
- Modify: `workflow/joint_compression.jl` (configuration only)
- Create: `notes/05-bilan-phase-B-smoke.md`

**Goal:** Run the script with N=15 (8 centered, 7 π-bond) and K=5, then write up the comparison with the phase A greedy baseline (8.6% H¹ at 11 functions).

- [ ] **Step 1: Bump configuration to baseline**

Edit `workflow/joint_compression.jl`:
```julia
const N_CENTERED    = 8
const N_PIBOND      = 7
const K_RESTARTS    = 5
const MAX_ITER      = 200
```

- [ ] **Step 2: Run the baseline at Ecut=15**

```bash
julia --project=. -t auto workflow/joint_compression.jl 2>&1 | tee /tmp/joint_baseline.log
```
Expected: 30 min – 2 h, depending on Zygote backward pass speed. Final error: target < 8% to beat phase A.

- [ ] **Step 3: Write the bilan note**

Create `notes/05-bilan-phase-B-smoke.md` with:
- Configuration used
- Convergence curves (extracted from `summary.txt` / restart JSONs)
- Final best loss + H¹ relative error
- Comparison table: greedy@phase A (8.6%, 11 funcs) vs joint@phase B (?%, 15 funcs)
- Distribution of final log ζ across the basis (compute from `best.json`)
- Notes for future runs (Ecut=30, hyperparameter tuning)

- [ ] **Step 4: Commit**

```bash
git add workflow/joint_compression.jl notes/05-bilan-phase-B-smoke.md
git commit -m "Phase B: baseline run at Ecut=15 + first bilan"
```

---

## Self-review

1. **Spec coverage:** every section of the design spec maps to a task —
   §2 architecture → Task 4, §3 frame → Task 7-9, §4 paramétrisation
   → Task 3-4, §5 init → Task 5, §6 hyperparams → Task 11 config,
   §7 backend → Task 1+8, §8 outputs → Task 10, §9 tests → Task 11+12.

2. **Placeholder scan:** all code blocks are concrete, no TBD/TODO/FIXME.
   Each step has either a specific command or specific code to write.

3. **Type consistency:** `JointLayout`, `JointParams` (not used as struct in
   this plan, kept as flat vector), `LBFGSResult`, `JointOptimResult` are all
   referenced consistently. Method signatures `init_params`, `params_to_basis_functions`,
   `joint_loss`, `joint_inner_solve`, `run_lbfgs_once`, `run_joint_optim`
   match across tasks.

4. **Out-of-scope items deferred to v2:** auto-N, AdamW, Enzyme, GPU,
   penalty terms — none touched here, all flagged in the spec § "Hors scope".

---

## Execution

The plan is bite-sized, TDD-flavored where it matters, and uses small commits.
Each task ends with a commit so we can roll back to any milestone.
