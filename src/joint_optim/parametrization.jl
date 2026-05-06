# JointLayout, flatten/unflatten, sigmoid encoding, params → BasisFunctions

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
    # Non-mutating form: avoids push! so Zygote can differentiate through flat.
    # Centered functions
    xy_c, z_c = select_orders(max_xy_order, max_z_order, true)
    centered = [begin
        u, λ = get_centered(flat, layout, i)
        ζ = exp(log_ζ_from_u(u, log_ζ_min, log_ζ_max))
        SAGTOs = SAGTO_basis(zeros(eltype(flat), 3), ζ, xy_c, z_c)
        BasisFunction(collect(λ), SAGTOs)
    end for i in 1:layout.N_centered]

    # π-bond functions, D3 enforced by construction (3 rotational copies)
    xy_p, z_p = select_orders(max_xy_order, max_z_order, false)
    pibond = [begin
        r, u, λ = get_pibond(flat, layout, j)
        ζ = exp(log_ζ_from_u(u, log_ζ_min, log_ζ_max))
        center = r .* π_bond_unit
        SAGTOs = SAGTO_basis(collect(center), ζ, xy_p, z_p)
        Φ = BasisFunction(collect(λ), SAGTOs)
        enforce_D3_symmetry(Φ)
    end for j in 1:layout.N_pibond]

    vcat(centered, pibond)
end
