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
