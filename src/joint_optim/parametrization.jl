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
