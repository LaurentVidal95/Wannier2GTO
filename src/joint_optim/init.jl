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
