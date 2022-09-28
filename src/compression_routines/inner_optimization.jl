"""
J(ζ, λ(ζ), res) where λ(ζ) ∈ argmin λ -> J(ζ, λ, res)
"""
function inner_optimization(info,
                            # Spread and center (in polar coords) of SAGTOs
                            α0::Vector, r::T, θ::T, ζs::Vector{T};
                            # Returns leat square problem coefficients
                            in_linesearch=false,
                            ) where {T<:Real}

    # Extract info
    s = info.s
    basis_SC = info.basis_SC
    res = info.res

    # Avoid NaN when ζ is to small
    (findmin(ζs)[1] < info.ζ_min) && (return Inf)

    # Compute SAGTO basis for given polynomial part orders, center and exponents
    xy_orders, z_orders = info.pol_orders
    α(r, θ) = polar_to_cartesian_coords(α0, r, θ)
    Χs = SAGTOs_basis(basis_SC, α(r, θ), ζs, xy_orders, z_orders)
    
    # Compute λ_opti s.t. Φ = ∑λ_opti_μ*Χμ ∈ argmin λ -> J(ζ, λ, res)
    Γ = ThreadsX.map(Χ -> Hs_scalar_prod(basis_SC, info.res, Χ; s), Χs)
    S = Hs_overlap(basis_SC, Χs; s)
    λ_opti = (Symmetric(S)\Γ)
    
    αΦ = [α(r, θ)]
    # Add AOs at jx and j^2x if needed and update centers and coefficients
    if !(info.D3_sym)
        append!(Χs, SAGTOs_basis(basis_SC, α(r, (2π/3)+θ), ζs, xy_orders, z_orders))
        append!(Χs, SAGTOs_basis(basis_SC, α(r, (4π/3)+θ), ζs, xy_orders, z_orders))
        αΦ = [α(r, θ), α(r, (2π/3)+θ), α(r, (4π/3)+θ)]
        λ_opti = vcat(λ_opti, λ_opti, λ_opti)
    end

    # Compute J with optimal linear combination of Χs
    Φ_opti =  sum(λ .* Χ for (λ, Χ) in zip(λ_opti, Χs))

    # Actualize info for outer J! computation
    info = merge(info, (; Φ=Φ_opti))
    
    
    # Return only J (for Optim routine) or J and coefficient (for storage)
    !(in_linesearch) && (return info, dist_Hs(basis_SC, Φ_opti, res, s=s), λ_opti, αΦ)
    dist_Hs(basis_SC, Φ_opti, res, s=s)
end
