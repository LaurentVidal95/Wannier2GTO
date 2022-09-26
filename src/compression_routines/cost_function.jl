"""
J(ζ, λ(ζ), Res) where λ(ζ) ∈ argmin λ -> J(ζ, λ, Res)
"""
function J(info,
            # Spread and center (in polar coords) of SAGTOs
            α0::Vector{T}, r::T, θ::T, ζs::Vector{T};
            # Returns leat square problem coefficients
            in_linesearch=false,
            ) where {T<:Real}

    # Extract info
    s = info.s
    basis_SC = info.basis_SC
    Res = info.Res

    # Avoid NaN when ζ is to small
    (findmin(ζs)[1] < info.ζ_min) && (return Inf)

    # Compute SAGTO basis for given polynomial part orders, center and exponents
    xy_orders, z_orders = info.pol_orders
    α(r, θ) = polar_to_cartesian_coords(α0, r, θ)
    Χs = SAGTOs_basis(basis_SC, α(r, θ), ζs, xy_orders, z_orders)

    # Compute λ_opti s.t. Φ = ∑λ_opti_μ*Χμ ∈ argmin λ -> J(ζ, λ, Res)
    Γ = ThreadsX.map(Χ -> Hs_scalar_prod(basis_SC, info.Res, Χ; s), Χs)
    S = Hs_overlap(basis_SC, Χs; s)
    λ_opti = (Symmetric(S)\Γ)

    αΦ = [α(r, θ)]
    # Add AOs at jx and j^2x if needed and update centers and coefficients
    D3_sym = iszero(r);
    if !(D3_sym)
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
    !(in_linesearch) && (return info, dist_Hs(basis_SC, Φ_opti, Res, s=s), λ_opti, αΦ)
    dist_Hs(basis_SC, Φ_opti, Res, s=s)
end

# function test_J1(basis::PlaneWaveBasis, α::Vector{T}, ζs::Vector{T}, Res;
#                  xy_orders=[0,2,3], z_orders=[1,3,5],
#                  in_linesearch=false,
#                  ) where {T<:Real}
    
#     # Compute all SAGTOs in fourier for given centers and spreads
#     SAGTOs = SAGTOs_graphene_pz_wannier(basis, α(r, θ), ζs, xy_orders, z_orders)
    
#     # First compute optimal Λα
#     Γ = real.(ThreadsX.map(Χ -> Hs_scalar_prod(basis, Res, Χ, s=s), SAGTOs))
#     S = real.(Hs_overlap(basis, SAGTOs, s=s))
#     λζ = inv(Symmetric(S))*Γ
#     α_Φ = α(r, θ)
    
#     # Compute J with optimal Λα
#     Φ_opti =  sum( λ .* Χ for (λ, Χ) in zip(λζ, SAGTOs) )    
# end

# function ∇J(basis::PlaneWaveBasis, α::Vector{T}, ζs::Vector{T}, λζ::Vector{T},
#               SAGTOs, Res,
#               xy_orders, z_orders;
#               s=0) where {T<:Real}

#     ∇SAGTOs = ∇ζ_SAGTOs_graphene_pz_wannier(basis, α, ζs, xy_orders, z_orders)

#     grad_J = zeros(ComplexF64, length(ζs))
#     # Compute all scalar product
#     for (i,∇Χi) in enumerate(∇SAGTOs)
#         grad_J[i] += -2*λζ[i]*Hs_scalar_prod(basis, Res, ∇Χi, s=s)
#         for (j, Χj) in enumerate(SAGTOs)
#             δij = 1; (i==j) && (δij+=1)
#             grad_J[i] += (λζ[i]*λζ[j]*δij) * Hs_scalar_prod(basis, ∇Χi, Χj, s=s)
#         end
#     end
#     grad_J
# end
