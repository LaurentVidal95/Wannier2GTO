"""
Cost function to minimize at each greedy algorithm iteration
"""
J(basis::PlaneWaveBasis, Φi, Res; s=1) = dist_Hs(basis, Φi, Res, s=s)^2

"""
J(ζ, λ(ζ), Res) where λ(ζ) ∈ argmin λ -> J(ζ, λ, Res)
"""
function J!(basis::PlaneWaveBasis, α0::Vector{T}, ζs::Vector{T}, r::T, θ::T, Res, Φ;
            s=1, ζ_min=1e-3,
            xy_orders=[0,2,3], z_orders=[1,3,5],
            return_all=false,
            ) where {T<:Real}

    # Avoid NaN when ζ is to small
    (findmin(ζs)[1] < ζ_min) && (return 1e10)

    # Compute SAGTOs for given center and respective spreads
    α(r, θ) = polar_to_cartesian_coords(α0, r, θ)
    Χs = SAGTOs_graphene_pz_wannier(basis, α(r, θ), ζs, xy_orders, z_orders)
    
    # Compute λ_opti s.t. Φ = ∑λ_opti_μ*Χμ ∈ argmin λ -> J(ζ, λ, Res)
    Γ = real.(ThreadsX.map(Χ -> Hs_scalar_prod(basis, Res, Χ, s=s), Χs))
    S = real.(Hs_overlap(basis, Χs, s=s))
    λ_opti = inv(Symmetric(S))*Γ; len_λ_opti = copy(length(λ_opti));
    α_Φ = [α(r, θ)]

    # Add AOs at jx and j^2x if needed and update centers and coefficients
    D3_sym = iszero(r);
    if !(D3_sym)
        append!(Χs,
                SAGTOs_graphene_pz_wannier(basis, α(r, (2π/3)+θ), ζs, xy_orders, z_orders))
        append!(Χs,
                SAGTOs_graphene_pz_wannier(basis, α(r, (4π/3)+θ), ζs, xy_orders, z_orders))
        α_Φ = [α(r, θ), α(r, (2π/3)+θ), α(r, (4π/3)+θ)]
        λ_opti = vcat(λ_opti, λ_opti, λ_opti)
    end

    # Compute J with optimal linear combination of Χs
    Φ_opti =  sum(λ .* Χ for (λ, Χ) in zip(λ_opti, Χs))
    # Actualize "in place" storage for exterior computations
    for (μ, Φμ_opti) in enumerate(Φ_opti)
        Φ[μ] .= Φμ_opti
    end
    # Return only J (for Optim routine) or J and coefficient (for storage)
    return_all && (return J(basis, Φ_opti, Res, s=s), λ_opti[1:len_λ_opti], α_Φ)
    J(basis, Φ_opti, Res, s=s)
end

function test_J1(basis::PlaneWaveBasis, α::Vector{T}, ζs::Vector{T}, Res;
                 xy_orders=[0,2,3], z_orders=[1,3,5],
                 return_all=false,
                 ) where {T<:Real}
    
    # Compute all SAGTOs in fourier for given centers and spreads
    SAGTOs = SAGTOs_graphene_pz_wannier(basis, α(r, θ), ζs, xy_orders, z_orders)
    
    # First compute optimal Λα
    Γ = real.(ThreadsX.map(Χ -> Hs_scalar_prod(basis, Res, Χ, s=s), SAGTOs))
    S = real.(Hs_overlap(basis, SAGTOs, s=s))
    λζ = inv(Symmetric(S))*Γ
    α_Φ = α(r, θ)
    
    # Compute J with optimal Λα
    Φ_opti =  sum( λ .* Χ for (λ, Χ) in zip(λζ, SAGTOs) )    
end

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
