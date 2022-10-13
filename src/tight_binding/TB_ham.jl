function τ(basis_SC::PlaneWaveBasis{T}, ψ::Vector{ComplexF64}, R::Vector{T}) where {T<:Real}
    Γ = only(basis_SC.kpoints)
    # check that ψ is a supercell vector
    @assert(length(Γ.G_vectors) == length(ψ))
    # Shift in fourier [τ_R(ψ)]_G = [ψ]_G * exp(iG⋅R)
    ψ .*= -cis.(dot.(G_vectors_cart(basis_SC, Γ), Ref(R)))
end

# function tight_binding_ham(basis_SC::PlaneWaveBasis{T}, w1, w2, neighbours, ρ)
#     ham = Hamiltonian(basis_SC; ρ)
#     H_Γ = only(ham.blocks) 
    
