function apply_term(term::Symbol, W₁::CompressedWannier, W₂::CompressedWannier, scfres=nothing)
    KS_potential_term = (term==:PBE_local_potential)
    if isnothing(scfres) && KS_potential_term
        error("The KS potential is required to compute the :PBE_potential term contribution")
    end
    
    if !KS_potential_term
        return integral(f₁, f₂; type=term)
    else
        V = extract_local_KS_potential(scfres)
        return potential_scalar_prod(scfres.basis, V, W₁, W₂)
    end
    error("Well, that was inexpected")
end

function hamiltonian_scalar_prod(BG::BilayerGraphene, W₁::CompressedWannier,
                                 W₂::CompressedWannier; terms=[:kinetic, :PBE_local_potential])
    sum(apply_term(term, W₁, W₂, BG.scfres) for term in terms)
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{T},
                               X₁::GaussianPolynomial, X₂::GaussianPolynomial) where {T<:Complex}
    Ω = basis.model.unit_cell_volume
    X₃ = X₁*X₂
    X₃_four = ℱ.(Ref(X₃), G_vectors_cart(basis))
    (1/√Ω)*dot(X₃_four,V_KS_four)
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{TC},
                               Φ₁::BasisFunction, Φ₂::BasisFunction) where {TC<:Complex}
    output = zero(TC)
    for (λ₁, X₁) in zip(Φ₁.coeffs, Φ₁.SAGTOs)
        for (λ₂, X₂) in zip(Φ₂.coeffs, Φ₂.SAGTOs)
            @show output += λ₁*λ₂*potential_scalar_prod(basis, V_KS_four, X₁, X₂)
        end
    end
    output
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{TC},
                               Wc₁::CompressedWannier, Wc₂::CompressedWannier) where {TC<:Complex}
    output = zero(TC)
    for (λ₁, Φ₁) in zip(Wc₁.coefficients, Wc₁.basis_functions)
        for (λ₂, Φ₂) in zip(Wc₂.coefficients, Wc₂.basis_functions)
            output += λ₁*λ₂*potential_scalar_prod(basis, V_KS_four, Φ₁, Φ₂)
        end
    end
    output
end
