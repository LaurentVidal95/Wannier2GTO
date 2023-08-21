struct TightBindingModel{T<:Real}
    # unit cell data
    system::BilayerGraphene{T}
    basis_functions::AbstractVector
    # Neighbours information
    R_vectors
    tol::T
end

# Generic function that works for both graphene and stacked bilayer graphene.
geometry(TB::TightBindingModel) = (TB.model.lattice, TB.model.atoms, TB.model.positions)
R_vectors_cart(TB::TightBindingModel) = map(R->TB.system.lattice*R, TB.R_vectors)

function TightBindingModel(BG::BilayerGraphene, tol;
                           basis_function_file="compressed_wannier_55_H1.json")
    # Extract CompressedWannier center at zero from data file
    path_dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2])
    datadir = joinpath(path_dir, "data")
    W = CompressedWannier(joinpath(datadir, basis_function_file))

    # TODO: modify
    # Compute W₂, the second basis function in the unit cell of graphene.
    # If Rᵢ is the position of the i-th carbon atom in the unit cell, one has:
    # W₂ = τ(R₂-R₁) ∘ ℛ(π/3)(W₁)
    # where τ(R) is the translation of vector R and ℛ(π/3) the rotation of angle π/3
    # around the wannier center.
    R₁, R₂ = Ref(BG.lattice) .* BG.positions[1:2]
    W₁ = translate(W, R₁)
    W₂ = translate(rotate(W, π/3), R₂)
    
    periodic = length(BG.positions) == 4
    if !(periodic)
        error("Non periodic system to be implemented")
    end        
    
    TightBindingModel(BG, [W₁, W₂], [zeros(2)], tol)
end
