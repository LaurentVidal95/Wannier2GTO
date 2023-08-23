struct TightBindingModel{T<:Real}
    # unit cell data
    system::BilayerGraphene{T}
    basis_functions::AbstractVector
    hamiltonian_terms::Vector{Symbol}
    
    # Neighbours information
    R_vectors::AbstractArray
    atomic_positions::AbstractVector
    tol::T
end

# Generic function that works for both graphene and stacked bilayer graphene.
geometry(TB::TightBindingModel) = (TB.system.lattice, TB.system.atoms, TB.system.positions)
R_vectors_cart(TB::TightBindingModel) = map(R->TB.system.lattice*R, TB.R_vectors)

function TightBindingModel(BG::BilayerGraphene;
                           tol=1e-8,
                           basis_function_file="compressed_wannier_55_H1.json",
                           terms = [:kinetic, :atomic, :coulomb])

    # Extract CompressedWannier center at zero from data file
    path_dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2])
    datadir = joinpath(path_dir, "data")
    W = CompressedWannier(joinpath(datadir, basis_function_file))

    periodic = length(BG.positions) == 4
    if !(periodic)
        error("Non periodic system are still to be implemented")
    end        

    # Construct  W₁ and W₂, the two pz-like wannier functions of monolayer graphene from
    # the given function W, centered at zero, in CompressedWannier format.
    # If Rᵢ is the position of the i-th carbon atom in the unit cell, one has:
    # - W₁ = τ(R₁)(W)
    # - W₂ = τ(R₂) ∘ ℛ(π/3)(W)
    # where τ(R) is the translation of vector R and ℛ(π/3) the rotation of angle π/3
    # around the wannier center.
    R₁, R₂, R₃, R₄ = Ref(BG.lattice) .* BG.positions
    W₁ = translate(W, R₁)
    W₂ = translate(rotate(W, π/3), R₂)
    # Add a second layer
    W₃ = translate(W₁, R₃-R₁)
    W₄ = translate(W₂, R₄-R₂)    

    R_vectors = get_R_vectors(W, BG; tol)
    positions = [W.center for W in [W₁, W₂, W₃, W₄]] # only for periodic systems
    
    TightBindingModel(BG, [W₁, W₂, W₃, W₄], terms, R_vectors, positions, tol)
end

function get_R_vectors(W::CompressedWannier, BG::BilayerGraphene;
                       tol=1e-8,
                       N_R_vectors_max=500)
    # R vectors are of the form R(i,j) = BG.lattice*[i,j,0], for (i,j)∈ℤ.
    # Computes the maximum i so that ⟨W, τ(R(i,i))(W)⟩_L² < tol.
    # Just renaming for clarity
    N_R_vectors(i) = (2i+1)^2
    overlap(W, i) = integral(W, translate(W, BG.lattice*[i,i,0]))

    i_max = zero(Int64)
    while (abs(overlap(W, i_max)) > tol) && (N_R_vectors(i_max) < N_R_vectors_max)
        i_max += 1
        (i_max) > 20 && error("Residual error")
    end

    [[i,j,0] for i in -i_max:i_max, j in -i_max:i_max]
end
