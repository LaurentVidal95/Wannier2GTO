using ProgressMeter

"""
   real_hamiltonian_bloc(R::AsbtractVector, W₁::CompressedWannier, W₂::CompressedWannier,
                          positions::AbstractVector)

Construct the Nb×Nb tight-binding hamiltonian ``H_{\\bm{R}}`` for given vector ``R``
of the lattice. The Nb basis functions are pz-like Wannier functions of graphene with 
``Nb = 2n`` for n--layer graphene. 

The hamiltonian is given for all ``1 ≤ i,j ≤ Nb`` by:

```math
[H_{\\bm{R}}]_{ij} = \\langle Wᵢ(⋅ - \\bm{R})\\,|\\, H \\,|\\, Wⱼ \\rangle
```
"""
function real_hamiltonian(R::AbstractVector{T}, TB::TightBindingModel) where {T<:Real}

    Nb = length(TB.basis_functions)
    Hᴿ = zeros(T, Nb, Nb)
    Sᴿ = zeros(T, Nb, Nb)
    Z_carbon = 6 # TODO check

    # Assemble upper half
    for (i, Wᵢ) in enumerate(TB.basis_functions)
        Wᵢᴿ = translate(Wᵢ, -R)
        for (j, Wⱼ) in enumerate(TB.basis_functions[i:end])
            Hᴿ[i,(j-1)+i] = hamiltonian_scalar_prod(TB, Wᵢᴿ, Wⱼ)
            Sᴿ[i,(j-1+i)] = integral(Wᵢᴿ, Wⱼ; type=:overlap)
        end
    end

    # Deduce lower half by symmetry
    for i in 1:Nb
        for j in i+1:Nb
            Hᴿ[j,i] = Hᴿ[i,j]
            Sᴿ[j,i] = Sᴿ[i,j]
        end
    end

    Symmetric(Hᴿ), Symmetric(Sᴿ)
end

# Time consuming part
function real_hamiltonian(TB::TightBindingModel, atomic_positions::AbstractVector)
    Rs_cart = R_vectors_cart(TB)    
    T = eltype(Rs_cart[1])
    NR = length(Rs_cart)
    Nb = length(TB.basis_functions)

    ℍᴿ = zeros(T, Nb, Nb, NR)
    𝕊ᴿ = zeros(T, Nb, Nb, NR)

    progress = Progress(length(Rs_cart), desc="Computing all blocs")
    for (iR, R) in enumerate(Rs_cart)
        Hᴿ, Sᴿ = real_hamiltonian(R, TB.basis_functions, atomic_positions)
        ℍᴿ[:,:,iR] .= Hᴿ
        𝕊ᴿ[:,:,iR] .= Sᴿ
        next!(progress)
    end
    ℍᴿ, 𝕊ᴿ
end

function bloch_transform(kpoints::Vector{Vector{T}}, TB::TightBindingModel,
                         𝕆ᴿ::Array{T,3}) where {T<:Real}
    Nb = size(𝕆ᴿ)[1]
    Nk = length(kpoints)
    @assert size(𝕆ᴿ)[1] == size(𝕆ᴿ)[2]

    R_vectors = TB.R_vectors
    R_basis = [W.center for W in TB.basis_functions]

    𝕆ᵏ = zeros(Complex{T}, Nb, Nb, Nk)
    for (ik, k) in enumerate(kpoints)
        for (iR, R) in enumerate(R_vectors)
            prefac = [cis(-2π*dot(R .+ Rᵢ .-Rⱼ, k)) for Rᵢ in R_basis, Rⱼ in R_basis]
            𝕆ᵏ[:,:,ik] .+= prefac .* @view 𝕆ᴿ[:,:,iR]
        end
    end
    𝕆ᵏ
end
