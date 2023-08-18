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
function real_hamiltonian_bloc(R::AbstractVector{T}, basis_functions,
                               positions::AbstractVector) where {T<:Real}
    Nb = length(basis_functions)
    Hᴿ = zeros(T, Nb, Nb)
    Sᴿ = zeros(T, Nb, Nb)
    Z_carbon = 6 # TODO check

    # Assemble upper half
    for (i, Wᵢ) in enumerate(basis_functions)
        Wᵢᴿ = translate(Wᵢ, R)
        for (j, Wⱼ) in enumerate(basis_functions[i:end])
            # kinetic term + electron <-> ion + electron <-> electron integration
            Hᴿ[i,(j-1)+i] = integral(Wᵢᴿ, Wⱼ; type=:kinetic) +
                           -Z_carbon * sum(integral(Wᵢᴿ, Wⱼ, R_nuc; type=:atomic) for R_nuc in positions) +
                           integral(Wᵢᴿ, Wⱼ; type=:coulomb)
            # TODO: replace by zero when the upperbound is to low
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

# Use premade interpolation in ``Wannier/src/interp/fourier.jl``
