using JSON3

import Base.*, Base.-

mutable struct CompressedWannier{TR<:Real, TC<:Complex, T}
    # This data is always needed
    # Each basis_function is a name tuple with a list of coeffs
    # and corresponding SAGTOs.
    center          :: AbstractVector{TR}
    basis_functions :: Vector{BasisFunction}
    coefficients    :: AbstractVector{T}

    # The following is only needed during the compression process
    # but empty for tight binding or Hubbard.
    basis_supercell # PlaneWaveBasis type or nothing
    wannier         :: AbstractArray{TC}
    residual        :: AbstractArray{TC}
    error           :: TR
    error_norm
end
@inline function CompressedWannier(basis_supercell::PlaneWaveBasis{T}, Wn_pz::Vector{TC},
                                   center; error_norm=0) where {TC<:Complex, T}
    CompressedWannier(center, BasisFunction[], T[],
                      basis_supercell, Wn_pz, Wn_pz, NaN, error_norm # initial data for compression.
                      )
end

@inline function (Wc::CompressedWannier)(basis_supercell::PlaneWaveBasis)
    Wc.coefficients'*[Φ(basis_supercell) for Φ in Wc.basis_functions]
end

# Operations on Compressed Wannier
function translate(Wc::CompressedWannier, R::AbstractVector{T}) where T
    basis_functions = map(Wc.basis_functions) do Φ
        translated_SAGTOs = [GaussianPolynomial(X, X.center + R) for X in Φ.SAGTOs]
        BasisFunction(Φ.coeffs, translated_SAGTOs)
    end
    basis_functions = Vector{BasisFunction}(basis_functions)
    CompressedWannier(Wc.center + R, basis_functions, Wc.coefficients,
                      Wc.basis_supercell, Wc.wannier, Wc.residual, Wc.error, Wc.error_norm)
end

function (*)(λ::T, Wc::CompressedWannier) where T
    new_coeffs = λ .* Wc.coefficients
    CompressedWannier(Wc.center, Wc.basis_functions, new_coeffs,
                      Wc.basis_supercell, Wc.wannier, Wc.residual, Wc.error, Wc.error_norm)
end
(-)(Wc::CompressedWannier) = -1*Wc

"""
Apply rotation of angle θ to the given CompressedWannier around its center.
"""
function rotate(Wc::CompressedWannier, θ::T) where T
    α0 = Wc.center
    ℛθ(r) = rot(θ) * (r - α0) + α0
    basis_functions = map(Wc.basis_functions) do Φ
        rotated_SAGTOs = [GaussianPolynomial(X, ℛθ(X.center)) for X in Φ.SAGTOs]
        BasisFunction(Φ.coeffs, rotated_SAGTOs)
    end
    basis_functions = Vector{BasisFunction}(basis_functions)
    CompressedWannier(Wc.center, basis_functions, Wc.coefficients,
                      Wc.basis_supercell, Wc.wannier, Wc.residual, Wc.error, Wc.error_norm)
end

"""
TODO: doc
Store the compressed wannier. Residual is not stored as it can be recomputed easily.
"""
function store(Wc::CompressedWannier; file="compressed_wannier.json")
    data = Dict{String, Any}()
    data["basis_functions"] = []
    for Φ in Wc.basis_functions
        Φ_dict = Dict{String, Any}()
        Φ_dict["coeffs"] = Φ.coeffs
        SAGTO_as_array = []
        for X in Φ.SAGTOs
            push!(SAGTO_as_array, [pol_to_arrays(X.pol)..., X.center, X.spread])
        end
        Φ_dict["SAGTO"] = SAGTO_as_array
        push!(data["basis_functions"], Φ_dict)
    end
    data["center"] = Wc.center
    data["error"] = Wc.error
    data["error_norm"] = Wc.error_norm
    data["coefficients"] = Wc.coefficients

    # Store as JSON file
    open(io->JSON3.write(io, data, allow_inf=true), file, "w")
    nothing
end

function CompressedWannier(file)
    data = open(JSON3.read, file)
    basis_functions = BasisFunction[]
    T = eltype(data.coefficients)
    TC = complex(T)

    # Parse JSON3 arrays as BasisFunctions
    for Φ in data.basis_functions
        SAGTOs = []
        for X in Φ.SAGTO
            exps = Tuple.(X[1])
            coeffs = Vector(X[2])
            center = Vector(X[3])
            spread = X[4]
            push!(SAGTOs, GaussianPolynomial(exps, coeffs, center, spread))
        end
        coeffs = TC.(Φ.coeffs)
        push!(basis_functions, BasisFunction(coeffs, SAGTOs))
    end
    CompressedWannier(zeros(3), #T.(data.center),
                      basis_functions, T.(data.coefficients), nothing, TC[],
                      TC[], data.error, data.error_norm)
end
