using JSON3

mutable struct CompressedWannier{TR<:Real, TC<:Complex}
    # Data on the original Wannier function to compress
    basis_supercell :: PlaneWaveBasis
    wannier         :: AbstractArray{TC}
    center          :: Vector{TR}
    # Each MO is a name tuple with a list of coeffs and corresponding SAGTOs
    basis_functions :: Vector{BasisFunction}
    coefficients    :: Vector{TC}
    # Optim related
    residual        ::AbstractArray{TC}
    error           ::TR
    error_norm
end
@inline function CompressedWannier(basis_SC::PlaneWaveBasis, Wn_pz::Vector{TC},
                                   center; error_norm=0) where {TC<:Complex}
    CompressedWannier(basis_SC, Wn_pz, center, BasisFunction[], TC[],  Wn_pz, NaN, error_norm)
end

function pol_to_arrays(pol::Polynomial)
    exps = Tuple.(eachcol(StaticPolynomials.exponents(pol)))
    coeffs = pol.coefficients
    exps, coeffs
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
    data["error"] = Wc.error
    data["error_norm"] = Wc.error_norm
    data["coefficients"] = Wc.coefficients

    # Store as JSON file
    open(io->JSON3.write(io, data, allow_inf=true), file, "w")
    nothing
end

# function CompressedWannier(basis_supercell, wannier, center, file)
#     # TODO: complete the function ! In particular have to parse JSON arrays as BasisFunction
#     data = open(JSON3.read, file)
#     basis_functions = [BasisFunction(ComplexF64.(Φ.coeffs), Φ.SAGTO) for Φ in test.basis_functions]
# end
