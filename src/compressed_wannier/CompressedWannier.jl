mutable struct CompressedWannier{TR<:Real, TC<:Complex}
    # Original Wannier function to compress
    basis_SC :: PlaneWaveBasis{TR}
    wannier :: AbstractArray{TC}
    center :: Vector{TR}
    # Each MO is a name tuple with a list of coeffs and corresponding SAGTOs
    MOs
    # Optim related
    error
    error_norm
end

function CompressedWannier(basis_SC::PlaneWaveBasis, Wn_pz::AbstractArray, center;
                            error_norm=:L2)
    # Initiate other args
    TR = eltype(basis_SC)
    MOs = [(;coeffs=[], SAGTOs=[])]
    error = NaN
    CompressedWannier(basis_SC, Wn_pz, center, MOs, error, error_norm)
end

function pol_to_arrays(pol::Polynomial)
    exps = Tuple.(eachcol(StaticPolynomials.exponents(pol)))
    coeffs = pol.coefficients
    exps, coeffs
end

function store(Wc::CompressedWannier; file="compressed_wannier.json")
    data = Dict{String, any}
    data["center"] = Wc.center
    data["MOs"] = []
    for Φ in Wc.MOs
        Φ_dict = Dict{String, any}
        Φ_dict["coeffs"] = Φ.coeffs
        SAGTO_as_array = []
        for X in Φ.SAGTOs
            push!(SAGTO_as_array, [pol_to_arrays(X.pol)..., X.center, X.spread])
        end
        Φ_dict["SAGTO"] = SAGTO_as_array
        push!(data["MOs"], Φ_dict)
    end
    data["error"] = Wc.error
    data["error_norm"] = Wc.error_norm
        
    # Store as JSON file
    open(io->JSON3.write(io, data, allow_inf=true), file, "w")
end
