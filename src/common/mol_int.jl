mutable struct compressed_wannier{T <: Real}
    # Compression parameters
    max_xy_order :: Int64
    max_z_order :: Int64
    precision :: Float64
    # MOs as LC-SAGTOs tables
    ζs :: Vector{Vector{T}}
    αs :: Vector{Matrix{T}}
    λs :: Vector{Vector{T}}
    # Projection coefficients
    proj_coeffs :: Vector{T}
end

## TODO better files...
function compressed_wannier(file, basis::PlaneWaveBasis, w, max_xy_order, max_z_order;
                            s=0)
    ζs, αs, λs = read_compression_file(file)
    # Estimate proj coeffs and precision
    Φ_tmp = construct_compression_MO_basis(basis, file, max_xy_order, max_z_order)
    w_proj, proj_coeffs = projection_on_GTOBasis(basis, w, Φ_tmp, s=s)
    res_norm = norm(w .- w_proj)/norm(w)
    compressed_wannier(max_xy_order, max_z_order, res_norm, ζs, αs, λs, proj_coeffs)
end

"""
   Each SAGTO Χ is writen as a sum of Gaussian-polynomial functions with same center
   and spread : Χ = ∑_i λi gi(X - α)

   Hence ⟨Χ1 | H | Χ2⟩ = ∑_ij λi*λj ⟨gi(.-α1)| H |gj(.-α2)⟩

   Iα gives for each center the powers of the polynomial part of all g.
   λα gives all λs in the same order than Iα.
"""
function SAGTOs_Ham_scalar_product(ζ1, α1, Iα_1, λα_1, ζ2, α2, Iα_2, λα_2;
                              terms=["overlap"])
    @assert !isempty(intersect(terms,
                ["overlap","kinetic","coulomb"])) "Hamiltoninan term not handled"
    sp = zero(Float64)
    for ((nx1, ny1, nz1), λ1) in zip(Iα_1, λα_1)
        for ((nx2, ny2, nz2), λ2) in zip(Iα_2, λα_2)
            ("kinetic"∈ terms) && (sp += kinetic(ζ1, α1, nx1, ny1, nz1, ζ2, α2, nx2, ny2, nz2))
            ("coulomb"∈ terms) && (sp += c_coulomb_c(ζ1, α1, nx1, ny1, nz1, ζ2, α2, nx2, ny2, nz2))
            ("overlap"∈ terms) && (sp += overlap(ζ1, α1, nx1, ny1, nz1, ζ2, α2, nx2, ny2, nz2))
            sp *= λ1*λ2
        end
    end
    sp
end

"""
    Each MO Φ is a linear combination of SAGTOs with different spreads and centers.
    Φ = ∑_i λαi Χi(. - αi)
    Tables αs contains all centers of the Χi, ζs and λs corresponding spread and coefficient.
"""
function LC_SAGTOs_Ham_scalar_product(params_μ, ζs_μ, αs_μ, λs_μ, params_ν, ζs_ν, αs_ν, λs_ν;
                            terms=["kinetic", "coulomb"])
    # Compute all scalar products
    sp = zero(Float64)    
    # Loop over centers. One center = one SAGTO
    for αi in eachcol(αs_μ)
        for αj in eachcol(αs_ν)
            T = eltype(αi)
            # Extract corresponding spread, coefficient and parameters
            for (ζi, λi, (Iα_i, λα_i)) in zip(ζs_μ, λs_μ, params_μ)
                for (ζj, λj, (Iα_j, λα_j)) in zip(ζs_ν, λs_ν, params_ν)
                    # Compute λi*λj*⟨Χi|H|Χj⟩
                    sp += λi*λj*SAGTOs_Ham_scalar_product(ζi, Vector{T}(αi), Iα_i, λα_i,
                                                          ζj, Vector{T}(αj), Iα_j, λα_j,
                                                          terms=terms)
                end
            end
        end
    end
    sp
end
# LC_SAGTOs_Ham_norm(max_xy_order, max_z_order, ζs_μ, αs_μ, λs_μ) =
#     √(LC_SAGTOs_Ham_scalar_product(max_xy_order, max_z_order, ζs_μ, αs_μ, λs_μ, ζs_μ, αs_μ, λs_μ,
#                        terms=["overlap"]))

function Ham_scalar_product(w1::compressed_wannier{T}, w2::compressed_wannier{T};
                            terms=["kinetic","coulomb"]) where {T<:Real}
    
    Nb_1 = length(w1.ζs); Nb_2 = length(w2.ζs);
    Hc = zeros(Float64, Nb_1, Nb_2)

    # Assemble [Hc]μν = ⟨Φμ | H | Φν⟩
    for μ in 1:Nb_1
        # Chose parameters according to the current MO symmetry
        D3_sym_μ = (size(w1.αs[μ],2) > 1)
        params_μ = SAGTO_parameters(w1.max_xy_order, w1.max_z_order, D3_sym_μ)
        for ν in 1:Nb_2
            D3_sym_ν = (size(w2.αs[ν],2) > 1);
            params_ν = SAGTO_parameters(w2.max_xy_order, w2.max_z_order, D3_sym_ν)
            Hc[μ,ν] = LC_SAGTOs_Ham_scalar_product(params_μ, w1.ζs[μ], w1.αs[μ], w1.λs[μ],
                                                   params_ν, w2.ζs[ν], w2.αs[ν], w2.λs[ν],
                                                   terms=terms)
            println(Hc[μ,ν])
            (μ≠ν) && (Hc[ν,μ] = Hc[μ,ν])
        end
    end
    # Compute scalar product
    w1.proj_coeffs'Hc*w2.proj_coeffs
end
