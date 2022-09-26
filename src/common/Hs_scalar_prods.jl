"""
Computes the Hs scalar product of two vectors containing all (k+G) Fourier coefficients.
In DFTK conventions, these vectors are Bloch wave of the supercell Gamma point.
Conversion from cell to supercell is hence needed
"""
function Hs_scalar_prod(basis_SC::PlaneWaveBasis{T}, ψ1, ψ2; s=0) where {T<:Real}
    prefac_Hs = [(1+norm(kpG)^2)^s
                 for kpG in G_vectors_cart(basis_SC, only(basis_SC.kpoints))]
    real(dot(prefac_Hs .* ψ1, ψ2))
end

function Hs_norm(basis_SC::PlaneWaveBasis, ψ; s=0)
    square_norm = Hs_scalar_prod(basis_SC, ψ, ψ, s=s)
    @assert (square_norm ≥ 0) "negative Hs_norm, try raising ζ_min"
    sqrt(square_norm)
end

dist_Hs(basis_SC::PlaneWaveBasis, ψ1, ψ2; s=0) = Hs_norm(basis_SC, ψ1 .- ψ2, s=s)

function Hs_overlap(basis_SC::PlaneWaveBasis, Χs_fourier; s=0)
    num_aos = length(Χs_fourier)
    S = zeros(Float64, num_aos, num_aos)
    # Run over all GTOs and compute overlaps
    for μ in 1:num_aos
        for ν in μ:num_aos
            S[μ,ν] = Hs_scalar_prod(basis_SC, Χs_fourier[μ], Χs_fourier[ν], s=s)
            # Symmetrize
            (μ≠ν) && (S[ν,μ] = S[μ,ν])
        end
    end
    S
end

@doc raw"""
    Compute projection on a given AO basis ``Χs``. Format:
    ``Χs = [Χ1, ... , Χ_Nb]`` where each ``Χμ`` is a vector ``Χμ[k+G]``
    of fourier coefficients of ``Χμ`` at frequency ``k+Gs``.
    Kwarg ``s`` is the choice of Hs norm for the projection.
"""
function Hs_projection_on_AO_basis(basis_SC::PlaneWaveBasis, ψ, Χs; s=0)
    # Check that Ψ and all AOs have been converted to supercell conventions
    num_kpG = length(G_vectors(basis_SC, only(basis_SC.kpoints)))
    @assert( (length(ψ)==num_kpG) && (length(Χs[1])==num_kpG) )

    # Compute the coefficients of the projection
    S = Hs_overlap(basis_SC, Χs, s=s)
    Χ = [Hs_scalar_prod(basis_SC, ψ, Χμ, s=s) for Χμ in Χs]
    # Check for conditioning issues before inverting
    (cond(S) > 1e15) && (@warn "cond(S)>1e15"; return nothing, nothing)
    C_opti = (Symmetric(S)\Χ)
    # Assemble projection in Fourier and normalize
    sum( c .* Χμ for (c, Χμ) in zip(C_opti, Χs) ), C_opti
end
