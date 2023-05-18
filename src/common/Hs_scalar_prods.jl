"""
Computes the Hs scalar product of two vectors containing all (k+G) Fourier coefficients.
In DFTK conventions, these vectors are Bloch wave of the supercell Gamma point.
Conversion from cell to supercell is hence needed beforehand.
"""
function Hs_dot(basis_SC::PlaneWaveBasis{T}, ψ1, ψ2; s=0) where {T<:Real}
    # Handle L² norm case
    s==0 && return dot(ψ1, ψ2)
    # s ≥ 1
    prefac_Hs = [(1+norm(kpG)^2)^s
                 for kpG in G_vectors_cart(basis_SC, only(basis_SC.kpoints))]
    real(dot(prefac_Hs .* ψ1, ψ2))
end

function Hs_norm(basis_SC::PlaneWaveBasis, ψ; s=0)
    square_norm = HS_dot(basis_SC, ψ, ψ, s=s)
    @assert (square_norm ≥ 0) "negative Hs_norm, try raising minimal spread ζ_min"
    sqrt(square_norm)
end

function Hs_overlap(basis_SC::PlaneWaveBasis, Χs_fourier; s=0)
    num_aos = length(Χs_fourier)
    # Run over all GTOs and compute overlaps
    [Hs_dot(basis_SC, Χs_fourier[μ], Χs_fourier[ν]; s) for μ in 1:num_aos, ν in 1:num_aos]
end

@doc raw"""
Compute projection on a given AO basis ``Χs``. Format:
``Χs = [Χ1, ... , Χ_Nb]`` where each ``Χμ`` is a vector ``Χμ[k+G]``
of fourier coefficients of ``Χμ`` at frequency ``k+Gs``.
Kwarg ``s`` is the choice of Hs norm for the projection.
"""
function project_on_AO_basis(basis_SC::PlaneWaveBasis, ψ,
                             Χs::Vector{Vector{ComplexF64}}; s=0)
    Xs = normalize.(Xs) # renormalize AOs to avoid super big or small coeffs

    # Check that Ψ and all AOs have been converted to supercell conventions
    num_kpG = length(G_vectors(basis_SC, only(basis_SC.kpoints)))
    @assert( (length(ψ)==num_kpG) && (length(Χs[1])==num_kpG) )

    # Compute the coefficients of the projection
    S = Hs_overlap(basis_SC, Χs; s)
    Χ = [Hs_dot(basis_SC, ψ, Χμ; s) for Χμ in Χs]
    # Check for conditioning issues before inverting and stop if conditioning is to high
    (cond(S) > 1e8) && (error("cond(S)>1e8"))
    C_opti = (Symmetric(S)\Χ)
    # Assemble projection in Fourier
    (; func=sum( c .* Χμ for (c, Χμ) in zip(C_opti, Χs) ), coeffs=C_opti)
end
