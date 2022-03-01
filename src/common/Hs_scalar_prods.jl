using DFTK

"""
    Compute Hs scalar product in Fourier between two Block state
"""
function Hs_scalar_prod(basis::PlaneWaveBasis, kpt::Kpoint{T}, ψ1_k, ψ2_k; s=1) where {T<:Real}
    # Prefact for the Hs norm
    prefac_Hs = [(1 + norm(kpG_cart)^2)^s for kpG_cart in Gplusk_vectors_cart(basis, kpt)]
    dot(prefac_Hs .* ψ1_k, ψ2_k)
end

"""
    Compute the sum of Hs scalar product in Fourier between all given 
    Bloch states
"""
Hs_scalar_prod(basis::PlaneWaveBasis, ψ1, ψ2; ns=(1,1), s=1) = sum(
    Hs_scalar_prod(basis, kpt, ψ1[ik][:,ns[1]], ψ2[ik][:,ns[2]], s=s)
    for (ik, kpt) in enumerate(basis.kpoints) )

function Hs_norm(basis, ψ; n=1, s=1)
    square_norm = real(Hs_scalar_prod(basis, ψ, ψ, ns=(n,n), s=s))
    @assert (square_norm ≥ 0) "negative Hs_norm, try raising ζ_min"
    sqrt(square_norm)
end

function dist_Hs(basis::PlaneWaveBasis, ψ1, ψ2;
                 ns=(1,1), s=1)
    ψ1_minus_ψ2 = [ψ1k[:,ns[1]] .- ψ2k[:,ns[2]] for (ψ1k, ψ2k) in zip(ψ1, ψ2)]
    Hs_norm(basis, ψ1_minus_ψ2, n=1, s=s)
end

function Hs_overlap(basis::PlaneWaveBasis, Χs; s=1, n=1)
    num_aos = length(Χs)
    S = zeros(ComplexF64, num_aos, num_aos)
    # Run over all GTOs and compute overlaps
    for μ in 1:num_aos
        for ν in μ:num_aos
            S[μ,ν] = Hs_scalar_prod(basis, Χs[μ], Χs[ν], s=s, ns=(n,n))
            # Symmetrize
            (μ≠ν) && (S[ν,μ] = conj(S[μ,ν]))
        end
    end
    S
end

@doc raw"""
    Compute projection on a given AO basis ``Χs``. Format:
    ``Χs = [Χ1, ... , Χ_Nb]`` where each ``Χμ`` is a table ``[Χμ][ik]``
    of fourier coefficients of ``Χμ`` at frequency ``k+Gs``.
    Kwarg ``s`` is the choice of Hs norm for the projection.
"""
function Hs_projection_on_AO_basis(basis::PlaneWaveBasis, ψ, Χs;
                                   s=0)
    # Check that Ψ and all AOs have ben computed on the given basis kpoint grid
    num_k = length(basis.kpoints)
    @assert( (length(Χs[1])==num_k) && (length(ψ)==num_k) )

    # Compute the coefficients of the projection
    S = real.(Hs_overlap(basis, Χs, s=s))
    Χ = real.([Hs_scalar_prod(basis, ψ, Χμ, s=s) for Χμ in Χs])
    # Check for conditioning issues before inverting
    (cond(S) > 1e15) && (@warn "cond(S)>1e15"; return nothing, nothing)
    C_opti = (Symmetric(S)\Χ)
    # Assemble projection in Fourier and normalize
    sum( c .* Χμ for (c, Χμ) in zip(C_opti, Χs) ), C_opti
end
