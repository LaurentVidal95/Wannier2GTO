"""
Computes the Hs scalar product of two vectors containing all (k+G) Fourier coefficients.
In DFTK conventions, these vectors are Bloch wave of the supercell Gamma point.
Conversion from cell to supercell is hence needed beforehand.
"""
function Hs_dot(basis_SC::PlaneWaveBasis, ψ1::AbstractVector{T1},
                ψ2::AbstractVector{T2}; s=0) where {T1, T2<:Complex}
    # Handle trivial case (L² norm)
    iszero(s) && return dot(ψ1, ψ2)
    # Hˢ norm (s ≠ 0)
    prefac_Hs = [(1+norm(Gpk)^2)^s
                 for Gpk in G_vectors_cart(basis_SC, only(basis_SC.kpoints))]
    dot(prefac_Hs .* ψ1, ψ2)
end

function Hs_norm(basis_SC::PlaneWaveBasis, ψ; s=0, tol=1e-10)
    square_norm = safereal(Hs_dot(basis_SC, ψ, ψ; s); tol)
    (square_norm < 0) && error("negative Hs_norm, try raising minimal spread ζ_min")
    sqrt(square_norm)
end

function Hs_overlap(basis_SC::PlaneWaveBasis, Xs_fourier; s=0)
    num_aos = length(Xs_fourier)
    S = zeros(eltype(Xs_fourier[1]), num_aos, num_aos)
    for μ in 1:num_aos
        for ν in μ:num_aos
            S[μ, ν] = Hs_dot(basis_SC, Xs_fourier[μ], Xs_fourier[ν]; s)
            (μ≠ν) && (S[ν, μ] = S[μ, ν])
        end
    end
    Hermitian(S)
end

@doc raw"""
TODO
"""
function project_wannier_on_basis(Wc::CompressedWannier)
    s = Wc.error_norm
    basis_SC = Wc.basis_supercell

    # Check that Ψ and all AOs have been converted to supercell conventions
    Φs_Four = [Φ(basis_SC) for Φ in Wc.basis_functions]
    normalize!.(Φs_Four) # Avoids weird coefficients

    # Compute the coefficients of the projection
    S = Hs_overlap(basis_SC, Φs_Four; s)
    Χ = [Hs_dot(basis_SC, Wc.wannier, Φ; s) for Φ in Φs_Four]
    # Check for conditioning issues before inverting and stop if conditioning is to high
    # TODO
    # (cond(S) > 1e8) && (error("cond(S)>1e8"))

    C_opti = filter_small_coeffs.(Hermitian(S)\Χ)
    # Assemble projection in Fourier
    sum( c .* Φ for (c, Φ) in zip(C_opti, Φs_Four) ), C_opti
end
