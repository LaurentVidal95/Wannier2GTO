# TODO: Complete documentation.

"""
Computes the Hˢ scalar product of two vectors containing all (k+G) Fourier coefficients.
In DFTK conventions, these vectors are Bloch wave of the supercell Gamma point.
Conversion from cell to supercell is hence needed beforehand.
"""
function Hˢ_dot(basis_supercell::PlaneWaveBasis, ψ1::AbstractVector{T1},
                ψ2::AbstractVector{T2}; s=0) where {T1, T2<:Complex}
    # Handle trivial case (L² norm)
    iszero(s) && return dot(ψ1, ψ2)
    # Hˢ norm (s ≠ 0)
    prefac_Hˢ = [(1+norm(Gpk)^2)^s
                for Gpk in G_vectors_cart(basis_supercell, only(basis_supercell.kpoints))]
    dot(prefac_Hˢ .* ψ1, ψ2)
end

function Hˢ_norm(basis_supercell::PlaneWaveBasis, ψ; s=0, tol=1e-10)
    square_norm = safereal(Hˢ_dot(basis_supercell, ψ, ψ; s); tol)
    (square_norm < 0) && error("negative Hˢ_norm, try raising minimal spread ζ_min")
    sqrt(square_norm)
end

function Hˢ_overlap(basis_supercell::PlaneWaveBasis, Xs_fourier; s=0)
    num_aos = length(Xs_fourier)
    S = zeros(eltype(Xs_fourier[1]), num_aos, num_aos)
    for μ in 1:num_aos
        for ν in μ:num_aos
            S[μ, ν] = Hˢ_dot(basis_supercell, Xs_fourier[μ], Xs_fourier[ν]; s)
            (μ≠ν) && (S[ν, μ] = S[μ, ν])
        end
    end
    Symmetric(S)
end

function Hˢ_overlap(Ms::T; s=0) where T # T is a vector of GaussianPolynomial or BasisFunctions
    S = [integral(Mi, Mj; type=:overlap) for Mi in Ms, Mj in Ms]
    (s==1) && (S .+= [integral(Mi, Mj; type=:laplacian) for Mi in Ms, Mj in Ms])
    S
end

function project_wannier_on_basis(Wc::CompressedWannier)
    s = Wc.error_norm
    isempty(Wc.basis_supercell) && error("The planewave basis is missing.")
    basis_supercell = Wc.basis_supercell

    # Compute the Fourier coefficients of all basis functions on the supercell
    # planewave basis.
    Φs = Wc.basis_functions
    Φs_Four = [Φ(basis_supercell) for Φ in Φs]

    # Compute the projection
    S = Hˢ_overlap(basis_supercell, Φs_Four; s) # Replace by fast analytic overlap if compatible with AD
    Χ = [Hˢ_dot(basis_supercell, Wc.wannier, Φ; s) for Φ in Φs_Four] # slow numeric dot prod

    # TODO: add cure for conditioning if needed.
    # Check for conditioning issues before inverting and stop if conditioning is to high
    # (cond(S) > 1e8) && (error("cond(S)>1e8"))
    C_opti = filter_small_coeffs.(S\Χ)

    # Assemble projection in Fourier
    sum( c .* Φ for (c, Φ) in zip(C_opti, Φs_Four) ), C_opti
end

function _check_L²_dot_precision(basis_supercell::PlaneWaveBasis; tol=1e-5)
    dummy_SAGTO = GaussianPolynomial([(0,0,1)], [1.], zeros(3), 1.)(basis_supercell)
    # Typical problematic num dot prod
    X1 = GaussianPolynomial([(1,0,1)], [1.], zeros(3), 3.)
    X2 = GaussianPolynomial([(1,0,1)], [1.], zeros(3), 1/2.)
    # Check difference between analytic and discrete dot product
    error = abs(integral(X1,X2) - dot(X1(basis_supercell), X2(basis_supercell)))
    if error > tol
        @warn "The numerical L² dot product is off by $error. You"*
            " want to raise the Ecut used for wannierization"
    end
    nothing
end
