import LinearAlgebra.normalize

struct BasisFunction{T<:Real}
    coeffs :: Vector{T}
    SAGTOs # Vector of Gaussian polynomials
end

function (Φ::BasisFunction)(basis_supercell::PlaneWaveBasis; normalize_SAGTO=true)
    SAGTOs_Fourier = [X(basis_supercell; normalize_SAGTO) for X in Φ.SAGTOs]
    sum(λ * SAGTOs_Fourier[i] for  (i, λ) in enumerate(Φ.coeffs))
end

function enforce_D3_symmetry(Φ::BasisFunction)
    SAGTOs = Φ.SAGTOs
    center = SAGTOs[1].center

    # Add SAGTOs with rotated centers at +2π/3 and +4π/3
    SAGTOs_ω = [GaussianPolynomial(X, rot(2π/3)*center) for X in SAGTOs]
    SAGTOs_ω_square = [GaussianPolynomial(X, rot(4π/3)*center) for X in SAGTOs]
    # Return new basis function
    new_coeffs = vcat([Φ.coeffs for _ in 1:3]...)
    BasisFunction(new_coeffs, vcat(SAGTOs, SAGTOs_ω, SAGTOs_ω_square))
end

"""
Construct a linear combination of given SAGTOs that best
approximates the residual contained in the CompressedWannier struct.
"""
function optimal_basis_function(Wc, SAGTOs; tol=1e-5)
    # Extract needed data
    s = Wc.error_norm
    basis_supercell = Wc.basis_supercell
    center = SAGTOs[1].center

    # Compute optimal coefficients for given SAGTOs
    SAGTOs_Four = [X(basis_supercell) for X in SAGTOs]
    Γ = ThreadsX.map(X->Hˢ_dot(basis_supercell, Wc.residual, X; s), SAGTOs_Four)
    S = Hˢ_overlap(basis_supercell, SAGTOs_Four; s)
    optimal_coeffs = filter_small_coeffs.(S\Γ; tol)
    
    # Assemble optimal basis function and enforce D3 symmetry if needed
    Φ = BasisFunction(optimal_coeffs, SAGTOs)
    !iszero(center - Wc.center) &&  (Φ = enforce_D3_symmetry(Φ))
    Φ, Hˢ_norm(basis_supercell, Φ(basis_supercell) - Wc.residual; s=Wc.error_norm)^2
end

function normalize(Φ::BasisFunction)
    norm_Φ = √(integral(Φ, Φ; type=:overlap))
    normalized_coeffs = Φ.coeffs ./ norm_Φ
    BasisFunction(normalized_coeffs, Φ.SAGTOs)
end
