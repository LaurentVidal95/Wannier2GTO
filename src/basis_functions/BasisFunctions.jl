struct BasisFunction2{TR<:Real}
    coeffs :: Vector{TR}
    SAGTOs # Vector of Gaussian polynomials
end
BasisFunction = BasisFunction2

function (Φ::BasisFunction)(basis_SC::PlaneWaveBasis)
    SAGTOs_Fourier = [X(basis_SC) for X in Φ.SAGTOs]
    sum(λ * SAGTOs_Fourier[i] for  (i, λ) in enumerate(Φ.coeffs))
end

function enforce_D3_symmetry(Φ::BasisFunction)
    SAGTOs = Φ.SAGTOs
    center = SAGTOs[1].center
    # Add SAGTOs with rotated centers at +2π/3 and +4π/3
    SAGTOs_ω = [GaussianPolynomial(X, rot(2π/3)*center) for X in SAGTOs]
    SAGTOs_ω_square = [GaussianPolynomial(X, rot(4π/3)*center) for X in SAGTOs]
    # Return new basis function
    BasisFunction(Φ.coeffs, vcat(SAGTOs, SAGTOs_ω, SAGTOs_ω_square))
end

"""
Construct a linear combination of given SAGTOs that best
approximates the residual contained in the CompressedWannier struct.
"""
function optimal_basis_function(Wc::CompressedWannier, SAGTOs;
                                s=0)
    basis_SC = Wc.basis_supercell
    center = SAGTOs[1].center

    # Compute optimal coefficients for given SAGTOs
    SAGTOs_Four = [X(basis_SC) for X in SAGTOs]
    Γ = ThreadsX.map(X->Hs_dot(basis_SC, Wc.residual, X; s), SAGTOs_Four)
    S = Hs_overlap(basis_SC, SAGTOs_Four; s)
    optimal_coeffs = S\Γ
    
    # Assemble optimal basis function and enforce D3 symmetry if needed
    Φ = BasisFunction(optimal_coeffs, SAGTOs)
    !iszero(center) && (Φ = enforce_D3_symmetry(Φ))
    Φ, Hs_norm(Wc.basis_SC, Φ - Wc.residual; s=Wc.error_norm)^2
end
