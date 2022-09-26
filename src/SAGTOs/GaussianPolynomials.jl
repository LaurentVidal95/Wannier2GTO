using ForwardDiff
using ThreadsX

import Base.+
    
"""
Polynom described as a linear combination of monoms
"""
mutable struct PolynomialPart
    exps::AbstractVector   # Powers of the involved monoms
    coeffs::AbstractVector # Corresponding linear combinantion coefficients.
end
(pol::PolynomialPart)(X) = sum( prod(X .^ exps)*λ
                               for (exps, λ) in zip(pol.exps, pol.coeffs) )
(+)(pol1::PolynomialPart, pol2::PolynomialPart) =
    PolynomialPart(vcat(pol1.exps, pol2.exps), vcat(pol1.coeffs, pol2.coeffs))

"""
Contains all parameters to describe a gaussian-polynomial g such that:
``g(x) = ( ∑_{(n_x,n_z,n_z)} λ_{(n_x,n_y,n_z)}x^{n_x}y^{n_y}z^{n_z} ) 
                            * exp(-ζ*norm([x,y,z] - α)^2)``
Parameters are:
   • ``pol.exps = [(n_x,n_y,n_z) ... ]`` contains all 3 polynomials powers ``(n_x,n_y,n_z)``
   • ``pol.coeffs = [λ_{(n_x,n_y,n_z)}, ∀(n_x,n_y,n_z) ∈ pol_exps]``
   • ``ζ`` is the spread of the gaussian part
   • ``α = [α_x, α_y, α_z]`` in cartesian coordinates
"""

mutable struct GaussianPolynomial{T<:Real}
    pol::PolynomialPart # See below
    α::Vector{T}
    ζ::T
    Fourier_transform
end
(Χ::GaussianPolynomial)(R)= Χ.pol(R .- Χ.α) * exp(-Χ.ζ*norm(R .- Χ.α)^2)
F(Χ::GaussianPolynomial) = Χ.Fourier_transform

function ∂n(h, n, x)
    (n==0) && (return h(x))
    (n==1) && (return ForwardDiff.derivative(h, x))
    ∂n(y->ForwardDiff.derivative(h, y), n-1, x)
end

function GaussianPolynomial(exps, coeffs, α, ζ)
    g_hat(q) = exp(-(q^2)/(4ζ))
    Χ_hat(q) = cis(-dot(q,α)) *
        ThreadsX.sum( λμ * prod( (im^nj)*∂n(y->g_hat(y), nj, qj)
                                 for (nj, qj) in zip(exp_μ, q) )
                      for (exp_μ, λμ) in zip(exps, coeffs)
                    )
    GaussianPolynomial(PolynomialPart(exps, coeffs), α, ζ, Χ_hat)
end
GaussianPolynomial(pol::PolynomialPart, α, ζ) = GaussianPolynomial(pol.exps, pol.coeffs, α, ζ)

"""
    Define sum for gaussians with same centers and spreads
    PB: Takes more time than summing directly tables of Fourier coefficients
"""
function (+)(Χ1::GaussianPolynomial{T}, Χ2::GaussianPolynomial{T}) where {T<:Real}
    @assert ((Χ1.α == Χ2.α) &&
             (Χ1.ζ == Χ2.ζ)) "Only sum over gaussians with same centers and spreads."
    GaussianPolynomial(Χ1.pol + Χ2.pol, Χ1.α, Χ1.ζ, F(Χ1) + F(Χ2))
end

(Χ::GaussianPolynomial)(basis_SC::PlaneWaveBasis) = discrete_Bloch_transform(basis_SC, Χ)

"""
Evaluate Χ at all k+G vectors given a basis and returns in supercell convention
(i.e. single vector containing all k+G coefficient in order given by
G_vectors(basis_supercell, basis_supercell.Γ_point))
"""
function discrete_Bloch_transform(basis_SC::PlaneWaveBasis, Χ::GaussianPolynomial)
    Χ_fourier = F(Χ).(G_vectors_cart(basis_SC, only(basis_SC.kpoints)))
    Χ_fourier ./ norm(Χ_fourier)
end

# """
# Compute the Bloch decomposition (stored as Fourier coefficients as in scfres.ψ) of a given
# GaussianPolynomial using DFTK FFT routines.
# """
# function discrete_Bloch_transform(basis_SC::PlaneWaveBasis, Χ::GaussianPolynomial)
#     # Shift Χ to the center of the cell to avoid sampling issues
#     shift = sum(ei ./ 2 for ei in eachcol(basis_SC.model.lattice)) .- Χ.α
#     Χ_fourier = r_to_G(basis_SC, Χ.(collect(r_vectors_cart(basis_SC)) .- Ref(shift)))
#     # Shift back in fourier for plotting
#     for (iG, G) in enumerate(G_vectors_cart(basis_SC))
#         Χ_fourier[iG] *= cis(dot(G, shift))
#     end
#     # Return only the coefficient corresponding to non zero wannier coefficients.
#     Χ_vec = Χ_fourier[basis_SC.kpoints[1].mapping]
#     Χ_vec ./ norm(Χ_vec)
# end
