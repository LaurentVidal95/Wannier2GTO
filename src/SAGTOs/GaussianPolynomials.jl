using ForwardDiff
using ThreadsX

import Base.+
    
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
end
(g::GaussianPolynomial)(X)= g.pol(X .- g.α) * exp(-g.ζ*norm(X .- g.α)^2)
GaussianPolynomial(exps, coeffs, α, ζ) = GaussianPolynomial(PolynomialPart(exps, coeffs),
                                                            α, ζ)
"""
    Define sum for gaussians with same centers and spreads
    PB: Takes more time than summing directly tables of Fourier coefficients
"""
function (+)(g1::GaussianPolynomial{T}, g2::GaussianPolynomial{T}) where {T<:Real}
    @assert ((g1.α == g2.α) &&
             (g1.ζ == g2.ζ)) "Only sum over gaussians with same centers and spreads."
    GaussianPolynomial(g1.pol + g2.pol, g1.α, g1.ζ)
end

(g::GaussianPolynomial)(basis::PlaneWaveBasis)=discrete_Bloch_transform(basis, g)

"""
Compute the Bloch decomposition (stored as Fourier coefficients as in scfres.ψ) of a given
GaussianPolynomial using DFTK FFT routines.
TODO: use supercell instead of cell if cell is to small to avoid sampling issues.
"""
function discrete_Bloch_transform(basis::PlaneWaveBasis, g::GaussianPolynomial)
    # Shift g to the center of the cell to avoid sampling issues
    shift = sum(ei ./ 2 for ei in eachcol(basis.model.lattice)) .- g.α
    g_fourier = r_to_G(basis, g.(collect(r_vectors_cart(basis)) .- Ref(shift)))
    # Re-arrange as g[kpt][G] table
    ThreadsX.map(basis.kpoints) do kpt
        # prefactor to shift it back directly in Fourier.
        shift_prefac = cis.(dot.(Gplusk_vectors_cart(basis, kpt), .-Ref(shift)))
        gk = view(g_fourier, kpt.mapping) .* shift_prefac
        gk ./ norm(gk)
    end
end

"""
Polynom described as a linear combination of monoms
"""
mutable struct PolynomialPart
    exps::AbstractVector # Powers of the involved monoms
    coeffs::AbstractVector # Corresponding linear combinantion coefficients.
end
(pol::PolynomialPart)(X) = sum( prod(X .^ exps)*λ
                               for (exps, λ) in zip(pol.exps, pol.coeffs) )
(+)(pol1::PolynomialPart, pol2::PolynomialPart) =
    PolynomialPart(vcat(pol1.exps, pol2.exps), vcat(pol1.coeffs, pol2.coeffs))
