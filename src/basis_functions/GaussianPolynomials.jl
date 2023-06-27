"""
Contains all parameters to describe a gaussian-polynomial g such that:
``g(x) = ( ∑_{(n_x,n_z,n_z)} λ_{(n_x,n_y,n_z)}x^{n_x}y^{n_y}z^{n_z} ) 
                                * exp(-ζ*norm([x,y,z] - α)^2)``
Parameters are:
   • ``pol``, a StaticPolynomial.jl object that allow fast multivariate
   polynomial evalutation.
   • ``ζ`` is the spread of the gaussian part
   • ``center = [α_x, α_y, α_z]`` the center of the GTO in cartesian coordinates

 Types are specific for each parameter to be compatible with Forward Diff
"""
struct GaussianPolynomial{T1<:Real, T2<:Real}
    pol::Polynomial
    center::Vector{T1}
    spread::T2
end
function GaussianPolynomial(exps::Vector{Tuple{Int64, Int64, Int64}},
                            coeffs::Vector{T1}, center::Vector{T2},
                            spread::T3) where {T1, T2, T3 <: Real}
    # Construct polynomial part
    @polyvar x y z
    pol = Polynomial(sum( prod([x,y,z] .^ exps)*λ for (exps, λ) in zip(exps, coeffs) ))
    GaussianPolynomial(pol, center, spread)
end
function (X::GaussianPolynomial)(A::AbstractArray)
    @assert length(size(A)) == 3
    pol_part = evaluate.(Ref(X.pol), A)
    exp_part = ThreadsX.map(R->exp(-X.spread*norm(R .- X.center)^2), A)
    pol_part .* exp_part
end
(Χ::GaussianPolynomial)(basis_supercell::PlaneWaveBasis) = fft_supercell(basis_supercell, Χ)

GaussianPolynomial(X::GaussianPolynomial, center) = GaussianPolynomial(X.pol, center, X.spread)

"""
Compute the Bloch decomposition (stored as Fourier coefficients as in scfres.ψ) of a given
GaussianPolynomial using DFTK FFT routines.
"""
function fft_supercell(basis_supercell::PlaneWaveBasis, X::GaussianPolynomial)
    # Shift X to the center of the cell to avoid sampling issues
    T = eltype(basis_supercell)
    shift = sum(eachcol((basis_supercell.model.lattice .- X.center) ./ 2))

    # Compute FFT of X
    X_real = X(r_vectors_cart(basis_supercell) .- Ref(shift))
    X_fourier = fft(basis_supercell, X_real)[basis_supercell.kpoints[1].mapping]

    # Shift back to original center
    X_fourier .*= ThreadsX.map( Gpk_cart-> cis(dot(Gpk_cart, shift)),
                      G_vectors_cart(basis_supercell, only(basis_supercell.kpoints)) )
    # normalize!(X_fourier)
end
