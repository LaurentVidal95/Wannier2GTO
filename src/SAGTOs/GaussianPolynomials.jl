import Base.+

"""
Contains all parameters to describe a gaussian-polynomial g such that:
``g(x) = ( ∑_{(n_x,n_z,n_z)} λ_{(n_x,n_y,n_z)}x^{n_x}y^{n_y}z^{n_z} ) 
                                * exp(-ζ*norm([x,y,z] - α)^2)``
Parameters are:
   • ``pol.exps = [(n_x,n_y,n_z) ... ]`` contains all 3 polynomials powers ``(n_x,n_y,n_z)``
   • ``pol.coeffs = [λ_{(n_x,n_y,n_z)}, ∀(n_x,n_y,n_z) ∈ pol_exps]``
   • ``ζ`` is the spread of the gaussian part
   • ``center = [α_x, α_y, α_z]`` in cartesian coordinates

 Types are specific for each parameter to be compatible with Forward Diff

"""
struct GaussianPolynomial{T1<:Real, T2<:Real, F1, F2<:Function}
    pol::F1
    center::Vector{T1}
    spread::T2
end

function GaussianPolynomial(exps::Vector{Tuple{Int64, Int64, Int64}},
                            coeffs::Vector{T1}, center::Vector{T2},
                            spread::T3) where {T1, T2, T3 <: Real}
    # Construct polynomial part
    @polyvar x y z
    static_pol = Polynomial(sum( prod([x,y,z] .^ exps)*λ for (exps, λ) in zip(exps, coeffs) ))
    pol(X::AbstractVector) = evaluate(static_pol, X)
    GaussianPolynomial(pol, center, spread, X_hat) 
end
(X::GaussianPolynomial)(R)= X.pol(R .- X.center) * exp(-X.spread*norm(R .- X.center)^2)
(Χ::GaussianPolynomial)(basis_SC::PlaneWaveBasis) = fft_supercell(basis_SC, Χ)


"""
Compute the Bloch decomposition (stored as Fourier coefficients as in scfres.ψ) of a given
GaussianPolynomial using DFTK FFT routines. For now very slow somehow...
Maybe because of the convert part.
"""
function fft_supercell(basis_SC::PlaneWaveBasis, X::GaussianPolynomial)
    # Shift X to the center of the cell to avoid sampling issues
    T = eltype(basis_SC)
    shift = sum(eachcol((basis_SC.model.lattice .- X.center) ./ 2))

    # Compute FFT of X
    X_real = map(r_cart -> X(r_cart - shift), r_vectors_cart(basis_SC))
    X_fourier = fft(basis_SC, X_real)[basis_SC.kpoints[1].mapping]

    # Shift back to original center
    X_fourier .*= ThreadsX.map(Gpk_cart-> cis(dot(Gpk_cart, shift)),
                      G_vectors_cart(basis_SC, only(basis_SC.kpoints)))
    X_fourier ./ norm(X_fourier)
end
