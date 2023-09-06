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
struct GaussianPolynomial{T1<:Real, T2<:Real, F<:Function}
    pol::Polynomial
    center::AbstractVector{T1}
    spread::T2
    # Hack to avoid boundary problems
    Fourier_transform::F
end

function pol_to_arrays(pol::Polynomial)
    exps = Tuple.(eachcol(StaticPolynomials.exponents(pol)))
    coeffs = pol.coefficients
    exps, coeffs
end

function ∂n(h::F, n::Int64, x::T) where {T<:Real, F<:Function}
    (n==0) && (return h(x))
    (n==1) && (return ForwardDiff.derivative(h, x))
    ∂n(y->ForwardDiff.derivative(h, y), n-1, x)
end
function GaussianPolynomial(exps::Vector{Tuple{Int64, Int64, Int64}},
                            coeffs::AbstractVector{T1}, center::AbstractVector{T2},
                            spread::T3;
                            normalize_SAGTO=true) where {T1, T2, T3 <: Real}
    # ensure that the given Gaussian Polynomial is normalized
    prefac = normalize_SAGTO ? analytic_norm(exps, filter_dual.([coeffs, center, spread])...) : 1.
    coeffs = coeffs ./ prefac

    # Construct polynomial part
    @polyvar x y z
    pol = Polynomial( sum(prod([x,y,z] .^ exps)*λ for (exps, λ) in zip(exps, coeffs)) )

    # Fourier part
    X_hat = SAGTO_fourier_transform(exps, coeffs, center, spread)
    GaussianPolynomial(pol, center, spread, X_hat)
end
function SAGTO_fourier_transform(exps, coeffs, center, spread)
    g_hat(q) = exp(-(q^2)/(4*spread))
    prefac = (π/spread)^(3/2)
    X_hat(q) = prefac * cis(-dot(q,center)) *
        sum( λμ * prod( (im^nj)*∂n(y->g_hat(y), nj, qj)
                        for (nj, qj) in zip(exp_μ, q) )
             for (exp_μ, λμ) in zip(exps, coeffs)
             )
    X_hat
end

function (X::GaussianPolynomial)(A::AbstractArray)
    pol_part = evaluate.(Ref(X.pol), A)
    exp_part = ThreadsX.map(R->exp(-X.spread*norm(R .- X.center)^2), A)
    pol_part .* exp_part
end
(X::GaussianPolynomial)(basis_supercell::PlaneWaveBasis; normalize_SAGTO=true) =
    slow_fourier_transform_supercell(basis_supercell, X; normalize_SAGTO)
ℱ(X::GaussianPolynomial, x)= X.Fourier_transform(x)

function GaussianPolynomial(X::GaussianPolynomial, center)
    exps, coeffs = pol_to_arrays(X.pol)
    GaussianPolynomial(exps, coeffs, center, X.spread)
end

function (Base.:*)(X1::GaussianPolynomial, X2::GaussianPolynomial;
                   normalize_SAGTO=true)
    # Compute new spread and center
    ζ₁ = X1.spread; R₁ = X1.center
    ζ₂ = X2.spread; R₂ = X2.center
    ζ₃ = ζ₁ + ζ₂
    R₃ = (ζ₁*R₁ + ζ₂*R₂)/ζ₃

    # Compute new polynomial
    α = (ζ₁*norm(R₁)^2 + ζ₂*norm(R₂)^2) - norm(ζ₁*R₁+ζ₂*R₂)^2/ζ₃
    prefac = exp(-α)

    exp1, coeffs1 = pol_to_arrays(X1.pol)
    exp2, coeffs2 = pol_to_arrays(X2.pol)
    exp3 = eltype(exp1)[]
    coeffs3 = eltype(coeffs1)[]
    for ((nx1,ny1,nz1), λ1) in zip(exp1, coeffs1)
        for ((nx2,ny2,nz2), λ2) in zip(exp2, coeffs2)
            push!(exp3, (nx1+nx2, ny1+ny2, nz1+nz2))
            push!(coeffs3, λ1*λ2*prefac)
        end
    end

    GaussianPolynomial(exp3, coeffs3, R₃, ζ₃; normalize_SAGTO)
end

"""
Same as above but clumsy, hack with AD to avoid boundary issues.
This routines gives the proper coefficients up to a prefactor. The prefactor doesn't
matter since the basis functions are always normalized. The "normalize_SAGTO"
keyword is used for debugging. Beware: a small Ecut results in numerical errors.
(Ecut=15 ≈ 2e-4 precision).
"""
function slow_fourier_transform_supercell(basis_supercell::PlaneWaveBasis,
                                          X::GaussianPolynomial;
                                          normalize_SAGTO=true)
    X_fourier = ℱ.(Ref(X), G_vectors_cart(basis_supercell, only(basis_supercell.kpoints)))
    X_fourier ./= √(basis_supercell.model.unit_cell_volume)
    # Integrate a normalization step to the fft to avoid a lot of recomputations
    # Indeed, it is better to have X normalized during compression to avoid numerical
    # problems with very small norm and very big coefficients.
    # At the end of the procedure, all coefficients in X are recomputed so that it is
    # really normalized. See "fix_coefficients"
    (normalize_SAGTO) && normalize!(X_fourier)
    X_fourier
end

# Fast implementation using fft that suffers from issues when the spreads is to high.
#
# """
# Compute the Bloch decomposition (stored as Fourier coefficients as in scfres.ψ) of a given
# GaussianPolynomial using DFTK FFT routines.
# Suffers from boundary issues when the spread is to high (≥ 1 or so for Ecut=15)
# """
function fft_supercell(basis_supercell::PlaneWaveBasis, X::GaussianPolynomial;
                       normalize_SAGTO=true)
    # Shift X to the center of the cell to avoid sampling issues
    T = eltype(basis_supercell)
    shift = sum(eachcol((basis_supercell.model.lattice .- X.center) ./ 2))

    # Compute FFT of X
    X_real = X(r_vectors_cart(basis_supercell) .- Ref(shift))
    X_fourier = fft(basis_supercell, X_real)[basis_supercell.kpoints[1].mapping]

    # Shift back to original center
    X_fourier .*= ThreadsX.map( Gpk_cart-> cis(dot(Gpk_cart, shift)),
                      G_vectors_cart(basis_supercell, only(basis_supercell.kpoints)) )
    (normalize_SAGTO) && normalize!(X_fourier)
    X_fourier
end
