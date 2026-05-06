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

"""
Physicists' Hermite polynomial H_n(t), defined by the recurrence
  H_0 = 1,  H_1 = 2t,  H_n = 2t·H_{n-1} − 2(n−1)·H_{n-2}.
Returns an expression that is differentiable by Zygote w.r.t. t.
"""
function _hermite_phys(n::Int, t)
    (n == 0) && return one(t)
    (n == 1) && return 2 * t
    h_prev2 = one(t)
    h_prev1 = 2 * t
    for k in 2:n
        h_curr = 2 * t * h_prev1 - 2 * (k - 1) * h_prev2
        h_prev2 = h_prev1
        h_prev1 = h_curr
    end
    h_prev1
end

"""
Analytic n-th derivative of g(q) = exp(-q²/(4ζ)) w.r.t. q:
  g^(n)(q) = (-1/(2√ζ))^n · H_n(q/(2√ζ)) · exp(-q²/(4ζ))
where H_n is the physicists' Hermite polynomial.
Fully differentiable by Zygote w.r.t. both q and ζ.
"""
@inline function _dghat_dn(n::Int, q, spread)
    inv2sqrtζ = inv(2 * sqrt(spread))
    t = q * inv2sqrtζ
    (-inv2sqrtζ)^n * _hermite_phys(n, t) * exp(-q^2 / (4 * spread))
end

function GaussianPolynomial(exps::Vector{Tuple{Int64, Int64, Int64}},
                            coeffs::AbstractVector{T1}, center::AbstractVector{T2},
                            spread::T3;
                            normalize_SAGTO=true) where {T1, T2, T3 <: Real}
    # ensure that the given Gaussian Polynomial is normalized
    prefac = normalize_SAGTO ? analytic_norm(exps, filter_dual.([coeffs, center, spread])...) : 1.
    coeffs = coeffs ./ prefac

    # Construct polynomial part (symbolic — invisible to AD via Zygote.ignore)
    pol = Zygote.ignore() do
        @polyvar x y z
        Polynomial( sum(prod([x,y,z] .^ exp_μ)*λ for (exp_μ, λ) in zip(exps, coeffs)) )
    end

    # Fourier part: AD-friendly closure (uses _dghat_dn, no ForwardDiff)
    X_hat = SAGTO_fourier_transform(exps, coeffs, center, spread)
    GaussianPolynomial(pol, center, spread, X_hat)
end
function SAGTO_fourier_transform(exps, coeffs, center, spread)
    # Freeze the integer exponent data outside the AD tape.
    # `exps` are tuples of Ints — never differentiated.
    exps_frozen = Zygote.ignore(() -> collect(exps))  # Vector{Tuple{Int,Int,Int}}
    n_terms = length(exps_frozen)
    prefac = (π / spread)^(3/2)
    # AD-friendly closure: differentiates w.r.t. coeffs, center, spread.
    X_hat(q) = begin
        phase = cis(-dot(q, center))
        term_sum = sum(1:n_terms) do μ
            exp_μ = Zygote.ignore(() -> exps_frozen[μ])  # Tuple{Int,Int,Int}: constant
            gd = _dghat_dn(exp_μ[1], q[1], spread) *
                 _dghat_dn(exp_μ[2], q[2], spread) *
                 _dghat_dn(exp_μ[3], q[3], spread)
            im_factor = Zygote.ignore(() -> im^exp_μ[1] * im^exp_μ[2] * im^exp_μ[3])
            coeffs[μ] * im_factor * gd
        end
        prefac * phase * term_sum
    end
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
    # exps are integer tuples — freeze them outside the AD tape.
    exps, coeffs = Zygote.ignore(() -> pol_to_arrays(X.pol))
    GaussianPolynomial(exps, coeffs, center, X.spread)
end

function (Base.:*)(X1::GaussianPolynomial, X2::GaussianPolynomial;
                   normalize_SAGTO=false)
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
    # Non-mutating form: AD-friendly for Zygote (no .= or normalize!)
    Gs = G_vectors_cart(basis_supercell, only(basis_supercell.kpoints))
    vol_factor = √(basis_supercell.model.unit_cell_volume)
    X_fourier = ℱ.(Ref(X), Gs) ./ vol_factor
    normalize_SAGTO ? X_fourier ./ norm(X_fourier) : X_fourier
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
