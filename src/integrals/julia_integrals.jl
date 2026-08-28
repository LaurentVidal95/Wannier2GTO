#
# Pure-Julia analytic integrals over Gaussian-polynomial primitives.
#
# CONVENTION (matches the Fourier transform, `translate`, `rotate`,
# `enforce_D3_symmetry`, `SAGTO_basis`, and standard quantum chemistry):
# the polynomial is expressed in coordinates **relative to the center**,
#
#     X(r) = (x-α_x)^{n_x} (y-α_y)^{n_y} (z-α_z)^{n_z} · exp(-ζ‖r-α‖²).
#
# This matters: an absolute-coordinate reading agrees with the relative one only
# when α = 0, and diverges badly otherwise (⟨s@0, px@[1,0,0]⟩ differs in sign).
# Only the relative convention makes `translate`/`rotate` actual translations
# and rotations of the function, and makes symmetry-adapted polynomials
# symmetric about their own center.
#

@doc raw"""
1D Gaussian moment, centered:

```math
M_n(\zeta) = \int_{-\infty}^{+\infty} u^n\, e^{-\zeta u^2}\,du
= \begin{cases}
    0                                                    & n\text{ odd} \\
    \dfrac{(n-1)!!}{(2\zeta)^{n/2}}\,\sqrt{\dfrac{\pi}{\zeta}} & n\text{ even}
  \end{cases}
```
"""
@inline function gaussian_moment_1d(n::Int, ζ)
    isodd(n) && return zero(float(ζ))
    n == 0 && return sqrt(π / ζ)
    df = float(prod(1:2:(n - 1)))  # (n-1)!! for even n ≥ 2
    df / (2ζ)^(n ÷ 2) * sqrt(π / ζ)
end

@doc raw"""
1D two-center moment appearing in the overlap of two Gaussian-polynomial
primitives, with **both** monomials relative to their own center. Expanding
each about the Gaussian-product center ``P`` (binomial) reduces it to centered
moments:

```math
\int_{\mathbb{R}} (u-A)^i (u-B)^j\, e^{-\zeta_p (u-P)^2}\,du
= \sum_{k=0}^{i}\sum_{l=0}^{j} \binom{i}{k}\binom{j}{l}
  (P-A)^{i-k}\,(P-B)^{j-l}\, M_{k+l}(\zeta_p).
```
"""
@inline function two_center_moment_1d(i::Int, j::Int, PA, PB, ζp)
    T = float(promote_type(typeof(PA), typeof(PB), typeof(ζp)))
    s = zero(T)
    for k in 0:i, l in 0:j
        s += binomial(i, k) * binomial(j, l) *
             PA^(i - k) * PB^(j - l) * gaussian_moment_1d(k + l, ζp)
    end
    s
end

@doc raw"""
Pure-Julia analytic overlap of two Gaussian-polynomial primitives, monomials
taken relative to their respective centers. Drop-in replacement for
`GaIn.overlap` for the `:overlap` integral type.

```math
\langle X_1, X_2 \rangle = \int_{\mathbb{R}^3}
    \prod_j (r_j - A_j)^{n_{j1}} (r_j - B_j)^{n_{j2}}\,
    e^{-\zeta_1\|\mathbf{r}-\mathbf{A}\|^2}\,
    e^{-\zeta_2\|\mathbf{r}-\mathbf{B}\|^2}\, d^3\mathbf{r}.
```

Strategy: the Gaussian product theorem reduces the two-Gaussian factor to a
single Gaussian of spread ``\zeta_p = \zeta_1+\zeta_2`` centered at
``\mathbf{P} = (\zeta_1\mathbf{A} + \zeta_2\mathbf{B})/\zeta_p``, times
``K = \exp(-\zeta_1\zeta_2/\zeta_p\,\|\mathbf{A}-\mathbf{B}\|^2)``. What remains
separates over (x,y,z) into `two_center_moment_1d`.

ForwardDiff-compatible (only arithmetic + sqrt/exp on `ζ` and the centers).
"""
function overlap_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                       ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    ζp = ζ1 + ζ2
    Rp = (ζ1 .* R1 .+ ζ2 .* R2) ./ ζp
    K = exp(-ζ1 * ζ2 / ζp * sum(abs2, R1 .- R2))
    n1 = (nx1, ny1, nz1)
    n2 = (nx2, ny2, nz2)
    K * prod(two_center_moment_1d(n1[j], n2[j], Rp[j] - R1[j], Rp[j] - R2[j], ζp)
             for j in 1:3)
end

@doc raw"""
Per-axis monomial terms of ``\partial^2/\partial u^2`` acting on
``v^a\, e^{-\zeta v^2}`` with ``v = u - R`` the coordinate **relative to the
center**. In this convention the center drops out entirely and only three terms
survive:

```math
\partial_v^2\left[v^a e^{-\zeta v^2}\right]
  = \Big[a(a-1)\,v^{a-2} - 2\zeta(2a+1)\,v^{a} + 4\zeta^2\,v^{a+2}\Big] e^{-\zeta v^2}.
```

Returns `(coeff, new_exponent)` tuples, the exponents again relative.
"""
function _axis_laplacian_terms(a::Int, ζ::T) where {T}
    C = float(T)
    terms = Tuple{C, Int}[]
    a ≥ 2 && push!(terms, (C(a * (a - 1)), a - 2))
    push!(terms,           (-2ζ * (2a + 1), a    ))
    push!(terms,           (C(4ζ^2),        a + 2))
    terms
end

@doc raw"""
Pure-Julia positive Laplacian form of two Gaussian-polynomial primitives.
Drop-in replacement for `GaIn.laplacian`, matching the sign convention required
by `Hˢ_overlap(Ms; s=1)` (which adds this to the L² overlap to build the H¹
overlap):

```math
\texttt{laplacian\_julia}(X_1, X_2) = \langle \nabla X_1, \nabla X_2\rangle
    = -\langle X_1, \nabla^2 X_2\rangle .
```

``\nabla^2 X_2`` is a finite sum of monomials sharing the spread and center of
``X_2``, so each contribution is an exact `overlap_julia` call. No Boys
function, no McMurchie-Davidson. ForwardDiff-compatible.
"""
function laplacian_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                         ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    n2 = (nx2, ny2, nz2)
    acc = zero(float(promote_type(typeof(ζ1), typeof(ζ2),
                                  eltype(R1), eltype(R2))))
    for axis in 1:3                      # ∇²X₂ = Σ_axis ∂²_axis X₂
        for (c, e) in _axis_laplacian_terms(n2[axis], ζ2)
            m = ntuple(i -> i == axis ? e : n2[i], 3)
            acc += c * overlap_julia(ζ1, R1, nx1, ny1, nz1,
                                     ζ2, R2, m[1], m[2], m[3])
        end
    end
    -acc    # ⟨∇X₁,∇X₂⟩ = -⟨X₁,∇²X₂⟩
end

@doc raw"""
Kinetic-energy form of two Gaussian-polynomial primitives,

```math
\langle X_1, -\tfrac12\nabla^2 X_2\rangle = \tfrac12\,\texttt{laplacian\_julia}(X_1, X_2).
```

Drop-in replacement for `GaIn.kinetic`.
"""
function kinetic_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                       ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    laplacian_julia(ζ1, R1, nx1, ny1, nz1, ζ2, R2, nx2, ny2, nz2) / 2
end
