@doc raw"""
1D Gaussian moment

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
1D Gaussian moment around a non-zero center, by binomial expansion of ``(u-R+R)^n``:

```math
\int_{-\infty}^{+\infty} u^n\, e^{-\zeta(u-R)^2}\,du
= \sum_{k=0}^n \binom{n}{k}\, R^{\,n-k}\, M_k(\zeta).
```
"""
@inline function gaussian_moment_1d_centered(n::Int, R, ζ)
    T = float(promote_type(typeof(R), typeof(ζ)))
    s = zero(T)
    for k in 0:n
        s += binomial(n, k) * R^(n - k) * gaussian_moment_1d(k, ζ)
    end
    s
end

@doc raw"""
Pure-Julia analytic overlap of two unnormalized Gaussian-polynomial primitives.
Drop-in replacement for `GaIn.overlap` for the `:overlap` integral type.

```math
\langle X_1, X_2 \rangle = \int_{\mathbb{R}^3}
    x^{n_{x_1}} y^{n_{y_1}} z^{n_{z_1}}\, e^{-\zeta_1\|\mathbf{r}-\mathbf{R}_1\|^2}
    \cdot
    x^{n_{x_2}} y^{n_{y_2}} z^{n_{z_2}}\, e^{-\zeta_2\|\mathbf{r}-\mathbf{R}_2\|^2}\, d^3\mathbf{r}.
```

Strategy: Gaussian product theorem reduces the two-Gaussian integrand to a
single Gaussian of spread ``\zeta_p = \zeta_1+\zeta_2`` centered at
``\mathbf{R}_p = (\zeta_1 \mathbf{R}_1 + \zeta_2 \mathbf{R}_2)/\zeta_p``,
times a constant ``K = \exp(-\zeta_1\zeta_2/\zeta_p\,\|\mathbf{R}_1-\mathbf{R}_2\|^2)``.
The remaining integral separates over (x,y,z) into 1D centered Gaussian moments.

ForwardDiff-compatible (only arithmetic + sqrt/exp on `ζ` and `R`).
"""
function overlap_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                       ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    ζp = ζ1 + ζ2
    Rp = (ζ1 .* R1 .+ ζ2 .* R2) ./ ζp
    K = exp(-ζ1 * ζ2 / ζp * sum(abs2, R1 .- R2))
    K *
    gaussian_moment_1d_centered(nx1 + nx2, Rp[1], ζp) *
    gaussian_moment_1d_centered(ny1 + ny2, Rp[2], ζp) *
    gaussian_moment_1d_centered(nz1 + nz2, Rp[3], ζp)
end

@doc raw"""
Per-axis monomial terms of ``\partial^2/\partial u^2`` acting on
``u^a\, e^{-\zeta (u-R)^2}``, written in the absolute coordinate ``u`` (the
Gaussian factor is left unchanged). Returns `(coeff, new_exponent)` tuples:

```math
\partial_u^2 = a(a-1)\,u^{a-2}
             + \big(-2\zeta(1+2a) + 4\zeta^2 R^2\big)\,u^{a}
             + 4\zeta a R\, u^{a-1}
             - 8\zeta^2 R\, u^{a+1}
             + 4\zeta^2\, u^{a+2}.
```
"""
function _axis_laplacian_terms(a::Int, ζ::T, R) where {T}
    C = float(promote_type(T, typeof(R)))
    terms = Tuple{C, Int}[]
    a ≥ 2 && push!(terms, (C(a * (a - 1)),                  a - 2))
    push!(terms,          (-2ζ * (1 + 2a) + 4ζ^2 * R^2,     a    ))
    a ≥ 1 && push!(terms, (4ζ * a * R,                      a - 1))
    push!(terms,          (-8ζ^2 * R,                       a + 1))
    push!(terms,          (C(4ζ^2),                         a + 2))
    terms
end

@doc raw"""
Pure-Julia positive Laplacian form of two unnormalized Gaussian-polynomial
primitives. Drop-in replacement for `GaIn.laplacian`, matching the sign
convention required by `Hˢ_overlap(Ms; s=1)` (which adds this to the L²
overlap to build the H¹ overlap):

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
        for (c, e) in _axis_laplacian_terms(n2[axis], ζ2, R2[axis])
            m = ntuple(i -> i == axis ? e : n2[i], 3)
            acc += c * overlap_julia(ζ1, R1, nx1, ny1, nz1,
                                     ζ2, R2, m[1], m[2], m[3])
        end
    end
    -acc    # ⟨∇X₁,∇X₂⟩ = -⟨X₁,∇²X₂⟩
end

@doc raw"""
Kinetic-energy form of two unnormalized Gaussian-polynomial primitives,

```math
\langle X_1, -\tfrac12\nabla^2 X_2\rangle = \tfrac12\,\texttt{laplacian\_julia}(X_1, X_2).
```

Drop-in replacement for `GaIn.kinetic`.
"""
function kinetic_julia(ζ1, R1::AbstractVector, nx1::Int, ny1::Int, nz1::Int,
                       ζ2, R2::AbstractVector, nx2::Int, ny2::Int, nz2::Int)
    laplacian_julia(ζ1, R1, nx1, ny1, nz1, ζ2, R2, nx2, ny2, nz2) / 2
end
