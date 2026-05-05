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
