import LinearAlgebra.normalize

struct BasisFunction{T<:Real}
    coeffs :: Vector{T}
    SAGTOs # Vector of Gaussian polynomials
end

function (Φ::BasisFunction)(basis_supercell::PlaneWaveBasis; normalize_SAGTO=true)
    SAGTOs_Fourier = [X(basis_supercell; normalize_SAGTO) for X in Φ.SAGTOs]
    sum(λ * SAGTOs_Fourier[i] for  (i, λ) in enumerate(Φ.coeffs))
end

function enforce_D3_symmetry(Φ::BasisFunction)
    SAGTOs = Φ.SAGTOs
    center = SAGTOs[1].center

    # Add SAGTOs with rotated centers at +2π/3 and +4π/3
    SAGTOs_ω = [GaussianPolynomial(X, rot(2π/3)*center) for X in SAGTOs]
    SAGTOs_ω_square = [GaussianPolynomial(X, rot(4π/3)*center) for X in SAGTOs]
    # Return new basis function
    new_coeffs = vcat([Φ.coeffs for _ in 1:3]...)
    BasisFunction(new_coeffs, vcat(SAGTOs, SAGTOs_ω, SAGTOs_ω_square))
end

@doc raw"""
Build a closure `(S, Γ) -> coeffs` that solves ``S\,c = \Gamma`` with the chosen
regularization. Available strategies:

- `:none`         — plain `S \ Γ` (no regularization).
- `:tikhonov`     — solve ``(S + \varepsilon I)\,c = \Gamma`` with
  ``\varepsilon = \texttt{param}\cdot \mathrm{tr}(S)/n``. Default `param=1e-8`.
- `:svd_truncation` — eigendecompose ``S = V\Lambda V^*``, drop directions with
  ``\lambda < \texttt{param}\cdot \lambda_\max``, solve in the kept subspace and
  re-expand. Default `param=1e-10`.
- `:pivoted_cholesky` — pivoted Cholesky with rank revelation
  (`tol = param`, default `1e-10`). Truncates the rank to discard quasi-null
  directions, then back-solves in the kept subspace.

The returned closure is meant to replace `S \ Γ` in `optimal_basis_function`.
All strategies preserve the shape and scaling of the result; only the
treatment of near-null singular directions differs.
"""
function _make_inner_solver(reg::Symbol, param=nothing)
    if reg === :none
        return (S, Γ) -> S \ Γ

    elseif reg === :tikhonov
        ε_rel = something(param, 1e-8)
        return function (S, Γ)
            n = size(S, 1)
            ε = ε_rel * real(tr(S)) / n
            (S + ε * I) \ Γ
        end

    elseif reg === :svd_truncation
        rel_tol = something(param, 1e-10)
        return function (S, Γ)
            F = eigen(Hermitian(Matrix(S)))
            λ_max = maximum(F.values)
            keep = F.values .> rel_tol * λ_max
            Vk = F.vectors[:, keep]
            λk = F.values[keep]
            Vk * (Diagonal(1 ./ λk) * (Vk' * Γ))
        end

    elseif reg === :pivoted_cholesky
        rel_tol = something(param, 1e-10)
        return function (S, Γ)
            F = cholesky(Hermitian(Matrix(S)), RowMaximum(); check=false,
                         tol=rel_tol * real(tr(S)) / size(S, 1))
            r = F.rank
            P = F.p[1:r]                      # kept pivot indices
            L = F.L[1:r, 1:r]                 # r×r lower-triangular factor
            # Solve L Lᵀ y = Γ_P, then scatter back
            y = L' \ (L \ Γ[P])
            x = zero(Γ)
            x[P] = y
            x
        end

    else
        error("Unknown regularization strategy: $reg. Expected one of: " *
              ":none, :tikhonov, :svd_truncation, :pivoted_cholesky.")
    end
end

"""
Construct a linear combination of given SAGTOs that best
approximates the residual contained in the CompressedWannier struct.

`solve_S` is a callable `(S, Γ) -> coeffs`; defaults to plain `S \\ Γ`.
Use `_make_inner_solver(:tikhonov, ε_rel)` etc. to inject a regularized solve.
"""
function optimal_basis_function(Wc, SAGTOs; tol=1e-5,
                                solve_S=(S, Γ) -> S \ Γ)
    # Extract needed data
    s = Wc.error_norm
    basis_supercell = Wc.basis_supercell
    center = SAGTOs[1].center

    # Compute optimal coefficients for given SAGTOs
    SAGTOs_Four = [X(basis_supercell) for X in SAGTOs]
    Γ = ThreadsX.map(X -> Hˢ_dot(basis_supercell, Wc.residual, X; s), SAGTOs_Four)
    S = Hˢ_overlap(basis_supercell, SAGTOs_Four; s)
    raw_coeffs = solve_S(S, Γ)
    # Phase A diagnostic: log inner cond(S) + |imag(coeffs)|. After Tikhonov or
    # SVD truncation the imaginary part should stay around machine epsilon;
    # `safereal` then enforces this as a contract.
    let σ = svdvals(Hermitian(Matrix(S)))
        @printf("    [inner] n_AOs=%2d  cond(S)=%.3e  σ=[%.2e,%.2e]  |imag(coeffs)|=%.3e\n",
                length(SAGTOs), σ[1] / σ[end], σ[end], σ[1], norm(imag(raw_coeffs)))
    end
    optimal_coeffs = safereal(filter_small_coeffs.(raw_coeffs; tol))

    # Assemble optimal basis function and enforce D3 symmetry if needed
    Φ = BasisFunction(optimal_coeffs, SAGTOs)
    !iszero(center - Wc.center) && (Φ = enforce_D3_symmetry(Φ))
    Φ, Hˢ_norm(basis_supercell, Φ(basis_supercell) - Wc.residual; s=Wc.error_norm)^2
end

function normalize(Φ::BasisFunction)
    norm_Φ = √(integral(Φ, Φ; type=:overlap))
    normalized_coeffs = Φ.coeffs ./ norm_Φ
    BasisFunction(normalized_coeffs, Φ.SAGTOs)
end
