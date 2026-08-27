using Test
using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra

@doc raw"""
Closed-form positive Laplacian form of a *normalized* Cartesian Gaussian
``g = N\,x^{n_x}y^{n_y}z^{n_z} e^{-\zeta r^2}``. Separating per axis and using
the 1D moment ratios ``M_{2a\pm2}/M_{2a}`` gives

```math
\langle \nabla g, \nabla g\rangle
    = \zeta \sum_{i \in \{x,y,z\}} \frac{4n_i - 1}{2n_i - 1}.
```

This reduces to the familiar ``(2L+3)\zeta`` only when every exponent is 0 or 1;
for ``z^2``-type Cartesians (which mix L=2 and L=0) it does not, e.g.
``(0,0,2) \mapsto \tfrac{13}{3}\zeta`` rather than ``5\zeta``.
"""
closed_form_laplacian(exps::NTuple{3,Int}, ζ) =
    ζ * sum((4n - 1) / (2n - 1) for n in exps)

@testset "laplacian_julia: closed form, normalized Cartesian Gaussians" begin
    # Covers L = 0, 1, 2 and, crucially, exponent 2 — which exercises the
    # a ≥ 2 branch of `_axis_laplacian_terms` that an s-only test never reaches.
    for exps in [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 0, 1), (0, 0, 2), (1, 1, 1)]
        for ζ in (0.3, 1.0, 1.7)
            g = W2G.GaussianPolynomial([exps], [1.0], zeros(3), ζ)  # normalized
            expected = closed_form_laplacian(exps, ζ)
            @test W2G.integral(g, g; type=:laplacian) ≈ expected     rtol = 1e-10
            @test W2G.integral(g, g; type=:kinetic)   ≈ expected / 2 rtol = 1e-10
        end
    end
end

# ---------------------------------------------------------------------------
# Grid cross-validation. Heavy fixture: a Γ-point supercell plane-wave basis,
# built once for the testsets below.
# ---------------------------------------------------------------------------
using DFTK
using DFTK.Unitful

include(joinpath(@__DIR__, "..", "workflow", "monolayer_graphene.jl"))

const _BASIS_SC = let
    basis = Graphene(; d=10u"Å", kgrid=[5, 5, 1], Ecut=15).basis()
    DFTK.cell_to_supercell(basis)
end

"""
``\\langle\\nabla X_1,\\nabla X_2\\rangle`` on the plane-wave grid: ``\\nabla^2``
is multiplication by ``-|G|^2`` in Fourier, so the positive form is
`dot(X₁_four, |G|² .* X₂_four)`.
"""
function _grid_laplacian(X1::W2G.GaussianPolynomial, X2::W2G.GaussianPolynomial)
    kpt = only(_BASIS_SC.kpoints)
    G2 = [sum(abs2, q) for q in DFTK.G_vectors_cart(_BASIS_SC, kpt)]
    real(dot(X1(_BASIS_SC), G2 .* X2(_BASIS_SC)))
end

# Non-degenerate pairs: pairs whose integral vanishes by parity (e.g. xz × yz
# about a shared y-center) are useless here, since a relative error on 0 ≈ 0
# tests nothing. Spreads stay moderate (ζ ≲ Ecut/4) so the Gaussians remain
# resolved on the grid — the phase-A numerical-escape pathology otherwise
# makes the plane-wave representation meaningless.
const _GRID_PAIRS = [
    ("s ζ=1.0 @0 / s ζ=0.5 @0",
     W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), 1.0),
     W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), 0.5)),
    ("pz ζ=1.0 @0 / pz ζ=0.5 @0",
     W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), 1.0),
     W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), 0.5)),
    ("pz ζ=1.2 @0 / pz ζ=0.7 off-center",
     W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), 1.2),
     W2G.GaussianPolynomial([(0, 0, 1)], [1.0], [0.3, 0.0, 0.0], 0.7)),
    ("dxz ζ=1.1 @0 / pz ζ=0.6 off-center",
     W2G.GaussianPolynomial([(1, 0, 1)], [1.0], zeros(3), 1.1),
     W2G.GaussianPolynomial([(0, 0, 1)], [1.0], [0.2, 0.1, 0.0], 0.6)),
    ("dz2 ζ=0.9 @0 / dz2 ζ=0.9 @0",      # exercises the a ≥ 2 branch on the grid
     W2G.GaussianPolynomial([(0, 0, 2)], [1.0], zeros(3), 0.9),
     W2G.GaussianPolynomial([(0, 0, 2)], [1.0], zeros(3), 0.9)),
]

@testset "laplacian_julia: analytic vs grid" begin
    # Measured relative errors at Ecut=15 range from 7e-7 to 1.4e-4; 1e-3 keeps
    # a safety margin without being vacuous.
    for (label, X1, X2) in _GRID_PAIRS
        @testset "$label" begin
            @test W2G.integral(X1, X2; type=:laplacian) ≈ _grid_laplacian(X1, X2) rtol = 1e-3
        end
    end
end

@testset "laplacian_julia: symmetry" begin
    for (label, X1, X2) in _GRID_PAIRS
        @testset "$label" begin
            @test W2G.integral(X1, X2; type=:laplacian) ≈
                  W2G.integral(X2, X1; type=:laplacian) rtol = 1e-12
        end
    end
end

@testset "compare_onsite_kinetic: single-Gaussian triple agreement" begin
    # A lone normalized s-Gaussian playing both roles: its own grid transform is
    # the "reference" Wannier, its analytic integral the "Gaussian" side. Grid,
    # analytic and closed form (3ζ/2) must all agree.
    ζ = 1.0
    X = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)   # normalized
    Φ = W2G.BasisFunction([1.0], [X])
    # NB: the field is typed `Vector{BasisFunction}` (unparameterized), so the
    # element type must be annotated — `[Φ]` would be a Vector{BasisFunction{Float64}}.
    Wc = W2G.CompressedWannier(zeros(3), W2G.BasisFunction[Φ], [1.0],
                               _BASIS_SC, ComplexF64[], ComplexF64[], 0.0, 0.0)
    w_fourier = normalize(X(_BASIS_SC))

    out = W2G.compare_onsite_kinetic(Wc, w_fourier, _BASIS_SC)
    @test out.T_gto ≈ 3ζ / 2 rtol = 1e-10          # closed form
    @test out.rel_err < 1e-3                        # grid vs analytic
    @test out.rel_err_per_norm < 1e-3
    @test out.norm_ref ≈ 1.0 rtol = 1e-12
    @test out.norm_gto ≈ 1.0 rtol = 1e-6
end

@testset "Hˢ_overlap(Ms; s=1) is GaIn-free and SPD" begin
    Ms = [W2G.GaussianPolynomial([(0, 0, 1)], [1.0], zeros(3), ζ) for ζ in (0.5, 1.0, 2.0)]
    S = W2G.Hˢ_overlap(Ms; s=1)
    @test all(isfinite, S)
    @test S ≈ S' rtol = 1e-12
    @test isposdef(Symmetric(S))
end
