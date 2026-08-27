using Test
using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra

@testset "laplacian_julia: closed form (normalized s-Gaussian)" begin
    # g(r) = exp(-ζ r²), normalized to L² norm 1 by the GaussianPolynomial ctor.
    # ⟨∇g, ∇g⟩ = -⟨g, ∇²g⟩ = 3ζ ;  ⟨g, -½∇² g⟩ = 3ζ/2.
    for ζ in (0.3, 1.0, 2.5)
        g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)  # normalized
        @test W2G.integral(g, g; type=:laplacian) ≈ 3ζ      rtol = 1e-10
        @test W2G.integral(g, g; type=:kinetic)   ≈ 3ζ / 2  rtol = 1e-10
    end
end
