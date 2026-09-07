using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra
using ForwardDiff

function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "gto_hoppings: closed forms" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    for R in ([1.0, 0.0, 0.0], [0.5, -0.3, 0.2])
        R² = sum(abs2, R)
        h = W2G.gto_hoppings(Φs, c, R)
        @test h.S ≈ exp(-ζ * R² / 2)                     rtol = 1e-10
        @test h.T ≈ (ζ/2) * (3 - ζ*R²) * exp(-ζ*R²/2)   rtol = 1e-10
    end
end

# Synthetic targets around a single training R + one ortho index.
function synthetic_targets(; S1, T1, T0)
    W2G.HoppingTargets(["intra a1"], [:training], [[1.0, 0.0, 0.0]],
                       [S1], [T1], T0, [1])
end

@testset "hopping_penalty: zero at target, positive off target" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    h = W2G.gto_hoppings(Φs, c, [1.0, 0.0, 0.0])   # norm² = 1 here
    t_exact = synthetic_targets(S1=h.S, T1=h.T, T0=1.5)
    @test W2G.hopping_penalty(Φs, c, t_exact; μ=1.0, ν=0.0) ≈ 0.0  atol = 1e-14
    t_off = synthetic_targets(S1=h.S + 0.1, T1=h.T, T0=1.5)
    @test W2G.hopping_penalty(Φs, c, t_off; μ=1.0, ν=0.0) ≈ 0.1^2 / 1.5^2  rtol = 1e-10
    # ν penalty adds S² for ortho-indexed entries.
    @test W2G.hopping_penalty(Φs, c, t_exact; μ=0.0, ν=1.0) ≈ h.S^2  rtol = 1e-10
end

@testset "joint_loss: targets=nothing ≡ μ=0, gradient finite" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...)))
    w = w / norm(w)
    layout = W2G.JointLayout(N_centered=1, N_pibond=1)
    flat = W2G.init_params(W2G.MersenneTwister(11), layout)
    kw = (log_ζ_min=log(1e-2), log_ζ_max=log(4.0), π_bond_unit=[1.0, 0.0, 0.0])

    L0 = W2G.joint_loss(flat, layout, w, basis; kw...)
    targets = synthetic_targets(S1=0.02, T1=0.01, T0=1.5)
    @test W2G.joint_loss(flat, layout, w, basis; kw..., targets, μ=0.0, ν=0.0) ≈ L0  rtol = 1e-12
    Lμ = W2G.joint_loss(flat, layout, w, basis; kw..., targets, μ=10.0, ν=10.0)
    @test Lμ > L0

    f = θ -> W2G.joint_loss(θ, layout, w, basis; kw..., targets, μ=10.0, ν=10.0)
    G = ForwardDiff.gradient(f, flat)
    @test all(isfinite, G)
    G0 = ForwardDiff.gradient(θ -> W2G.joint_loss(θ, layout, w, basis; kw...), flat)
    @test !(G ≈ G0)   # the penalty actually contributes to the gradient
end

@testset "joint_loss with penalty: ForwardDiff gradient vs finite differences" begin
    # Regression for the `filter_dual` normalization bug (Sept 2026): the H¹
    # gradient was FD-validated in phase B, but the analytic-integral path of
    # the penalty reads normalized SAGTO coefficients directly, so its ∂/∂ζ
    # must be checked separately. Spread parameters are the ones that broke.
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...)))
    w = w / norm(w)
    layout = W2G.JointLayout(N_centered=1, N_pibond=1)
    flat = W2G.init_params(W2G.MersenneTwister(11), layout)
    kw = (log_ζ_min=log(1e-2), log_ζ_max=log(4.0), π_bond_unit=[1.0, 0.0, 0.0])
    targets = W2G.HoppingTargets(["intra a1", "inter AA"], [:training, :training],
                                 [[1.0, 0.0, 0.0], [0.3, 0.2, 1.2]],
                                 [0.0, -0.05], [0.01, -0.001], 1.5, [1])
    f = θ -> W2G.joint_loss(θ, layout, w, basis; kw..., targets, μ=10.0, ν=10.0)
    G = ForwardDiff.gradient(f, flat)
    FD_STEP = 1e-5
    for i in eachindex(flat)
        e = zeros(length(flat)); e[i] = FD_STEP
        fd = (f(flat + e) - f(flat - e)) / (2 * FD_STEP)
        @test isapprox(G[i], fd; rtol=1e-5, atol=1e-8)
    end
end
