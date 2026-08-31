using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra

# Synthetic compressed "Wannier": one normalized s-Gaussian at the origin.
# Closed forms for two normalized s-Gaussians with equal spread ζ, separation R
# (Szabo–Ostrund A.9/A.11 with α=β=ζ):
#   S(R) = exp(-ζ|R|²/2)
#   T(R) = ⟨g|-½∇²|g_R⟩ = (ζ/2)(3 - ζ|R|²) exp(-ζ|R|²/2)
function single_gaussian_Wc(ζ)
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)  # normalized
    Φ = W2G.BasisFunction([1.0], [g])
    W2G.CompressedWannier(zeros(3), W2G.BasisFunction[Φ], [1.0],
                          nothing, ComplexF64[], ComplexF64[], 0.0, 0)
end

# Minimal Γ-only plane-wave basis (kinetic-only model, no atoms needed).
function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr, anisotropic on purpose
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "compare_hopping: GTO side against closed forms" begin
    ζ = 0.7
    Wc = single_gaussian_Wc(ζ)
    basis = tiny_basis()
    # Any admissible w_fourier does for this testset: take the Gaussian itself
    # expanded on the basis, so the reference side runs without erroring.
    w = Wc(basis)
    w = w / norm(w)
    for R in ([1.0, 0.0, 0.0], [0.5, -0.3, 0.2], [0.0, 0.0, 1.3])
        out = W2G.compare_hopping(Wc, w, basis, R)
        R² = sum(abs2, R)
        @test out.S_gto ≈ exp(-ζ * R² / 2)                      rtol = 1e-10
        @test out.T_gto ≈ (ζ/2) * (3 - ζ*R²) * exp(-ζ*R²/2)    rtol = 1e-10
    end
end

@testset "compare_hopping: reference side identities" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    n_G = length(G_vectors(basis, kpt))
    # Real random real-space field -> Fourier coefficients with w_{-G} = conj(w_G),
    # i.e. a real "Wannier", as in production.
    w_real = randn(basis.fft_size...)
    w = DFTK.fft(basis, kpt, complex.(w_real))
    w = w / norm(w)
    @assert length(w) == n_G
    Wc = single_gaussian_Wc(0.7)

    out0 = W2G.compare_onsite_kinetic(Wc, w, basis)
    # R = 0 must reproduce the on-site diagnostics exactly.
    outR0 = W2G.compare_hopping(Wc, w, basis, zeros(3))
    @test outR0.S_ref ≈ 1.0            rtol = 1e-12
    @test outR0.T_ref ≈ out0.T_ref     rtol = 1e-12
    @test outR0.T_gto ≈ out0.T_gto     rtol = 1e-12

    # Translation by a lattice vector is the identity on the periodic reference.
    L = basis.model.lattice[:, 1]
    outL = W2G.compare_hopping(Wc, w, basis, L)
    @test outL.S_ref ≈ outR0.S_ref     rtol = 1e-10
    @test outL.T_ref ≈ outR0.T_ref     rtol = 1e-10

    # Hermiticity: S(R) = S(-R), T(R) = T(-R) for a real w.
    R = [1.1, -0.4, 0.6]
    outp = W2G.compare_hopping(Wc, w, basis, R)
    outm = W2G.compare_hopping(Wc, w, basis, -R)
    @test outp.S_ref ≈ outm.S_ref      rtol = 1e-9
    @test outp.T_ref ≈ outm.T_ref      rtol = 1e-9
end

@testset "translate(BasisFunction) and reference_hopping" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φ = W2G.BasisFunction([1.0], [g])
    R = [0.4, -0.2, 0.9]
    ΦR = W2G.translate(Φ, R)
    # Same closed form as the CompressedWannier path.
    @test W2G.integral(Φ, ΦR; type=:overlap) ≈ exp(-ζ * sum(abs2, R) / 2)  rtol = 1e-10

    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w_real = randn(basis.fft_size...)
    w = DFTK.fft(basis, kpt, complex.(w_real))
    w = w / norm(w)
    out0 = W2G.reference_hopping(w, basis, zeros(3))
    @test out0.S_ref ≈ 1.0  rtol = 1e-12
    @test out0.T_ref ≈ out0.T0_ref  rtol = 1e-12
    # Must agree with compare_hopping's reference side.
    Wc = single_gaussian_Wc(ζ)
    Rb = [1.1, -0.4, 0.6]
    @test W2G.reference_hopping(w, basis, Rb).S_ref ≈ W2G.compare_hopping(Wc, w, basis, Rb).S_ref  rtol = 1e-12
    @test W2G.reference_hopping(w, basis, Rb).T_ref ≈ W2G.compare_hopping(Wc, w, basis, Rb).T_ref  rtol = 1e-12
end
