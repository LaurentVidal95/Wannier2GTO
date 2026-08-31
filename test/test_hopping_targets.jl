using Test
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using LinearAlgebra

function tiny_basis()
    lattice = diagm([8.0, 9.0, 10.0])  # Bohr
    model = Model(lattice; terms=[Kinetic()], n_electrons=1,
                  spin_polarization=:spinless, symmetries=false)
    PlaneWaveBasis(model; Ecut=8, kgrid=(1, 1, 1))
end

@testset "hopping_R_sets: determinism and geometry" begin
    a₁ = [4.99, 0.0, 0.0]; a₂ = [-2.49, 4.32, 0.0]
    kw = (d_inter=6.33, z_val=(5.67, 6.99), seed=1234, n_random=5,
          r_min=3.78, r_max=7.56)   # Bohr ≈ (3.35, 3.0, 3.7, 2.0, 4.0) Å
    rs1 = W2G.hopping_R_sets(a₁, a₂; kw...)
    rs2 = W2G.hopping_R_sets(a₁, a₂; kw...)
    @test rs1.labels == rs2.labels
    @test all(rs1.Rs .≈ rs2.Rs)                      # same seed → identical
    @test count(==(:training), rs1.sets) == 5
    @test count(==(:validation), rs1.sets) == 4 + 5   # deterministic core + random
    @test rs1.ortho_idx == findall(l -> startswith(l, "intra"), rs1.labels)
    for i in eachindex(rs1.labels)
        startswith(rs1.labels[i], "rand") || continue
        @test kw.r_min ≤ norm(rs1.Rs[i]) ≤ kw.r_max
    end
    rs3 = W2G.hopping_R_sets(a₁, a₂; kw..., seed=4321)
    @test !all(rs1.Rs .≈ rs3.Rs)                     # seed changes randoms
end

@testset "build_hopping_targets + JSON round-trip" begin
    basis = tiny_basis()
    kpt = only(basis.kpoints)
    w = DFTK.fft(basis, kpt, complex.(randn(basis.fft_size...)))
    w = w / norm(w)
    a₁ = basis.model.lattice[:, 1]; a₂ = basis.model.lattice[:, 2]
    rsets = W2G.hopping_R_sets(a₁, a₂; d_inter=6.33, z_val=(5.67, 6.99),
                               seed=7, n_random=2, r_min=3.78, r_max=7.56)
    targets = W2G.build_hopping_targets(w, basis, rsets)
    @test targets.T0_ref ≈ W2G.reference_hopping(w, basis, zeros(3)).T0_ref
    i1 = findfirst(==("intra a1"), targets.labels)
    @test targets.S_ref[i1] ≈ W2G.reference_hopping(w, basis, a₁).S_ref  rtol = 1e-12

    file = joinpath(mktempdir(), "targets.json")
    W2G.store(targets; file)
    loaded = W2G.HoppingTargets(file)
    @test loaded.labels == targets.labels
    @test loaded.sets == targets.sets
    @test all(loaded.Rs .≈ targets.Rs)
    @test loaded.S_ref ≈ targets.S_ref
    @test loaded.T_ref ≈ targets.T_ref
    @test loaded.T0_ref ≈ targets.T0_ref
    @test loaded.ortho_idx == targets.ortho_idx
end
