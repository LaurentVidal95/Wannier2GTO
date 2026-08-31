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

@testset "evaluate_hoppings + hopping_criteria" begin
    ζ = 0.7
    g = W2G.GaussianPolynomial([(0, 0, 0)], [1.0], zeros(3), ζ)
    Φs = W2G.BasisFunction[W2G.BasisFunction([1.0], [g])]
    c = [1.0]
    R1 = [1.0, 0.0, 0.0]; R2 = [0.0, 0.0, 1.5]
    h1 = W2G.gto_hoppings(Φs, c, R1); h2 = W2G.gto_hoppings(Φs, c, R2)
    # Targets: entry 1 = training "intra a1" met exactly; entry 2 = validation
    # with a deliberately wrong reference sign.
    targets = W2G.HoppingTargets(["intra a1", "val x"], [:training, :validation],
                                 [R1, R2], [h1.S, h2.S], [h1.T, -h2.T],
                                 1.5, [1])
    rows = W2G.evaluate_hoppings(Φs, c, targets)
    @test length(rows) == 2
    @test rows[1].rel_err_T ≈ 0.0  atol = 1e-12
    @test rows[1].sign_ok
    @test !rows[2].sign_ok                      # wrong sign detected

    crit = W2G.hopping_criteria(rows, targets)
    @test length(crit) == 3
    c1 = crit[findfirst(x -> x.name == "relerr_T(a1) ≤ 5%", crit)]
    @test c1.pass
    @test c1.value ≈ 0.0  atol = 1e-12
    c2 = crit[findfirst(x -> x.name == "validation signs", crit)]
    @test !c2.pass                              # |T_ref| = h2.T > floor, sign wrong
    c3 = crit[findfirst(x -> x.name == "|S(a1)| ≤ 1e-3", crit)]
    @test c3.value ≈ abs(h1.S)  rtol = 1e-10    # here S(a1) is large → FAIL
    @test !c3.pass
end
