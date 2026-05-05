using Test
using Random
using DFTK
using DFTK.Unitful
using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra

# Heavy fixture: bring `Graphene` into scope and build the supercell basis
# once for all testsets in this file.
include(joinpath(@__DIR__, "..", "workflow", "monolayer_graphene.jl"))

const _BASIS_SC = let
    basis = Graphene(; d=10u"Å", kgrid=[5,5,1], Ecut=15).basis()
    DFTK.cell_to_supercell(basis)
end

const _DATA = W2G.read_wannier_function(joinpath(@__DIR__, "..",
                            "workflow/wannier_functions/wannier_pz_Ecut-15.json"))
const _W_FOURIER = normalize(_DATA.wannier)
const _π_BOND_DIR = let
    θ = _DATA.π_bond.θ
    [cos(θ), sin(θ), 0.0]
end

@testset "joint_loss: positive and finite" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = init_params(MersenneTwister(0), layout)
    L = joint_loss(flat, layout, _W_FOURIER, _BASIS_SC;
                   log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                   π_bond_unit=_π_BOND_DIR)
    @test isfinite(L)
    @test L > 0
end

@testset "joint_loss: different params → different losses" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat1 = init_params(MersenneTwister(0), layout)
    flat2 = init_params(MersenneTwister(1), layout)
    L1 = joint_loss(flat1, layout, _W_FOURIER, _BASIS_SC;
                    log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                    π_bond_unit=_π_BOND_DIR)
    L2 = joint_loss(flat2, layout, _W_FOURIER, _BASIS_SC;
                    log_ζ_min=log(1e-2), log_ζ_max=log(4.0),
                    π_bond_unit=_π_BOND_DIR)
    @test L1 != L2
end
