using Test
using Random
using Wannier2GTO

@testset "init_params: determinism" begin
    layout = JointLayout(N_centered=3, N_pibond=2)
    flat1 = init_params(MersenneTwister(42), layout)
    flat2 = init_params(MersenneTwister(42), layout)
    @test flat1 == flat2
    flat3 = init_params(MersenneTwister(43), layout)
    @test flat1 != flat3
end

@testset "init_params: dimension and finiteness" begin
    layout = JointLayout(N_centered=4, N_pibond=3)
    flat = init_params(MersenneTwister(1), layout)
    @test length(flat) == flat_dim(layout)
    @test all(isfinite, flat)
end

@testset "init_params: r bounds" begin
    layout = JointLayout(N_centered=2, N_pibond=4)
    flat = init_params(MersenneTwister(7), layout; r_min=0.5, r_max=5.0)
    for j in 1:layout.N_pibond
        r, _, _ = get_pibond(flat, layout, j)
        @test 0.5 ≤ r ≤ 5.0
    end
end
