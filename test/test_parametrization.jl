using Test
using Wannier2GTO

@testset "JointLayout" begin
    layout = JointLayout(N_centered=3, N_pibond=2,
                         poly_dim_centered=6, poly_dim_pibond=4)
    # 3*(1+6) + 2*(1+1+4) = 21 + 12 = 33
    @test flat_dim(layout) == 33
end

@testset "flat ↔ structured" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = collect(1.0:flat_dim(layout))   # 2*(1+6) + 1*(1+1+4) = 20

    # Centered Φ_1: u = 1.0, λ = [2, 3, 4, 5, 6, 7]
    u_c, λ_c = get_centered(flat, layout, 1)
    @test u_c == 1.0
    @test λ_c == [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    # Centered Φ_2: u = 8.0, λ = [9, 10, 11, 12, 13, 14]
    u_c2, _ = get_centered(flat, layout, 2)
    @test u_c2 == 8.0

    # π-bond Φ_1: r = 15.0, u = 16.0, λ = [17, 18, 19, 20]
    r_p, u_p, λ_p = get_pibond(flat, layout, 1)
    @test r_p == 15.0 && u_p == 16.0
    @test λ_p == [17.0, 18.0, 19.0, 20.0]
end
