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

@testset "log ζ encoding round-trip" begin
    log_ζ_min = log(1e-2)
    log_ζ_max = log(4.0)
    for ζ in (0.05, 0.5, 1.0, 3.0, 3.99)
        log_ζ = log(ζ)
        u = log_ζ_to_u(log_ζ, log_ζ_min, log_ζ_max)
        log_ζ_back = log_ζ_from_u(u, log_ζ_min, log_ζ_max)
        @test log_ζ_back ≈ log_ζ rtol=1e-12
    end
    # Out-of-bounds u stays in the box
    log_ζ_huge = log_ζ_from_u(100.0, log_ζ_min, log_ζ_max)
    @test log_ζ_huge ≈ log_ζ_max rtol=1e-10
    log_ζ_neg = log_ζ_from_u(-100.0, log_ζ_min, log_ζ_max)
    @test log_ζ_neg ≈ log_ζ_min rtol=1e-10
end

@testset "params_to_basis_functions: shape" begin
    layout = JointLayout(N_centered=2, N_pibond=1)
    flat = zeros(flat_dim(layout))    # all zeros
    log_ζ_min, log_ζ_max = log(1e-2), log(4.0)
    π_bond_dir = [cos(-π/6), sin(-π/6), 0.0]

    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_dir,
                                    max_xy_order=3, max_z_order=3)
    @test length(Φs) == 3
    # Centered Φ_1, Φ_2: 6 SAGTOs each
    @test length(Φs[1].SAGTOs) == 6
    @test length(Φs[2].SAGTOs) == 6
    # π-bond Φ_3: 4 SAGTOs × 3 D3 rotations = 12 SAGTOs
    @test length(Φs[3].SAGTOs) == 12
end

@testset "params_to_basis_functions: centers and spreads" begin
    layout = JointLayout(N_centered=1, N_pibond=1)
    flat = zeros(flat_dim(layout))
    flat[1] = 0.0                  # centered u = 0 -> log ζ at midpoint
    flat[8] = 2.0                  # π-bond r = 2.0
    flat[9] = 1.0                  # π-bond u = 1.0
    log_ζ_min, log_ζ_max = log(1e-2), log(4.0)
    π_bond_dir = [cos(-π/6), sin(-π/6), 0.0]

    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_dir,
                                    max_xy_order=3, max_z_order=3)
    # Φ_1 centered at origin
    @test Φs[1].SAGTOs[1].center ≈ zeros(3)
    # Φ_2 first SAGTO at r * π_bond_dir, the others at the rotated positions
    @test Φs[2].SAGTOs[1].center ≈ 2.0 .* π_bond_dir rtol=1e-12
    # Spread reflects sigmoid(u=1.0) ≈ 0.731
    expected_log_ζ = log(1e-2) + (log(4.0) - log(1e-2)) * (1 / (1 + exp(-1.0)))
    @test log(Φs[2].SAGTOs[1].spread) ≈ expected_log_ζ rtol=1e-10
end
