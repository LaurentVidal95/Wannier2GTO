"""
Routines to build an SAGTO adapted to the symmetries of pz-type Wannier functions of graphene.
    • If AOs are centered on the Wannier center, the chosen symmetry is D3 in the (x,y) plane.
    • In the other case, only the symmetry w.r. to the (x,y) and (x,z) planes are retained.
    (i.e. order 0 and 2 for xy orders)
"""

"""
Generate polynom parameters so that the associated AOs have the same symmetry
as the pz-like Wannier function of graphene.
"""
function symmetry_adapted_polynoms(nxy::Int64, nz::Int64)
    @assert( nxy ∈ (0,2,3,4,6,9) )
    (nxy==0) && (return [(0,0,nz)], [1.])
    (nxy==2) && (return [(2,0,nz), (0,2,nz)], [1., 1.])
    (nxy==3) && (return [(3,0,nz), (1,2,nz)], [1.,-3.])
    (nxy==4) && (return [(4,0,nz), (2,2,nz), (0,4,nz)], [1.,2.,1.])
    (nxy==6) && (return [(6,0,nz), (4,2,nz), (2,4,nz), (0,6,nz)], [1.,-15.,15.,-1.])
    (nxy==9) && (return [(9,0,nz), (7,2,nz), (5,4,nz), (3,6,nz), (1,8,nz)],
                      [1.,-36., 126., -84., 9.])
    nothing
end

"""
Select among above parameters given a maximum order in (x,y) and z.
"""
function select_orders(max_xy_order, max_z_order, D3_sym::Bool)
    !(D3_sym) && (max_xy_order=min(max_xy_order,2))
    xy_orders = [0,2,3,6,9]
    z_orders = [1,3,5,7,9]
    xy_orders[xy_orders .≤ max_xy_order], z_orders[z_orders .≤ max_z_order]
end

# """
# Attention: the GaIn lib only handles gaussian integrals when (nx+ny+nz)≤6
# """
@inline function symmetry_adapted_polynoms(xy_orders, z_orders)
    [symmetry_adapted_polynoms(nxy, nz) for nz in z_orders
     for nxy in xy_orders if (nz+nxy ≤ 6)]
end

"""
Wraps the two functions above
"""
function symmetry_adapted_polynoms(max_xy_order::Int64, max_z_order::Int64, D3_sym::Bool)
    xy_orders, z_orders = select_orders(max_xy_order, max_z_order, D3_sym)
    symmetry_adapted_polynoms(xy_orders, z_orders)
end

# function ∇SAGTO_parameters(xy_orders, z_orders)
#     SAGTO_params = SAGTO_parameters(xy_orders, z_orders)
#     ∇_SAGTO_params = []
#     for i in 1:3
#         update = zeros(Int64, 3); update[i] = 2;
#         raise_powers(param) = map(powers -> powers .+ update, param)
#         for params in SAGTO_params
#             push!(∇_SAGTO_params, [raise_powers(params[1]), .- params[2]])
#         end
#     end
#     ∇_SAGTO_params
# end

# function ∇SAGTOs_graphene_pz_wannier(basis::PlaneWaveBasis, α::Vector{T}, ζs::Vector{T},
#                                        xy_orders, z_orders) where {T<:Real}
#     ∇SAGTO_params = ∇SAGTO_parameters(xy_orders, z_orders)
#     N_∇Χ = length(ζs)
#     # Compute x^2, y^2 and z^2 terms as separate SAGTOs
#     AOs = [GaussianPolynomial(Iα, λα, α, ζ) for ((Iα, λα), ζ) in zip(∇SAGTO_params,
#                                                                      vcat(ζs, ζs, ζs))]
#     AOs_fourier = Bloch_decomposition.(basis, AOs)
#     AOs[1:N_∇Χ] .+ AOs[N_∇Χ+1:2*N_∇Χ] .+ AOs[2*N_∇Χ+1:end]
#     Bloch_decomposition.(basis, AOs[1:N_∇Χ] .+ AOs[N_∇Χ+1:2*N_∇Χ] .+ AOs[2*N_∇Χ+1:end])
#     # AOs_fourier[1:N_∇Χ] .+ AOs_fourier[N_∇Χ+1:2*N_∇Χ] .+ AOs_fourier[2*N_∇Χ+1:end]
# end
