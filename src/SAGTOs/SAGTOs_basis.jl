using ThreadsX

@doc raw"""
Return the table of Symmetry Adatped Gaussian polynomial with same center and
respective spread in ζs, where the polynomial parts orders are given by xy_orders 
and z_orders.
"""
function SAGTOs_basis(α, ζs, xy_orders, z_orders)
    @assert( length(ζs) == length(xy_orders)* length(z_orders) )
    # Import all polynomial parameters in a single table and create corresponding AOs
    pol_parts = SAGTOs_polynoms(xy_orders, z_orders)
    [GaussianPolynomial(pol, α, ζ) for (pol, ζ) in zip(pol_parts, ζs)]
end

"""
Returns the table of Bloch decomposition of SAGTOs given by the above routine.
"""
function SAGTOs_basis(basis_SC::PlaneWaveBasis, α, ζs, xy_orders, z_orders)
    Χs = SAGTOs_basis(α, ζs, xy_orders, z_orders)
    ThreadsX.map(Χ->Χ(basis_SC), Χs)
end
