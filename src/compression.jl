"""
Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
on an optimized symmetry adapted AO basis.
The basis is constructed by greedy iterations described in Ref[].
"""
function compress_graphene_pz_wannier(Wc::CompressedWannier, π_bond_axis; 
                  ### Kwargs ...
                  callback=default_callback(),
                  ## Outer loop parameters
                  max_xy_order=3, # Maximal order in (x,y) plane for the SAGTO polynomial
                  max_z_order=3,  # Maximal order in z axis for the SAGTO polynomial
                  ## Inner loop parameters
                  ζ_min=1e-3,     # Minimal possible spread to avoid NAN
                  max_iter=20,    # Maximal number of compression iterations
                  tol=1e-2,       # maximum attaigned for now is ≈ 5%
                  optim_method=ConjugateGradient(),
                  optim_options=Optim.Options(g_abstol=5e-1, show_trace=true),
                  )

    # Extract compressed wannier initial data
    wannier = Wc.wannier
    basis_SC = Wc.basis_supercell
    center = Wc.center
    residual = Wc.residual
    error = Wc.error
    s = Wc.error_norm
    residual_norm = Hs_norm(basis_SC, residual; s)
    MOs = Wc.MOs
    coefficients = Wc.coefficients

    @assert  size(wannier,2)==1 "The wannier function is to be"*
        " given as a single supercell vector"

    # Populate info and print header
    converged = false
    n_iter = zero(Int64)
    info = (; Wc, converged, n_iter)
    callback(info)

    # Objective function to optimize at each iteration
    function f(spreads; center, xy_orders, z_orders, in_linesearch=true)
        # Hack to avoid conditioning problems with small spreads
        findmin(spreads)[1] < ζ_min && return Inf
        
        # Compute optimal basis function for given spreads and residual (in Wc)
        SAGTOs = SAGTO_basis(center, spreads, xy_orders, z_orders)
        Φ, error = optimal_basis_function(Wc, SAGTOs; s)

        # Only return the Hs distance between Φ and the residual
        # for linesearch purposes
        (in_linesearch) && (return error)

        # Otherwise return Φ
        Φ, error
    end

    while ( (info.n_iter < max_iter) && !info.converged )
        n_iter += 1;
        
        # Inner optimization: find new optimal MO
        D3_sym = iszero(center - Wc.center)
        xy_orders, z_orders = select_orders(max_xy_order, max_z_order, D3_sym)
        n_SAGTOs = length(xy_orders) * length(z_orders)

        spreads_init = ones(n_SAGTOs) .* 1/2        
        res = NaN

        # Optimize w.r. to spreads or center and spreads if non centered SAGTOs.
        if D3_sym
            res = optimize(spreads -> f(spreads; Wc.center, xy_orders, z_orders),
                           spreads_init,
                           optim_method,
                           optim_options,
                           autodiff=:forward) # make it backward...
        else
            res = optimize(center_and_spreads -> f(center_and_spreads[3:end];
                                                  center=center_and_spreads[1:2],
                                                   xy_orders, z_orders),
                           [polar_to_cartesian_coords(Wc.center,
                                                      π_bond_axis...)
                            , ones(n_SAGTOs)...],
                           optim_method,
                           optim_options,
                           autodiff=:forward) # make it backward...
        end

        Φ_opti, _ = f(res.minimizer; center, xy_orders, z_orders, in_linesearch=false)

        # Add new MO to the MO basis and project Wn
        # Compute new res and check CV
        push!(MOs, Φ_opti)
        compressed_wannier, coeffs = project_wannier_on_basis(Wc; s)
        residual= wannier .- compressed_wannier
        error = Hs_norm(basis_SC, res; s) / Hs_norm(basis_SC, wannier; s)
        (error < tol) && (converged = true)

        # Actualize compressed wannier structure and info
        Wc.MOs = MOs
        Wc.coefficients = coeffs
        Wc.residual = residual
        Wc.error = error

        # Change symmetry requirement to explore both D3_sym and non D3_sym part of Wn
        (n_iter ≥ 1) && (center = iszero(center) ? rand(3) : zeros(3))

        # Actualize info
        info = merge(info, (;Wc, converged, n_iter))
        callback(info)
    end
    info 
end
