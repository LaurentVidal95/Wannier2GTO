using ReverseDiff

function objective_f_and_g!(Wc, ζ_min)
    # Objective function to optimize at each iteration. The function produce an optimal
    # basis function (which best approximate the current residual) as linear combination of
    # SAGTOs with given spreads and center, respecting the symmetries encoded in xy_orders and z_orders.
    # When "in_linesearch", the function only returns the Hˢ error of approximation of the
    # residual by the optimal basis function.
    function f(spreads; center, xy_orders, z_orders, in_linesearch=true)
        # Hack to avoid conditioning problems with small spreads
        findmin(spreads)[1] < ζ_min && return Inf
        
        # Compute optimal basis function for given spreads, center and residual (stored in Wc)
        SAGTOs = SAGTO_basis(center, spreads, xy_orders, z_orders)
        Φ, error = optimal_basis_function(Wc, SAGTOs)

        (in_linesearch) && (return error)
        Φ, error
    end

    # Define optimized autodiff gradient with ReverseDiff
    # n_AOs_max = prod(length.(select_orders(max_xy_order, max_z_order, true)))
    # g_tape = ReverseDiff.GradientTape(f, 

    function g!(G, input; D3_sym, kwargs...)
        foo = similar(G)
        if D3_sym
            foo = ReverseDiff.gradient(X->f(X; kwargs...), input)
        else
            foo = ReverseDiff.gradient(X->f(X[4:end]; center=X[1:3], kwargs...), input)
        end
        G .= foo
    end
    f, g!
end


"""
Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
on an optimized symmetry adapted basis.
The basis functions are constrained to have their centers somewhere on the π-bonds axes.
The basis is constructed by greedy iterations described in arXiv:1712.02996.
"""
function compress_graphene_pz_wannier(Wc::CompressedWannier, π_bond;
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
                  optim_options=Optim.Options(g_abstol=1e-4, show_trace=true),
                  )

    # Extract compressed wannier data
    wannier = Wc.wannier
    basis_SC = Wc.basis_supercell
    center = copy(Wc.center)
    π_bond_center = polar_to_cartesian_coords(Wc.center, π_bond.r, π_bond.θ)
    residual = Wc.residual
    error = Wc.error
    s = Wc.error_norm
    residual_norm = Hs_norm(basis_SC, residual; s)
    basis_functions = Wc.basis_functions
    coefficients = Wc.coefficients

    @assert  size(wannier,2)==1 "The wannier function is to be"*
        " given as a single supercell vector"

    # Populate info and print header
    converged = false
    n_iter = zero(Int64)
    info = (; Wc, converged, n_iter)
    callback(info)

    # Define objective function to optimize at each iteration
    f, g! = objective_f_and_g!(Wc, ζ_min)
    
    while ( (info.n_iter < max_iter) && !info.converged )
        n_iter += 1;
        
        # Inner optimization: find new optimal MO
        D3_sym = iszero(center - Wc.center)
        xy_orders, z_orders = select_orders(max_xy_order, max_z_order, D3_sym)
        n_SAGTOs = length(xy_orders) * length(z_orders)

        spreads_init = ones(n_SAGTOs) .* 1/2        
        optim_res = NaN

        # Optimize w.r. to spreads or center and spreads if non centered SAGTOs.
        if D3_sym
            tmp_kwargs = (; center, xy_orders, z_orders)
            optim_res = optimize(spreads -> f(spreads; tmp_kwargs...),
                                 (G, spreads) -> g!(G, spreads; D3_sym, tmp_kwargs...),
                                 spreads_init,
                                 optim_method,
                                 optim_options).minimizer
        else
            tmp_kwargs=(;xy_orders, z_orders)
            optim_res = optimize(center_and_spreads -> f(center_and_spreads[4:end];
                                                         center=center_and_spreads[1:3],
                                                         tmp_kwargs...),
                                 (G, center_and_spreads) -> g!(G, center_and_spreads; D3_sym, tmp_kwargs...),
                                 [π_bond_center..., spreads_init...],
                                 optim_method,
                                 optim_options).minimizer
        end
        tmp_kwargs = (;xy_orders, z_orders, in_linesearch=false)
        Φ_opti, _ = D3_sym ? f(optim_res; center, tmp_kwargs...) : f(optim_res[4:end]; center=optim_res[1:3], tmp_kwargs...)

        # Add new MO to the MO basis and project Wn
        push!(basis_functions, Φ_opti)
        compressed_wannier, coeffs = project_wannier_on_basis(Wc)

        # Compute new res and check CV
        residual = wannier .- compressed_wannier
        error = Hs_norm(basis_SC, residual; s) / Hs_norm(basis_SC, wannier; s)
        (error < tol) && (converged = true)

        # TODO: Add storage of the CompressedWannier
        # TODO: Add step that cures conditioning problems

        # Actualize compressed wannier structure and info
        Wc.basis_functions = basis_functions
        Wc.coefficients = coeffs
        Wc.residual = residual
        Wc.error = error

        # Change symmetry requirement to explore both D3_sym and non D3_sym part of Wn        
        (n_iter ≥ 1) && (center = D3_sym ? rand(3) : copy(Wc.center))

        # Actualize info
        info = merge(info, (;Wc, converged, n_iter))
        callback(info)
    end
    info 
end
