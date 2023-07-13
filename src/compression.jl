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
                  max_iter=15,    # Maximal number of compression iterations
                  tol=1e-2,       # maximum attaigned for now is ≈ 5%
                  optim_method=ConjugateGradient(),
                  optim_options=Optim.Options(g_abstol=1e-5, show_trace=false),
                  file="compressed_wannier.json")

    # Extract compressed wannier data
    wannier = Wc.wannier
    basis_supercell = Wc.basis_supercell
    π_bond_center = polar_to_cartesian_coords(Wc.center, π_bond.r, π_bond.θ)
    residual = Wc.residual
    error = Wc.error
    s = Wc.error_norm
    residual_norm = Hˢ_norm(basis_supercell, residual; s)
    basis_functions = Wc.basis_functions
    coefficients = Wc.coefficients

    # Check that the initial wannier is normalized and in superclel convention
    @assert  size(wannier,2)==1 "The wannier function is to be"*
        " given as a single supercell vector"
    @assert norm(wannier) ≈ 1 "The Wannier function must be normalized before compression"
    # Check that the plane-wave basis is big enough for good num dot products
    _check_L²_dot_precision(basis_supercell)

    # Center of the next basis function to be built
    center = copy(Wc.center)

    # Populate info and print header
    converged = false
    n_iter = zero(Int64)
    info = (; Wc, converged, n_iter)
    callback(info)

    # Objective function to optimize at each iteration. The function produce an optimal
    # basis function (which best approximate the current residual) as linear combination of
    # SAGTOs with given spreads and center, respecting the symmetries encoded in xy_orders and z_orders.
    # When "in_linesearch", the function only returns the Hˢ error of approximation of the
    # residual by the optimal basis function.
    function f(spreads::Vector{T1}; center::Vector{T2}, xy_orders, z_orders,
               in_linesearch=true) where {T1, T2<:Real}
        # Hack to avoid conditioning problems with small spreads
        findmin(spreads)[1] < ζ_min && return Inf

        # Compute optimal basis function for given spreads, center and residual (stored in Wc)
        SAGTOs = SAGTO_basis(center, spreads, xy_orders, z_orders)
        Φ, error = optimal_basis_function(Wc, SAGTOs)

        (in_linesearch) && (return error)
        Φ, error
    end

    while ( (info.n_iter < max_iter) && !info.converged )
        n_iter += 1;

        # Inner optimization: find new optimal basis function
        # Check if the current basis function has D3_sym or only parity w.r.t. graphene plane
        D3_sym = iszero(center - Wc.center)

        # Select polynomial orders that respect the current basis function's symmetry.
        xy_orders, z_orders = select_orders(max_xy_order, max_z_order, D3_sym)
        n_SAGTOs = length(xy_orders) * length(z_orders)

        # Optimize w.r. to spreads or center and spreads if non centered SAGTOs.
        spreads_init = ones(n_SAGTOs) .* 1/2
        optim_res = NaN
        if D3_sym
            tmp_kwargs = (; center, xy_orders, z_orders)
            optim_res = optimize(spreads -> f(spreads; tmp_kwargs...),
                                 spreads_init,
                                 optim_method,
                                 optim_options,
                                 autodiff=:forward).minimizer # make it :reverse
        else
            tmp_kwargs=(;xy_orders, z_orders)
            optim_res = optimize(center_and_spreads -> f(center_and_spreads[4:end];
                                                         center=center_and_spreads[1:3],
                                                         tmp_kwargs...),
                                 [π_bond_center..., spreads_init...],
                                 optim_method,
                                 optim_options,
                                 autodiff=:forward).minimizer # make it :reverse
        end
        tmp_kwargs = (;xy_orders, z_orders, in_linesearch=false)
        Φ_opti, _ = D3_sym ? f(optim_res; center, tmp_kwargs...) : f(optim_res[4:end]; center=optim_res[1:3], tmp_kwargs...)

        # Add new MO to the MO basis and project Wn
        # Normalize before projection to cure conditioning
        Φ_opti = normalize(basis_supercell, Φ_opti)
        push!(basis_functions, Φ_opti)
        compressed_wannier, coeffs = project_wannier_on_basis(Wc)

        # Compute new res and check CV
        residual = wannier .- compressed_wannier
        error = Hs_norm(basis_supercell, residual; s) / Hs_norm(basis_supercell, wannier; s)
        (error < tol) && (converged = true)

        # TODO: Add step that cures conditioning problems

        # Actualize compressed wannier structure and info
        Wc.basis_functions = basis_functions
        Wc.coefficients = coeffs
        Wc.residual = residual
        Wc.error = error
        store(Wc; file)

        # Hack: change symmetry requirement to explore both D3_sym and non D3_sym part of the wannier
        (n_iter ≥ 1) && (center = D3_sym ? rand(3) : copy(Wc.center))

        # Actualize info
        info = merge(info, (;Wc, converged, n_iter))
        callback(info)
    end
    # Store compressed wannier with true coefficients
    # See "fft_supercell" in "GaussianPolynomials.jl" for further details.
    Wc = fix_coefficients!(basis_supercell, info.Wc)
    store(Wc; file)
    info = merge(info, (;Wc))
end

# Prototype for ReverseDiff gradient. For now produced memory errors.

function objective_f_and_g!(Wc, ζ_min)
    # Objective function to optimize at each iteration. The function produce an optimal
    # basis function (which best approximate the current residual) as linear combination of
    # SAGTOs with given spreads and center, respecting the symmetries encoded in xy_orders and z_orders.
    # When "in_linesearch", the function only returns the Hˢ error of approximation of the
    # residual by the optimal basis function.
    function f(spreads::Vector{T}; center::Vector{T}, xy_orders, z_orders,
               in_linesearch=true) where {T<:Real}
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
    # g_tape = ReverseDiff.GradientTape(f, etc...) see ReverseDiff doc
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
