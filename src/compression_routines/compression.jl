"""
Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
on an optimized symmetry adapted AO basis.
The basis is constructed by greedy iterations described in Ref[].
"""
function compress_graphene_pz_wannier(basis_SC::PlaneWaveBasis{TR},
                  Wn::AbstractArray{TC},
                  wannier_center::AbstractVector{TR},
                  π_bond_axis;
                  ### Kwargs ...
                  s=0,            # Choice of Hs scalar product
                  res_init = Wn,  # Initial residual for restarts
                  prefix="compressed_wannier",
                  callback=default_callback(),
                  ## Outer loop parameters
                  max_xy_order=3, # Maximal order in (x,y) plane for the SAGTO polynomial
                  max_z_order=3,  # Maximal order in z axis for the SAGTO polynomial
                  ## Inner loop parameters
                  ζ_min=1e-3,     # Minimal possible spread to avoid NAN
                  max_iter=20,    # Maximal number of compression iterations
                  tol=1e-2,       # maximum attaigned for now is ≈ 5%
                  optim_method=ConjugateGradient(),
                  optim_options=Optim.Options(g_abstol=1e-5, show_trace=true),
                  ) where {TR <: Real, TC <: Complex}

    @assert  size(Wn,2)==1 "The wannier function is to be given as a single supercell vector"

    # Compute initial values and populate info
    J = NaN
    res = res_init
    res_norm = NaN
    Wn_norm = Hs_norm(basis_SC, Wn; s)
    Wn_proj = []
    Φ = empty(Wn) # Last computed MO
    Φs = Vector{ComplexF64}[] # All MOs
    D3_sym = true
    pol_orders = select_orders(max_xy_order, max_z_order, D3_sym)
    converged = false
    n_iter = zero(Int64)

    info = (; basis_SC, J, res, res_norm, converged, n_iter,  # info for general optimization
            Φ, Φs, Wn_norm, pol_orders, D3_sym,               # info for inner optimization
            s, ζ_min)                                         # general parameters
    callback(info)

    output = Dict{Int64, Any}()

    # Function for the inner optimzization loop.
    # in_linesearch=false returns J AND the coefficients obtained by solving the least square
    # problem in J.
    function f(X; D3_sym=true, in_linesearch=true)
        X[1:2] = D3_sym ? zeros(2) : vcat(X[1], θ)
        inner_optimization(info, wannier_center, X[1], X[2], X[3:end]; in_linesearch)
    end

    while ( (info.n_iter < max_iter) && !info.converged )
        n_iter += 1;

        Φs = info.Φs
        
        # Inner optimization: find new optimal MO
        N_AOs = prod(length.(pol_orders)) # Number of AOs
        res = optimize(X -> f(X; D3_sym), vcat([1.0, π_bond_axis], ones(N_AOs)),
                       optim_method, optim_options, autodiff=:forward)
        info, J, Λ, α = f(res.minimizer; D3_sym, in_linesearch=false)
        ζs = res.minimizer[3:end]

        # Add new MO to the MO basis and project Wn
        # Compute new res and check CV
        push!(Φs, info.Φ)
        Wn_proj = Hs_projection_on_AO_basis(basis_SC, Wn, Φs; s)
        res = Wn .- Wn_proj.func
        res_norm = Hs_norm(basis_SC, res; s) / Wn_norm
        (res_norm < tol) && (converged = true)

        # Save compression data for current iter
        output[zero(Int64)] = Wn_proj.coeffs
        output[n_iter] = [Λ, α, ζs, pol_orders]
        open(io->JSON3.write(io, output, allow_inf=true), prefix*".json", "w")

        # Change symmetry requirement to explore both D3_sym and non D3_sym part of Wn
        (n_iter ≥ 1) && (D3_sym = D3_sym ? false : true)
        pol_orders = select_orders(max_xy_order, max_z_order, D3_sym)

        # Actualize info
        info = merge(info, (;Φs=Φs,  res=res, res_norm=res_norm, D3_sym=D3_sym,
                            n_iter=n_iter, J=J, converged=converged, pol_orders=pol_orders))
        callback(info)
    end
    (;res=info.res_norm, projection=Wn_proj, n_iter, converged)
end
