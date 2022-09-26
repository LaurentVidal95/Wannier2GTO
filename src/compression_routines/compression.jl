"""
Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
on an optimized symmetry adapted AO basis.
The basis is constructed by greedy iterations described in Ref[].
"""
function compress_graphene_pz_wannier(basis, Wn, α, D3_sym_angle;
                  basis_SC = DFTK.cell_to_supercell(basis),
                  s=0, # Choice of Hs scalar product
                  Res_init = Wn,  # Initial residual for restarts
                  save_file="compression_output.dat",
                  ## Outer loop parameters
                  max_xy_order=3, # Maximal order in (x,y) plane for the SAGTO polynomial
                  max_z_order=3,  # Maximal order in z axis for the SAGTO polynomial
                  ## Inner loop parameters
                  ζ_min=1e-3,     # Minimal possible spread to avoid NAN
                  max_iter=100,
                  thresh=1e-2, 
                                      optim_method=ConjugateGradient(linesearch=BackTracking(order=3)),
                  callback=default_callback()
                  )

    @assert  size(Wn,2)==1 "Wn is to be given as a single supercell vector"

    # Precomputation
    basis_SC = DFTK.cell_to_supercell(basis)
    Wn_norm = Hs_norm(basis_SC, Wn; s)

    # Populate info with outer and inner loop objects
    converged = false
    n_iter = zero(Int64)
    Wn_proj = zero(Wn)
    Res = Res_init
    
    Φ = zero(Wn)
    Φs = GaussianPolynomial{Real}[]     # Φs = Vector{Vector{ComplexF64}}[]
    proj_coeffs = ComplexF64[]

    # Initialize inner optimization loop parameters (for optimizing J)
    D3_sym = true
    pol_orders = select_orders(max_xy_order, max_z_order, D3_sym)

    info = (; s, converged, n_iter, Wn_proj, Res, Φ, Φs, proj_coeffs, D3_sym,
            pol_orders, basis_SC, ζ_min)
    callback(info)

    # Function for the inner optimzization loop.
    # return_coeffs=true returns J AND the coefficients obtained by solving the least square
    # problem in J.
    function f(X; D3_sym=true, in_linesearch=true)
        X[1:2] = D3_sym ? zeros(2) : [X[1], D3_sym_angle]
        J(info, α, X[1], X[2], X[3:end]; in_linesearch)
    end

    # while ( (info.n_iter < max_iter) && !info.converged )
    #     n_iter += 1;
    
    # Find optimal combination and spread.
    N_ζs = prod(length.(pol_orders)) # Number of spreads to optimize
    res_opti = optimize(X -> f(X; D3_sym), ones(N_ζs+2), optim_method,
             Optim.Options(g_abstol=1e-6, show_trace=true), autodiff=:forward)
    info, J_opti, Λ_opti, α_opti = f(res_opti.minimizer, in_linesearch=false)
    ζs_opti = res_opti.minimizer[3:end]
    
    # Actualise info
    Φs = push!(info.Φs, info.Φ)
    info = merge(info, (;Φs=Φs))
    info
        # Store MOs information
        # push!(Φs_coeffs, [Λ_opti, α_opti, ζs_opti])
        # !(isempty(save_file)) && (writedlm(save_file, Φs_coeffs))

    #     # Add new MO to the MO basis and project Wn
    #     push!(Φs, Φ_current)
    #     w_proj, proj_coeffs = Hs_projection_on_AO_basis(basis_SC, Wn, Φs, s=s)
    #     # Check conditioning issues
    #     isnothing(w_proj) && (pop!(Φs); return w_proj, Φs, proj_coeffs)
        
    #     # Actualize Res and check convergence.
    #     Res = Wn .- w_proj
    #     Res_Hs_norm = Hs_norm(basis_SC, Res, s=s) / Wn_norm
    #     (Res_Hs_norm < thresh) && (converged = true)
        
    #     # Print current info
    #     infos = [n_iter, J_opti, Res_Hs_norm]
    #     @printf "%-6i %-10.6f %-10.6f \n" infos...
        #     flush(stdout)
        # (n_iter > 2) && (D3_sym = D3_sym ? false : true)
    # end
    
    # (converged) && (@info "Converged"); (!converged) && (@warn "Not converged")
    # w_proj, Φs, proj_coeffs, Φs_coeffs
end
