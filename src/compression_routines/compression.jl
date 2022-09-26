"""
Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
on an optimized symmetry adapted AO basis.
The basis is constructed by greedy iterations described in Ref[].
"""
function compress_graphene_pz_wannier(basis, Wn, α; s=0, # Choice of Hs norm.
               Res_init = Wn,  # initial residual
               ## Outer loop parameters
               max_xy_order=3, # Maximal order in (x,y) plane for the SAGTO polynomial
               max_z_order=3,  # Maximal order in z axis for the SAGTO polynomial
               D3_sym_axis=zero(Float64), # polar coord of one axis of D3 (x,y)-plane symmetry
               save_file="compression_output.dat",
               ## Inner loop parameters
               ζ_min=1e-3,     # Minimal possible spread to avoid NAN
               max_iter=100, thresh=1e-2, 
               optim_method=ConjugateGradient(linesearch=BackTracking(order=3)),
               )
    @assert  size(Wn,2)==1 "Wn is to be given as a single supercell vector"

    # Initialise outer loop objects and parameters
    converged = false
    iter = zero(Int64); Wn_proj = zero(Wn); Res = copy(Wn);
    basis_SC = DFTK.cell_to_supercell(basis, touch_atoms=false)
    Wn_norm = Hs_norm(basis_SC, Wn, s=s)

    # SAGTOs storage
    proj_coeffs = []; Φs = []; Φs_coeffs = [];
    Φ_current = copy(Wn)

    # Initialize inner optimization loop parameters (for optimizing J)
    opti_init_point = [];
    pol_orders = [];
    D3_sym = true

    optimize_D3_sym(ζs; return_all=false) =
        J!(basis_SC, Res, α, 0., 0., ζs, Φ_current,
           s=s, ζ_min=ζ_min, pol_orders=pol_orders, return_all=return_all)
    optimize_no_D3_sym(X; return_all=false) =
        J!(basis_SC, Res, α, X[end], D3_sym_axis, X[1:end-1], Φ_current,
           s=s, ζ_min=ζ_min, pol_orders=pol_orders, return_all=return_all)
    f_opti = optimize_D3_sym

    # Compute blobs angle
    iszero(D3_sym_axis) && (@show D3_sym_axis = find_D3_sym_axis(basis_SC, Wn, α))

    # Print header
    header = ["Iter", "J", "||Residual||"]
    println("-"^40); @printf "%-10s %-10s %10s \n" header...; println("-"^40)
    flush(stdout)
    
    while( (iter < max_iter) && !converged )
        iter += 1;

        # Todo : better strategy than alternating blobs and not blobs ?
        (rem(iter-1,2)==0) && (D3_sym=true)
        (rem(iter-1,2)==1) && (D3_sym=false)

        # Choose polynomial part orders and function to optimze w.r. to D3_sym
        pol_orders = select_orders(max_xy_order, max_z_order, D3_sym)
        N_ζs = prod(length.(pol_orders))
        (D3_sym) && (f_opti = optimize_D3_sym; opti_init_point = ones(N_ζs);)
        !(D3_sym) && (f_opti = optimize_no_D3_sym;
                      opti_init_point = vcat(ones(N_ζs), 1.);)

        # Find optimal combination and spread. Resulting MO is stored in Φ_current
        Φ_current = copy(Wn)
        res_opti = optimize(X -> f_opti(X), opti_init_point, optim_method,
                            Optim.Options(g_abstol=1e-6, show_trace=true))
        # Call J! again to fix issues with Φ_current and recover coefficients
        J_opti, Λ_opti, α_opti = f_opti(res_opti.minimizer, return_all=true)

        # Store MOs information
        ζs_opti = res_opti.minimizer; !(D3_sym) && (_ = pop!(ζs_opti))
        push!(Φs_coeffs, [Λ_opti, α_opti, ζs_opti])
        !(isempty(save_file)) && (writedlm(save_file, Φs_coeffs))

        # Add new MO to the MO basis and project Wn
        push!(Φs, Φ_current)
        w_proj, proj_coeffs = Hs_projection_on_AO_basis(basis_SC, Wn, Φs, s=s)
        # Check conditioning issues
        isnothing(w_proj) && (pop!(Φs); return w_proj, Φs, proj_coeffs)
        
        # Actualize Res and check convergence.
        Res = Wn .- w_proj
        Res_Hs_norm = Hs_norm(basis_SC, Res, s=s) / Wn_norm
        (Res_Hs_norm < thresh) && (converged = true)
        
        # Print current info
        infos = [iter, J_opti, Res_Hs_norm]
        @printf "%-6i %-10.6f %-10.6f \n" infos...
        flush(stdout)
    end
    
    (converged) && (@info "Converged"); (!converged) && (@warn "Not converged")
    w_proj, Φs, proj_coeffs, Φs_coeffs
end
