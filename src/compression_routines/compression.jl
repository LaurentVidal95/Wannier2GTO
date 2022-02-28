using DFTK
using Optim; using LineSearches
using Printf

"""
Provides the angle between axis x and the first axis of the (x,y)-plane
D3 symmetry of pz-like wannier functions.
"""
function find_D3_sym_axis(basis::PlaneWaveBasis, W_pz, α0)
    α(λ, θ) = polar_to_cartesian_coords(α0, λ, θ)
    res = optimize(X->norm(W_pz .- Bloch_decomposition(basis, s_orb(α(X[1],X[2]), X[3]))),
                   [1., -1/2, 1/2], # Guess roughly close to wanted axis by experience.
                   ConjugateGradient(linesearch=BackTracking(order=3)),
                   Optim.Options(show_trace=true))
    res.minimizer[2]
end

"""
    Compress pz-like Wannier function of graphene (given as a Fourier coefficient table)
    on an optimized symmetry adapted AO basis.
    The basis is constructed by greedy iterations described in Ref[].
"""
function compress_graphene_pz_wannier(basis::PlaneWaveBasis, wn, α;
                                      s=0,            # Choice of Hs norm. Default is L2.
                                      ζ_min=1e-3,     # Minimal possible spread to avoid NAN
                                      max_xy_order=3, # Maximal order in (x,y) plane for the SAGTO polynomial
                                      max_z_order=3,  # Maximal order in z axis for the SAGTO polynomial
                                      Res_init = wn,  # initial residual
                                      save_file="SAGTOs_coeffs_H$(s).dat", # File in which to store the SAGTOs coeffs
                                      # J! optim parameters
                                      max_iter=100, thresh=1e-2,
                                      optim_method=ConjugateGradient(
                                          linesearch=BackTracking(order=3)),
                                      θ=zero(Float64),
                                      )
    # Initialise needed objects
    converged = false; iter = zero(Int64);
    w_proj = [zero(wn_k) for wn_k in wn]; Res = Res_init .- w_proj
    wn_norm = Hs_norm(basis, wn, s=s)

    # SAGTOs storage
    proj_coeffs = []
    Φs = []; Φs_coeffs = [];
    
    # Initialize global parameters going in optimization loop
    opti_init_point = [];
    Φ_current = deepcopy(wn)
    D3_sym=true
    xy_orders = zero(Float64); z_orders = zero(Float64)

    optimize_ζ(ζs; return_all=false) =
        J!(basis, α, ζs, zero(Float64), θ, Res, Φ_current,
                        s=s, ζ_min=ζ_min,
                        xy_orders=xy_orders, z_orders=z_orders,
                        return_all=return_all)
    optimize_ζr(X; return_all=false) =
        J!(basis, α, X[1:end-1], X[end], θ, Res, Φ_current,
                             s=s, ζ_min=ζ_min,
                             xy_orders=xy_orders, z_orders=z_orders,
                             return_all=return_all)
    f_opti = optimize_ζ

    # Compute blobs angle
    (iszero(θ)) && ( @show θ = find_blobs_angle(basis, wn, α) )
    
    # Print header
    header = ["Iter", "J", "||Residual||"]
    println("-"^40); @printf "%-10s %-10s %10s \n" header...; println("-"^40)
    flush(stdout)
    
    while( (iter < max_iter) && !converged )
        iter += 1;

        # Todo : better strategy than alternating blobs and not blobs ?
        (rem(iter-1,2)==0) && (D3_sym=true)
        (rem(iter-1,2)==1) && (D3_sym=false)

        # Choose the type of polynomial used for this iteration
        xy_orders, z_orders = compute_orders(max_xy_order, max_z_order, D3_sym=D3_sym)
        N_ζs = prod(length.([xy_orders, z_orders]))

        # Choose between optimization of the central part and the "blobs".
        (D3_sym) && (f_opti = optimize_ζ; opti_init_point = ones(N_ζs);)
        !(D3_sym) && (f_opti = optimize_ζr; opti_init_point = vcat(ones(N_ζs), 1.);)

        # Find optimal combination and spread. Resulting SAGTO is stored in Φ_current
        Φ_current = deepcopy(wn)
        res_opti = optimize(X -> f_opti(X), opti_init_point, optim_method,
                           Optim.Options(g_abstol=1e-6, show_trace=true))
        # Call J! again to fix issues with Φ_current
        J_opti, Λ_opti, α_opti = f_opti(res_opti.minimizer, return_all=true)

        # Store MOs information
        ζs_opti = res_opti.minimizer; !(D3_sym) && (_ = pop!(ζs_opti))
        push!(Φs_coeffs, [Λ_opti, α_opti, ζs_opti])
        !(isempty(save_file)) && (writedlm(save_file, Φs_coeffs))

        # Add new MOs to the MO basis and project  
        push!(Φs, Φ_current)
        w_proj, proj_coeffs = Hs_projection_on_AO_basis(basis, wn, Φs, s=s)
        # Check conditioning issues
        isnothing(w_proj) && (pop!(Φs); return w_proj, Φs, proj_coeffs)
        
        # Actualize Res and check convergence.
        Res = wn .- w_proj
        Res_Hs_norm = Hs_norm(basis, Res, s=s) / wn_norm
        (Res_Hs_norm < thresh) && (converged = true)
        
        # Print current info
        infos = [iter, J_opti, Res_Hs_norm]
        @printf "%-6i %-10.6f %-10.6f \n" infos...
        flush(stdout)
    end
    
    (converged) && (@info "Converged"); (!converged) && (@warn "Not converged")
    w_proj, Φs, proj_coeffs, Φs_coeffs
end
