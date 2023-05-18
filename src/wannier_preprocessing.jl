# This file contains the routines to prepare the output of the wannierization
# done with Wannier.jl to be given to the compression algorithm.

using Wannier

@doc raw"""
For a given set of bands ``\{Ψ_{nk}\}_{1\leq n \leq N_in}``
and a Wannier90 file containing unitary matrices (for disentanglement or
wannierization) returns a set of band ``{Ψ_out_{nk}}_{1\leq n \leq N_out}``
with:
``Ψ_out_{nk} = \sum\limits_{m=1}^{N_in} U_{mn}Ψ_{nk}``
``N_in > N_out`` for disentanglement.
``N_in = N_out`` for wannierization.
"""
function apply_U_matrices(wann_model, ψ)
    Uks = wann_model.U
    [ψ[k][:,1:wann_model.n_bands]*Uks[:,:,k] for k in 1:wann_model.n_kpts]
end

function convert_wannier_to_supercell(Wn, basis, basis_SC)
    !(norm(Wn) ≈ 1) && @warn "The given wannier is not normalized"
    Wn_SC = cell_to_supercell(Wn, DFTK.unfold_bz(basis), basis_SC)
    sum(eachcol(Wn_SC))
end

@doc raw"""
Construct one s orbital of center ``α`` and spread ``ζ``.
"""
s_orb(α, ζ) = GaussianPolynomial([(0,0,1)], [1.], α, ζ)

"""
TODO: better routine not based on shady optimization ?
Provides the angle between axis x and the first axis of the (x,y)-plane
D3 symmetry of pz-like wannier functions.
"""
function find_π_bond_axis(basis_SC::PlaneWaveBasis, W_pz, wannier_center)
    @info "Computing the π_bond axis in polar coordinates"
    α(λ, θ) = polar_to_cartesian_coords(wannier_center, λ, θ)
    res = optimize(X->norm(W_pz + s_orb(α(X[1],X[2]), X[3])(basis_SC)),
                   [1., -1/2, 1/2], # Guess roughly close to wanted axis by experience.
                   ConjugateGradient(linesearch=BackTracking(order=3)))
    r, θ = res.minimizer[1:2]
    (norm(r) < 1e-2) && (error("r ≈ 0.. Try changing s sign in optimization"))
    r, θ
end

@doc raw"""
Extract from scf computation and wannierization all the data needed for compression.
To be fed directly to the compression routine.
"""
function prepare_for_compression(wann_model, scfres)
    # Compute supercell basis
    basis = scfres.basis
    basis_SC = cell_to_supercell(basis)

    # Compute Wannier functions from Bloch states and unitary matrices in wann_model
    # Unfold eventual k-point symmetries
    scfres_UN = DFTK.unfold_bz(scfres)
    Wns = apply_U_matrices(wann_model, scfres_UN.ψ)

    # Select a pz wannier (they have wider spread that π-bonds)
    wann_res = Wannier.omega(wann_model)
    i_pz = findmax(wann_res.ω)[2]
    center = wann_res.r[:,i_pz]
    Wn_pz = convert_wannier_to_supercell([Wns_k[:,i_pz] for Wns_k in Wns], basis, basis_SC)
    
    # Identify π-bond axis in polar coordinates
    r, θ = find_π_bond_axis(basis_SC, Wn_pz, center)
    (; basis_supercell=basis_SC, wannier=Wn_pz, center, π_axis=(r,θ))
end

#
# Old code adapted to the use of Wannier90
#

# using JSON3

# # Extract Wannier90 output
# function read_w90_output_mat(file, n_bands_tot)
#     Uks = readdlm(file); num_kpts, num_ψout, num_ψin = Uks[2,1:3];
#     len_Uk = num_ψin*num_ψout

#     list_Uk = []
#     for k in 1:num_kpts
#         i_Uk = 4 + (k-1)*(len_Uk+1)
#         Uk = reshape([x+im*y for (x,y) in eachrow(Uks[i_Uk:i_Uk+(len_Uk-1), 1:2])],
#                      num_ψin, num_ψout)
#         (norm(Uk'Uk - I) > 1e-8) && (error("U matrix has non-orthogonal columns"))
#         # Complete U matrices by zero to ease the writing of `apply_U_matrices`
#         (n_bands_tot ≠ num_ψin) &&
#             ( Uk = vcat(Uk, zeros(ComplexF64, n_bands_tot - num_ψin, num_ψout)) )
#         push!(list_Uk, Uk)
#     end
#     list_Uk
# end

# @doc raw"""
# Apply wannierization U matrices to return Wannier functions as a table ``W[k][G,n]``
# The prefix in argument is the one returned by the `run_wannier90` routine above.
# As in the wannierization process, the scfres has to be unfolded.
# """
# function extract_wannier_functions(prefix::String, scfres::NamedTuple)
#     scfres_unfold = DFTK.unfold_bz(scfres); ψ_unfold = scfres_unfold.ψ
#     # Apply disentanglement if needed
#     isfile("$(prefix)_u_dis.mat") && (ψ_unfold =
#                        apply_U_matrices("$(prefix)_u_dis.mat", ψ_unfold))
#     # Construct Wannier functions
#     apply_U_matrices("$(prefix)_u.mat", ψ_unfold)
# end

