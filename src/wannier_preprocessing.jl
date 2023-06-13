# This file contains the routines to prepare the output of the wannierization
# done with Wannier.jl to be given to the compression algorithm.

using Wannier

@doc raw"""
Apply the unitary matrices contained in wann_model to obtain optimized
Bloch states from initial ones.
"""
function apply_U_matrices(wann_model, ψ)
    Uks = wann_model.U
    [ψ[k][:,1:wann_model.n_bands]*Uks[:,:,k] for k in 1:wann_model.n_kpts]
end
@doc raw"""
Goes from a list of vectors `Wn[kpt][iG]` in the unit cell to a single vector at 
Γ point of the supercell.
TODO: Keeping the structure as vect of kpt might be faster. Try benchmarking.
"""
function convert_wannier_to_supercell(Wn, basis, basis_supercell)
    Wn_supercell = cell_to_supercell(Wn, DFTK.unfold_bz(basis), basis_supercell)
    sum(eachcol(Wn_supercell))
end

"""
Small hack to provide the angle between axis x and the first axis of the
(x,y)-plane D3 symmetry of pz-like wannier functions.
"""
function find_π_bond_axis(basis_supercell::PlaneWaveBasis, wannier, wannier_center;
                          debug_sign=1)

    @info "Computing the π_bond axis in polar coordinates"
    s_orb(α, ζ) = GaussianPolynomial([(0,0,1)], [1.], α, ζ)
    α(λ, θ) = polar_to_cartesian_coords(wannier_center, λ, θ)

    # Try to fit a s-like orbital to the part of the wannier function located
    # on the π-bond axis.
    res = optimize(X->norm(wannier + debug_sign*s_orb(α(X[1],X[2]), X[3])(basis_supercell)),
                   # Guess roughly close to wanted axis by experience.
                   [1/3., -1/2, 2] .+ wannier_center,
                   ConjugateGradient(linesearch=BackTracking(order=3)))
    r, θ = res.minimizer[1:2]

    # If r is approximately zero the hack failed or the s orbitals have wrong signs.
    # Launch optimization again by inverting the sign of the s orbital.
    if norm(r) < 1e-2
        (debug_sign==-1) && (error("r≈0, try to manualy chose the other pz wannier orbital"))
        @warn "Failed optimization, trying with opposite sign"
        return find_π_bond_axis(basis_supercell, wannier, wannier_center; debug_sign=-1)
    end

    (;r, θ)
end

@doc raw"""
Extract from scf computation and wannierization all the data needed for compression.
To be fed directly to the compression routine.
"""
function prepare_for_compression(wann_model, scfres; wannier_manual_selection=nothing)
    # Compute supercell basis
    basis = scfres.basis
    basis_supercell = cell_to_supercell(basis)

    # Compute Wannier functions from Bloch states and unitary matrices in wann_model
    # Unfold eventual k-point symmetries
    scfres_UN = DFTK.unfold_bz(scfres)
    Wns = apply_U_matrices(wann_model, scfres_UN.ψ)

    # Select a pz wannier (they have wider spread that π-bonds)
    wann_res = Wannier.omega(wann_model)
    i_pz = isnothing(wannier_manual_selection) ? findmax(wann_res.ω)[2] : wannier_manual_selection
    center = wann_res.r[:,i_pz] # TODO: * ang_to_bohr
    @info "Selected Wannier function: n°$(i_pz)"

    # Convert given Wannier from unit cell to supercell convention
    wannier = convert_wannier_to_supercell([Wns_k[:,i_pz] for Wns_k in Wns], basis, basis_supercell)
    normalize!(wannier) # Renormalize just in case

    # Identify π-bond axis in polar coordinates
    π_bond = find_π_bond_axis(basis_supercell, wannier, center)

    (; basis_supercell, wannier, center, π_bond)
end

#
# Old code adapted to the use of Wannier90
#

# using DelimitedFiles
#
# Extract Wannier90 output
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

# function apply_U_matrices_w90(file, ψ)
#     # read .mat file
#     n_bands_tot = size(ψ[1], 2)
#     num_kpts = length(ψ)
#     Uks = read_w90_output_mat(file, n_bands_tot)
#     [ψ[k]*Uks[k] for k in 1:num_kpts]
# end

# @doc raw"""
# Apply wannierization U matrices to return Wannier functions as a table ``W[k][G,n]``
# The prefix in argument is the one returned by the `run_wannier90` routine above.
# As in the wannierization process, the scfres has to be unfolded.
# """
# function extract_wannier_functions_w90(prefix::String, scfres::NamedTuple)
#     scfres_unfold = DFTK.unfold_bz(scfres); ψ_unfold = scfres_unfold.ψ
#     # Apply disentanglement if needed
#     isfile("$(prefix)_u_dis.mat") && (ψ_unfold =
#                        apply_U_matrices_w90("$(prefix)_u_dis.mat", ψ_unfold))
#     # Construct Wannier functions
#     apply_U_matrices_w90("$(prefix)_u.mat", ψ_unfold)
# end
