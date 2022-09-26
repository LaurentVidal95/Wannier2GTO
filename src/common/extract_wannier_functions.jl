using JSON3

# Extract Wannier90 output
function read_w90_output_mat(file, n_bands_tot)
    Uks = readdlm(file); num_kpts, num_ψout, num_ψin = Uks[2,1:3];
    len_Uk = num_ψin*num_ψout

    list_Uk = []
    for k in 1:num_kpts
        i_Uk = 4 + (k-1)*(len_Uk+1)
        Uk = reshape([x+im*y for (x,y) in eachrow(Uks[i_Uk:i_Uk+(len_Uk-1), 1:2])],
                     num_ψin, num_ψout)
        (norm(Uk'Uk - I) > 1e-8) && (error("U matrix has non-orthogonal columns"))
        # Complete U matrices by zero to ease the writing of `apply_U_matrices`
        (n_bands_tot ≠ num_ψin) &&
            ( Uk = vcat(Uk, zeros(ComplexF64, n_bands_tot - num_ψin, num_ψout)) )
        push!(list_Uk, Uk)
    end
    list_Uk
end

@doc raw"""
For a given set of bands ``\{Ψ_{nk}\}_{1\leq n \leq N_in}``
and a Wannier90 file containing unitary matrices (for disentanglement or
wannierization) returns a set of band ``{Ψ_out_{nk}}_{1\leq n \leq N_out}``
with:
``Ψ_out_{nk} = \sum\limits_{m=1}^{N_in} U_{mn}Ψ_{nk}``
``N_in > N_out`` for disentanglement.
``N_in = N_out`` for wannierization.
"""
function apply_U_matrices(file, ψ)
    # read .mat file
    n_bands_tot = size(ψ[1], 2)
    num_kpts = length(ψ)
    Uks = read_w90_output_mat(file, n_bands_tot)
    [ψ[k]*Uks[k] for k in 1:num_kpts]
end

@doc raw"""
Apply wannierization U matrices to return Wannier functions as a table ``W[k][G,n]``
The prefix in argument is the one returned by the `run_wannier90` routine above.
As in the wannierization process, the scfres has to be unfolded.
"""
function extract_wannier_functions(prefix::String, scfres::NamedTuple)
    scfres_unfold = DFTK.unfold_bz(scfres); ψ_unfold = scfres_unfold.ψ
    # Apply disentanglement if needed
    isfile("$(prefix)_u_dis.mat") && (ψ_unfold =
                       apply_U_matrices("$(prefix)_u_dis.mat", ψ_unfold))
    # Construct Wannier functions
    apply_U_matrices("$(prefix)_u.mat", ψ_unfold)
end
