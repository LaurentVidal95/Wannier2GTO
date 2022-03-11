using WriteVTK

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
function extract_wannier_functions(prefix, scfres::NamedTuple)
    scfres_unfold = DFTK.unfold_bz(scfres); ψ_unfold = scfres_unfold.ψ
    # Apply disentanglement if needed
    isfile("$(prefix)_u_dis.mat") && (ψ_unfold =
                       apply_U_matrices("$(prefix)_u_dis.mat", ψ_unfold))
    # Construct Wannier functions
    apply_U_matrices("$(prefix)_u.mat", ψ_unfold)
end

"""
   Produce vts files that are readable with Paraview to plot Wannier functions.
   The wannier functions are given in Bloch decomposition format.
"""
function Bloch_to_vtk(basis::PlaneWaveBasis, w_fourier, prefix, n_band)
    basis_supercell = DFTK.cell_to_supercell(basis, touch_atoms=false)
    # Extract real supercell grid
    r_vec = collect(r_vectors_cart(basis_supercell))
    x = [r[1] for r in r_vec]; y = [r[2] for r in r_vec]; z = [r[3] for r in r_vec];
    # Perfom ifft
    wn_real = G_to_r_supercell(basis, [wk[:, n_band] for wk in w_fourier])
    vtk_grid(prefix*"_wannier_$(n_band)", x, y, z) do vtk
        vtk["value"]=real.(wn_real)
    end
    nothing
end
plot_wannier_function(scfres::NamedTuple, prefix::String, n_band) =
    Bloch_to_vtk(DFTK.unfold_bz(scfres.basis),
                 extract_wannier_functions(prefix, scfres), prefix, n_band)

# Extract stored data
function read_wannier_dir(dir)
    readdlm_complex(file) = readdlm(file, Float64) |> x -> Complex.(x[:,1], x[:,2])
    wn = [readdlm_complex(joinpath(dir,"wn_$(k).dat")) for k in 1:64];
    αn = vec(readdlm(joinpath(dir,"wn_center.dat")))
    θn = only(readdlm(joinpath(dir,"wn_blob_angle.dat")))
    wn, αn, θn
end
