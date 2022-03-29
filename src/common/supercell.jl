function cell_to_supercell_vector(basis_unfold, basis_SC, W, n_band)
    num_kpG = length(G_vectors(basis_SC, only(basis_SC.kpoints)))
    W_vec = zeros(ComplexF64, num_kpG)
    cell_supercell_mapping(kpt) = DFTK.index_G_vectors.(basis_SC, Ref(basis_SC.kpoints[1]),
                                  DFTK.Gplusk_vectors_in_supercell(basis_unfold, kpt))
    for (ik, kpt) in enumerate(basis_unfold.kpoints)
        W_vec[cell_supercell_mapping(kpt)] .= W[ik][:,n_band]
    end
    W_vec ./ sqrt(length(basis_unfold.kpoints))
end
