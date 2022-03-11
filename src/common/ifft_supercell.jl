"""
Convert a given Bloch decomposition into a 3 dimensional tensor adapted to 
ifft in the supercell.
"""
function Bloch_to_tensor(basis::PlaneWaveBasis, basis_supercell::PlaneWaveBasis,
                           ψ_fourier)
    ψ_fourier_supercell = zeros(ComplexF64, Tuple(basis_supercell.fft_size))
    for (ik, kpt) in enumerate(basis.kpoints)
        ikpG = DFTK.index_G_vectors.(basis_supercell,
                                     DFTK.Gplusk_vectors_in_supercell(basis, kpt))
        ψ_fourier_supercell[ikpG] .= ψ_fourier[ik]
    end
    ψ_fourier_supercell
end

function G_to_r_supercell(basis::PlaneWaveBasis, ψ_fourier)
    basis_unfold = DFTK.unfold_bz(basis)
    basis_supercell = DFTK.cell_to_supercell(basis_unfold, touch_atoms=false)
    ψ_fourier_supercell = Bloch_to_tensor(basis_unfold, basis_supercell, ψ_fourier)
    G_to_r(basis_supercell, ψ_fourier_supercell)
end
