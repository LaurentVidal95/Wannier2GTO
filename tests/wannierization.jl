ℛ_red = [0 -1 0; 1 1 0; 0 0 1] # Rotation of π/3 in reduced coordinates
prefac = cis.( -map(G->-dot(basis.model.lattice * basis.model.positions[2],G),
                   G_vectors_cart(basis_supercell,basis_supercell.kpoints[1]))
               )

