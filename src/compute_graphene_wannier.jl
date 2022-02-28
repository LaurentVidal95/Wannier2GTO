using DFTK
using wannier90_jll
using DelimitedFiles

function scf_graphene(n_bands; kgrid=[8,8,1], Ecut=15)
    # Lattice
    ang_to_bohr = 1.88973; a = ang_to_bohr*2.641; # from wiki
    a_1 = [a; 0; 0];
    rot_120_deg = [-1/2   -√3/2   0;
                   √3/2    -1/2   0;
                   0         0    0]
    a_2 = rot_120_deg*a_1; a_3 = [0; 0; 20]
    lattice = hcat(a_1, a_2, a_3)
    
    C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
    atoms = [C => [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]]
    
    model = model_PBE(lattice, atoms)
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)
    self_consistent_field(basis, n_bands=n_bands);
end

run_wannier90_graphene(scfres; prefix="w90_output/graphene") =
    run_wannier90(scfres, 
                  fileprefix=prefix,
                  n_bands=13, n_wannier=5,
                  write_u_matrices=".TRUE.",
                  write_xyz=".TRUE.",
                  # plot
                  # wannier_plot=true, wannier_plot_supercell=4,
                  # disentanglement
                  dis_win_max=19.0, dis_froz_max=0.19, dis_num_iter=300)
