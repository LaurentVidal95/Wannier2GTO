using DFTK
using Plots
using DelimitedFiles

include("../include_src.jl")

# 2D Graphene in a 3D box
d = 10
a = 2.641*1.88973  # Graphene Lattice constant
lattice = [a  -a/2    0;
           0  √3*a/2  0;
           0     0    d]

C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
atoms  = [C => [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]]
model  = model_PBE(lattice, atoms)
basis  = PlaneWaveBasis(model; Ecut=15, kgrid=[5, 5, 1])
scfres = self_consistent_field(basis; n_bands=15, tol=1e-8);

# Wannierization
using wannier90_jll  # Needed to make run_wannier90 available
prefix = run_wannier90(scfres;
              fileprefix="w90_output_no_blobs/graphene",
              n_wannier=5,
              num_print_cycles=25,
              num_iter=200,
              ##
              dis_win_max=19.0,
              dis_froz_max=0.1,
              dis_num_iter=300,
              dis_mix_ratio=1.0,
              ##
              wannier_plot=true,
              wannier_plot_format="cube",
              wannier_plot_supercell=5,
              write_xyz=true,
              translate_home_cell=true,
              ##
              write_u_matrices=true,         
              );


# Extract wannier Fourier coefficients
scfres_unfold = DFTK.unfold_bz(scfres);
basis_unfold = scfres_unfold.basis;
w = extract_wannier(prefix, scfres_unfold);

# Extract Wannier centers
ang_to_bohr = 1.88973;
αs = (readdlm(prefix*"_centres.xyz")[3:3+6, 2:4]) .* ang_to_bohr
