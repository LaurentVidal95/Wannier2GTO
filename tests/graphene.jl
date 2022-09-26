"""
Scf and wannierization routine for Graphene
"""

using DFTK
using wannier90_jll

struct System{F}
    scf::F
    wannierize::F
end

function Graphene(; n_bands=10, kgrid=[8,8,1], Ecut=15, kshift=zeros(Float64, 3))
    function scf(; n_bands=n_bands, kgrid=kgrid, Ecut=Ecut, kshift=kshift, kwargs...)
        # Lattice
        ang_to_bohr = 1.88973; a = ang_to_bohr*2.641; # from wiki
        a_1 = [a; 0; 0];
        rot_120_deg = [-1/2   -√3/2   0;
                       √3/2    -1/2   0;
                       0         0    0]
        a_2 = rot_120_deg*a_1; a_3 = [0; 0; 20]
        lattice = hcat(a_1, a_2, a_3)
        
        C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
        atoms = [C,C]
        positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
        
        model = model_PBE(lattice, atoms, positions)
        basis = PlaneWaveBasis(model; Ecut, kgrid, kshift)
        self_consistent_field(basis, n_bands, kwargs...);
    end

    function wannierize(scfres; prefix="w90_output/graphene", plot_wannier=false)
        kwargs = (fileprefix=prefix,
                  n_bands=13, n_wannier=5,
                  # Data to recover Wannier as tensor in julia
                  length_unit="Bohr"
                  write_u_matrices=".TRUE.",
                  write_xyz=".TRUE.",
                  # Disentanglement
                  dis_win_max=19.0, dis_froz_max=0.19, dis_num_iter=300)
        # Plot if needed
        (plot_wannier) && (kwargs = merge(kwargs, (; wannier_plot=true, wannier_plot_supercell=4)))
        run_wannier90(scfres; kwargs...)
    end
    System(scf, wannierize)
end

