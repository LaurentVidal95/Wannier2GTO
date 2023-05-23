using Wannier2GTO
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using WannierIO
using Wannier

"""
Routines that handle the SCF and wannierization of Graphene with fixed
parameters.
"""

function Graphene(; kgrid=[5,5,1], Ecut=15, kshift=zeros(Float64, 3))
    function basis()
        d = 10u"Å"
        a = 2.641u"Å"  # Graphene Lattice constant
        lattice = [a  -a/2    0;
                   0  √3*a/2  0;
                   0     0    d]
        
        C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
        atoms     = [C, C]
        positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
        model  = model_PBE(lattice, atoms, positions)
        basis  = PlaneWaveBasis(model; Ecut, kgrid)
    end
    function scf()
        d = 10u"Å"
        a = 2.641u"Å"  # Graphene Lattice constant
        lattice = [a  -a/2    0;
                   0  √3*a/2  0;
                   0     0    d]
        
        C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
        atoms     = [C, C]
        positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
        model  = model_PBE(lattice, atoms, positions)
        basis  = PlaneWaveBasis(model; Ecut, kgrid)
        
        nbandsalg = AdaptiveBands(basis.model; n_bands_converge=15)
        self_consistent_field(basis; nbandsalg, tol=1e-5);
    end
    function wannierize(scfres; wannier_plot=false, prefix="wannier_output/graphene")
        exclude_bands = DFTK._default_exclude_bands(scfres)
        basis, ψ, eigenvalues = DFTK.unfold_scfres_wannier(scfres, exclude_bands)
        
        n_wann = 5
        centers = Vector{Vector{Float64}}()
        [push!(centers, basis.model.positions[1]) for _ in 1:2]
        [push!(centers, basis.model.positions[2]) for _ in 1:2]
        l = [0, 1, 0, 1]
        guess = DFTK.guess_amn_hydrogenic(basis, centers, l)
        
        A = DFTK.compute_amn(basis, ψ, guess)

        wann_model = only(run_wannier(
            scfres;
            fileprefix=prefix,
            n_wann,
            A,
            dis_froz_max=0.1,
            wannier_plot,
        ));

        (wannier_plot) && (WannierIO.write_amn("$(prefix).amn", A))
        
        # Remove tmp files
        !(wannier_plot) && rm(splitpath(prefix)[1], recursive=true)
        wann_model
    end
    function vesta_plot(wann_model, prefix)
        unkdir = splitpath(prefix)[1]
        Wannier.write_realspace_wf(prefix, wann_model; n_supercells=3, unkdir)        
        nothing
    end
    (;basis, scf, wannierize, vesta_plot)
end
