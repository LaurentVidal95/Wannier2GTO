using Wannier2GTO
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using WannierIO
using Wannier

function Graphene(; kgrid=[5,5,1], Ecut=15, kshift=zeros(Float64, 3))
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

    function wannierize(scfres; prefix="wannier_output/graphene", plot_wannier=false)
        exclude_bands = DFTK._default_exclude_bands(scfres)
        basis, ψ, eigenvalues = DFTK.unfold_scfres_wannier(scfres, exclude_bands)
        μ, σ = 0.0, 0.01
        f = DFTK.scdm_f_erfc(basis, eigenvalues, μ, σ)
        n_wann = 5
        A = DFTK.compute_amn_scdm(basis, ψ, n_wann, f)
        wann_model, _ = run_wannier(
            scfres;
            fileprefix=prefix,
            n_wann,
            A,
            dis_froz_max=0.1,
        );
    end
    (;scf, wannierize)
end
