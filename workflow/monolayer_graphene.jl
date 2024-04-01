using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using WannierIO
using Wannier

### TODO: Setup parameters

function Graphene(; d=20u"Å", kwargs...)
    function basis()
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
    function scf(n_bands_converge=15)
        a = 2.641u"Å"  # Graphene Lattice constant
        lattice = [a  -a/2    0;
                   0  √3*a/2  0;
                   0     0    d]        
        C = ElementPsp(:C, psp=load_psp("hgh/pbe/c-q4"))
        atoms     = [C, C]
        positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
        model  = model_PBE(lattice, atoms, positions)
        basis  = PlaneWaveBasis(model; kwargs...)
        nbandsalg = AdaptiveBands(basis.model; n_bands_converge)
        self_consistent_field(basis; nbandsalg, tol=1e-5);
    end
    function wannierize(scfres; wannier_plot=false, fileprefix="wannier_output/graphene")
        # Initial guesses for the 5 valence bands of graphene
        C_Z = 6
        s_guess(center) = DFTK.HydrogenicWannierProjection(center, 2, 0, 0, C_Z)
        pz_guess(center) = DFTK.HydrogenicWannierProjection(center, 2, 1, 0, C_Z)
        positions = scfres.basis.model.positions
        projections = [
            # Note: fractional coordinates for the centers!
            # 3 bond-centered 2s hydrogenic orbitals to imitate σ bonds
            s_guess((positions[1] + positions[2]) / 2),
            s_guess((positions[1] + positions[2] + [0, -1, 0]) / 2),
            s_guess((positions[1] + positions[2] + [-1, -1, 0]) / 2),
            # 2 atom-centered 2pz hydrogenic orbitals
            pz_guess(positions[1]),
            pz_guess(positions[2]),
        ]
        wannier_model = Wannier.Model(scfres;
                                      fileprefix,
                                      n_bands=scfres.n_bands_converge,
                                      n_wannier=5,
                                      projections,
                                      dis_froz_max=ustrip(auconvert(u"eV", scfres.εF))+1
                                      ) # maximum frozen window, for example 1 eV above Fermi level
        U = disentangle(wannier_model, max_iter=200);
        wannier_model, U
    end
    function vesta_plot(wann_model, prefix)
        unkdir = splitpath(prefix)[1]
        Wannier.write_realspace_wf(prefix, wann_model; n_supercells=3, unkdir)       
        nothing
    end

    (;basis, scf, wannierize, vesta_plot)
end
