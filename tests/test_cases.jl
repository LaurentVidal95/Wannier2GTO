using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using WannierIO
using Wannier
using wannier90_jll

const meV_to_hartree = ustrip(u"Eh_au", 1u"meV")
const hartree_to_meV = ustrip(u"meV", 1u"Eh_au")

default_lattice_length() = 4.66  # in Bohr
default_z_axis_length() = 30

struct GrapheneSystem{T<:Real}
    # Geometry
    lattice::Matrix{T}
    atoms::AbstractVector{DFTK.ElementPsp}
    positions::AbstractVector{AbstractVector{T}}
    # Tight Binding
    R_vectors::AbstractVector{AbstractVector{T}}
    tol::T
end

# Generic function that works for both graphene and stacked bilayer graphene.
geometry(G::GrapheneSystem) = (G.lattice, G.atoms, G.positions)
R_vectors_cart(G::GrapheneSystem) = map(R->G.lattice*R, G.R_vectors)

"""
    GrapheneSystem(atoms, positions; lattice_length, z_axis_length)

Creates a monolayer or stacked bilayer graphene system through respectively
[`MonolayerGraphene`](@ref) or [`StackedBilayerGraphene`](@ref).
"""
function GrapheneSystem(atoms, positions;
                        lattice_length=default_lattice_length(),
                        tol = 1e-8)
    @assert length(positions) == length(atoms)

    # Assemble Graphene system
    lz = z_axis_length / lattice_length
    lattice = lattice_length .* [  1/2     1/2  0;
                                 -√3/2    √3/2  0;
                                     0     0   lz]
    R_vectors = compute_R_vectors(tol)
    GrapheneSystem(lattice, atoms, position, tol,)
end

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
    function scf(n_bands_converge=15)
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
        
        nbandsalg = AdaptiveBands(basis.model; n_bands_converge)
        self_consistent_field(basis; nbandsalg, tol=1e-5);
    end
    function wannierize(scfres; wannier_plot=false, prefix="wannier_output/graphene")
        exclude_bands = DFTK._default_exclude_bands(scfres)
        basis, ψ, eigenvalues = DFTK.unfold_scfres_wannier(scfres, exclude_bands)

        n_wann = 8
        centers = Vector{Vector{Float64}}()
        [push!(centers, basis.model.positions[1]) for _ in 1:2]
        [push!(centers, basis.model.positions[2]) for _ in 1:2]
        l = [0, 1, 0, 1]
        guess = DFTK.guess_amn_hydrogenic(basis, centers, l)

        A = DFTK.compute_amn(basis, ψ, guess)
        (wannier_plot) && (WannierIO.write_amn("$prefix.amn", A))
        wann_model = only(run_wannier(
            scfres;
            fileprefix=prefix,
            n_wann,
            A,
            dis_froz_max=0.1,
            wannier_plot,
        ));

        wann_model
    end
    function vesta_plot(wann_model, prefix)
        unkdir = splitpath(prefix)[1]
        Wannier.write_realspace_wf(prefix, wann_model; n_supercells=3, unkdir)        
        nothing
    end
    function wannierize_w90(scfres; prefix="w90_output/graphene")
        run_wannier90(scfres;
                      fileprefix=prefix,
                      n_wann=5,  # Number of MLWFs
                      # following are inputs to wannier90
                      num_print_cycles=25,
                      num_iter=200,
                      ##
                      dis_win_max=19.0,
                      dis_froz_max=0.1,
                      dis_num_iter=300,
                      dis_mix_ratio=1.0,
                      ##
                      bands_plot=true,
                      ##
                      wannier_plot=true,
                      wannier_plot_format="cube",
                      wannier_plot_supercell=5,
                      write_u_matrices=true,
                      write_xyz=true,
                      translate_home_cell=true,
                      );
    end

    (;basis, scf, wannierize, vesta_plot, wannierize_w90)
end

# TODO: const of function ? Chose a side
const default_interlayer_dist = 6.45 #bohr TODO type

@doc raw"""
Define a system with two layers of graphene A and B, aligned (no twist) with the top one (B)
shifted in the longitudinal direction.
For a given interlayer distance ``d ≥ 0`` and a disregistry ``y ∈ ℝ²``, the layer A is
shifted by ``[0, 0, -d/2]`` with respect to the origin and the B layer is shifted by
``[y₁, y₂, d/2]``.
 """
function StackedBilayerGraphene(d_red; #interlayer distance in reduced coordinates
                                disregistry = zeros(eltype(d),2),
                                kwargs...)
    function model()
        G = Graphene(;)
        model = G.basis().model
        z_axis_length = model.lattice[end]
        d_red = d/z_axis_length
        shift = [disregistry..., d_red]
        
        # Shifted monolayer positions at z=-d/2
        lattice = model.lattice
        positions = model.positions .- Ref([0., 0., shift[3]/2])
        atoms = model.atoms
        # Add second layer at z=+d/2
        append!(atoms, atoms)
        append!(positions, positions .+ Ref(shift))
        
        DFTK.Model(lattice, atoms, positions)
    end
end
