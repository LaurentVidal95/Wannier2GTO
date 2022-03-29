using Pkg
Pkg.activate("/home/lvidal/Documents/CERMICS/these/wannier_and_GTOs/graphene/GrapheneTB")
using GrapheneTB
using Plots
using DFTK

################ Element with custom potential
struct ElementGaussian <: DFTK.Element
    symbol
    α  # Prefactor
    L  # Width of the Gaussian nucleus
    Z  # Needed so that pymatgen doesn't raise errors
end

function DFTK.local_potential_real(el::ElementGaussian, r::Real)
    -el.α / (√(2π) * el.L) * exp(- (r / el.L)^2 / 2)
end

function DFTK.local_potential_fourier(el::ElementGaussian, q::Real)
    ## = ∫ V(r) exp(-ix⋅q) dx
    -el.α * exp(- (q * el.L)^2 / 2)
end

DFTK.charge_nuclear(el::ElementGaussian) = el.Z
nucleus = ElementGaussian(:X, 1.0, 0.5, 1)

################ SCF and Wannierization
function scf_1d_chain(b; n_atoms=2, n_bands=2,
                   kgrid=[10,1,1], Ecut=500, fft_grid=[550,1,1],
                   kshift=[0,0,0])
    # lattice
    lat = 10 .* [[1 0 0]; [0 b 0]; [0 0 b]];

    # Atoms
    atoms = [nucleus => [[0.5, 0, 0]]]
    (n_atoms==2) && (atoms = [nucleus => [[0.2, 0, 0], [0.8, 0, 0]]])
    
    # Model and basis
    terms = [Kinetic(), AtomicLocal()]
    model = Model(lat, atoms=atoms, n_electrons=n_atoms,
                  terms=terms, spin_polarization=:spinless);
    basis = PlaneWaveBasis(model, Ecut=Ecut, kgrid=kgrid, fft_size=fft_grid, kshift=kshift)
    self_consistent_field(basis, tol=1e-8, n_bands=n_bands, n_ep_extra=0);
end

using wannier90_jll
function run_wannie90_1d_chain(scfres; prefix="w90_output/1d_chain")
    atom_chain_wannierization_kwargs = (n_bands=2, n_wannier=2,
                                      write_u_matrices=true,
                                      skip_B1_tests=true, shell_list=1,
                                      write_xyz=true)
    DFTK.run_wannier90(scfres; atom_chain_wannierization_kwargs...)
end

################### Tight Binding for 2D atom chain
"""
Apply τ_R toa given function, periodic of the supercell.
"""
function τR_1D(basis_unfold, R, ψ::Vector{Matrix{T}}) where {T<:Complex}
    τ_ψ = []
    for (ik, k) in enumerate(basis_unfold.kpoints)
        prefac_k = cis.(-dot(G[1], R) for G in Gplusk_vectors_cart(basis_unfold, k))
        push!(τ_ψ, prefac_k .* ψ[ik])
    end
    τ_ψ
end

"""
Warning, w are to be normalized over the supercell before hand. Otherwise, divide result by √(num_k)
"""
function tight_binding_ham(basis_unfold, ρ_unfold, w; R_neigh=[10,20,30,40])
    ham = Hamiltonian(basis_unfold, ρ=ρ_unfold)
    # Extract both wannier at R=0
    extract_wn(w) = ([wk[:,1] for wk in w], [wk[:,2] for wk in w])
    Wn_R_list = [extract_wn(w)...]
    for R in R_neigh
        w1_R, w2_R = extract_wn(τR_1D(basis_unfold, R, w))
        push!(Wn_R_list, w1_R, w2_R)
    end
    N_neigh = length(R_neigh)
    [w1_R'*(ham*w2_R) for w1_R in Wn_R_list, w2_R in Wn_R_list]
end

tight_binding_ham(scfres_unfold, R_neigh, w)= H_RR(scfres_unfold.basis, scfres_unfold.ρ, R_neigh, w)

################## Finding back H_k
"""
2×2 matrix s.t. [H_R]_μν = ⟨W_μ0 | H | W_νR⟩
"""
function H_R(ham_unfold::Hamiltonian, w, R)
    w_R = τR_1D(ham_unfold.basis, R, w)
    w_R'*(ham*w)
end

"""
Recover H(k) by Bloch transform of all H_R
"""
function Hk_tight_binding(k, ham_unfold, w, R1, R2; R_neigh=[10,20,30,40])
    HR_list = H_R.(Ref(ham_unfold), Ref(w), R_neigh)
    dists(R) = .- [dot(k, R + Ra - Rb) for Ra in [R1,R2], Rb in  [R1,R2]]
    sum(HR .* cis.(dists(R)) for (HR, R) in zip(HR_list, R_neigh))
end
