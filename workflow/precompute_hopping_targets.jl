#
# Precompute hopping-loss targets (design note 12 §4): reference S(R), T(R)
# from the true plane-wave Wannier, for the training/validation R sets.
# Stored once in data/, so the optimization loop never touches plane waves
# for the hopping terms.
#
# Usage:
#   julia --project=. workflow/precompute_hopping_targets.jl          # Ecut = 50
#   W2G_ECUT=15 julia --project=. workflow/precompute_hopping_targets.jl
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using LinearAlgebra
using Printf

const ECUT = parse(Int, get(ENV, "W2G_ECUT", "50"))
const D_SLAB = 10.0      # Å, interlayer vacuum used for the wannierization
const KGRID = [5, 5, 1]
const SEED = 1234        # validation random-shell draws (reproducibility)
const N_RANDOM = 5

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
include(joinpath(ROOT, "workflow", "monolayer_graphene.jl"))

reference_file(Ecut) = joinpath(ROOT, "workflow", "wannier_functions",
                                "wannier_pz_Ecut-$(Ecut).json")
target_file(Ecut) = joinpath(ROOT, "data", "hopping_targets_Ecut$(Ecut).json")
isfile(reference_file(ECUT)) || error("missing data file: $(reference_file(ECUT))")

@info "Building bases" ECUT D_SLAB KGRID
basis_uc = Graphene(; d=D_SLAB * u"Å", kgrid=KGRID, Ecut=ECUT).basis()
basis_sc = DFTK.cell_to_supercell(basis_uc)
a₁ = basis_uc.model.lattice[:, 1]
a₂ = basis_uc.model.lattice[:, 2]

ref = W2G.read_wannier_function(reference_file(ECUT))
w_fourier = ref.wannier
n_G = length(DFTK.G_vectors(basis_sc, only(basis_sc.kpoints)))
@assert length(w_fourier) == n_G "Wannier/basis mismatch: len(w)=$(length(w_fourier)) vs n_G=$n_G"
@assert norm(w_fourier) ≈ 1 "reference Wannier is not normalized"

rsets = W2G.hopping_R_sets(a₁, a₂;
    d_inter=austrip(3.35u"Å"), z_val=(austrip(3.0u"Å"), austrip(3.7u"Å")),
    seed=SEED, n_random=N_RANDOM,
    r_min=austrip(2.0u"Å"), r_max=austrip(4.0u"Å"))
targets = W2G.build_hopping_targets(w_fourier, basis_sc, rsets)
W2G.store(targets; file=target_file(ECUT))

println()
println("Hopping targets (Ecut = $ECUT, seed = $SEED)  →  $(target_file(ECUT))")
println("  T(0) = ", @sprintf("%.10f", targets.T0_ref))
@printf("  %-18s %-11s %8s │ %12s %12s\n", "label", "set", "|R| (Å)", "S_ref", "T_ref")
for i in eachindex(targets.labels)
    @printf("  %-18s %-11s %8.3f │ %12.4e %12.4e\n",
            targets.labels[i], targets.sets[i],
            norm(targets.Rs[i]) / austrip(1u"Å"),
            targets.S_ref[i], targets.T_ref[i])
end
