#
# On-site one-body kinetic validation of the tight-binding chain.
#
# Compares ⟨w|−½∇²|w⟩ computed from the true plane-wave Wannier (grid) against
# the same integral computed from its Gaussian-compressed counterpart
# (analytic, via the native `kinetic_julia`).
#
# See notes/06-design-validation-TB-onsite.md. The V_KS term is deliberately
# out of scope here: it needs an scfres consistent with the supercell basis and
# is the subject of the follow-up step (L-2b).
#
# Usage:
#   julia --project=. workflow/validate_tb_onsite.jl          # Ecut = 50
#   W2G_ECUT=15 julia --project=. workflow/validate_tb_onsite.jl   # local smoke
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using LinearAlgebra
using Printf

const ECUT = parse(Int, get(ENV, "W2G_ECUT", "50"))
const D_SLAB = 10.0      # Å, interlayer vacuum used for the wannierization
const KGRID = [5, 5, 1]

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
include(joinpath(ROOT, "workflow", "monolayer_graphene.jl"))

# Compressed (Gaussian) Wannier: the Ecut=15 run kept the historical file name.
compressed_file(Ecut) = joinpath(ROOT, "data",
    Ecut == 15 ? "compressed_wannier_55_H1.json" :
                 "compressed_wannier_55_H1_Ecut$(Ecut).json")
reference_file(Ecut) = joinpath(ROOT, "workflow", "wannier_functions",
                                "wannier_pz_Ecut-$(Ecut).json")

for f in (compressed_file(ECUT), reference_file(ECUT))
    isfile(f) || error("missing data file: $f")
end

@info "Building supercell basis" ECUT D_SLAB KGRID
basis_sc = let
    basis = Graphene(; d=D_SLAB * u"Å", kgrid=KGRID, Ecut=ECUT).basis()
    DFTK.cell_to_supercell(basis)
end

# Reference: true plane-wave Wannier (supercell Γ point, Fourier coefficients).
ref = W2G.read_wannier_function(reference_file(ECUT))
w_fourier = ref.wannier

# Contract: the reference Wannier must live on exactly this basis, otherwise the
# grid integral below is meaningless.
n_G = length(DFTK.G_vectors(basis_sc, only(basis_sc.kpoints)))
@assert length(w_fourier) == n_G "Wannier/basis mismatch: len(w)=$(length(w_fourier)) vs n_G=$n_G — check (d, Ecut, kgrid)"
@assert norm(w_fourier) ≈ 1 "reference Wannier is not normalized (‖w‖ = $(norm(w_fourier)))"
@assert all(isfinite, w_fourier) "non-finite coefficients in the reference Wannier"

# Gaussian-compressed Wannier (phase-A greedy output, H¹ metric).
Wc = W2G.CompressedWannier(compressed_file(ECUT))
@assert !isempty(Wc.basis_functions) "compressed Wannier carries no basis function"

out = W2G.compare_onsite_kinetic(Wc, w_fourier, basis_sc)

println()
println("On-site kinetic ⟨w|-½∇²|w⟩   (Ecut = $ECUT, d = $D_SLAB Å, kgrid = $KGRID)")
println("  n basis functions (GTO)   = ", length(Wc.basis_functions))
println("  compression error (H¹)    = ", Wc.error)
println("  ‖w‖  reference / Gaussian = ", @sprintf("%.6f / %.6f", out.norm_ref, out.norm_gto))
println()
println("  raw:")
println("    T_ref  (grid, Wannier)  = ", @sprintf("%.10f", out.T_ref))
println("    T_gto  (analytic, GTO)  = ", @sprintf("%.10f", out.T_gto))
println("    relative error          = ", @sprintf("%.3e", out.rel_err))
println()
println("  norm-corrected (physical on-site TB element):")
println("    T_ref / ‖w‖²            = ", @sprintf("%.10f", out.T_ref_per_norm))
println("    T_gto / ‖w‖²            = ", @sprintf("%.10f", out.T_gto_per_norm))
println("    relative error          = ", @sprintf("%.3e", out.rel_err_per_norm))
println()
