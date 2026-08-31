#
# Hopping validation of the tight-binding chain (approach A, design note 09).
#
# Compares S(R) = ⟨w₀|w_R⟩ and T(R) = ⟨w₀|−½∇²|w_R⟩ from the true plane-wave
# Wannier (Fourier phase translation + grid) against the Gaussian-compressed
# counterpart (rigid translation + analytic integrals).
#
# R sampling (design note 09): intralayer a₁, a₁+a₂, 2a₁; pseudo-interlayer
# z-shift of 3.35 Å with lateral offset δ ∈ {0, (a₁+a₂)/3, a₁/2}. R = 0 is
# included as a consistency check against note 08.
#
# Usage:
#   julia --project=. workflow/validate_tb_hoppings.jl          # Ecut = 50
#   W2G_ECUT=15 julia --project=. workflow/validate_tb_hoppings.jl   # local smoke
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic
using LinearAlgebra
using Printf

const ECUT = parse(Int, get(ENV, "W2G_ECUT", "50"))
const D_SLAB = 10.0        # Å, interlayer vacuum used for the wannierization
const KGRID = [5, 5, 1]
const D_INTERLAYER = austrip(3.35u"Å")   # Bernal-like interlayer distance

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
include(joinpath(ROOT, "workflow", "monolayer_graphene.jl"))

compressed_file(Ecut) = joinpath(ROOT, "data",
    Ecut == 15 ? "compressed_wannier_55_H1.json" :
                 "compressed_wannier_55_H1_Ecut$(Ecut).json")
reference_file(Ecut) = joinpath(ROOT, "workflow", "wannier_functions",
                                "wannier_pz_Ecut-$(Ecut).json")

for f in (compressed_file(ECUT), reference_file(ECUT))
    isfile(f) || error("missing data file: $f")
end

@info "Building bases" ECUT D_SLAB KGRID
basis_uc = Graphene(; d=D_SLAB * u"Å", kgrid=KGRID, Ecut=ECUT).basis()
basis_sc = DFTK.cell_to_supercell(basis_uc)

# Unit-cell lattice vectors, Cartesian Bohr (DFTK stores them as columns).
a₁ = basis_uc.model.lattice[:, 1]
a₂ = basis_uc.model.lattice[:, 2]
ẑd = [0.0, 0.0, D_INTERLAYER]

hoppings = [
    ("R=0 (check)",        zeros(3)),
    ("intra a1",           a₁),
    ("intra a1+a2",        a₁ + a₂),
    ("intra 2a1",          2a₁),
    ("inter AA (δ=0)",     ẑd),
    ("inter AB-like",      ẑd + (a₁ + a₂) / 3),
    ("inter mid (δ=a1/2)", ẑd + a₁ / 2),
]

ref = W2G.read_wannier_function(reference_file(ECUT))
w_fourier = ref.wannier
n_G = length(DFTK.G_vectors(basis_sc, only(basis_sc.kpoints)))
@assert length(w_fourier) == n_G "Wannier/basis mismatch: len(w)=$(length(w_fourier)) vs n_G=$n_G — check (d, Ecut, kgrid)"
@assert norm(w_fourier) ≈ 1 "reference Wannier is not normalized (‖w‖ = $(norm(w_fourier)))"

Wc = W2G.CompressedWannier(compressed_file(ECUT))
@assert !isempty(Wc.basis_functions) "compressed Wannier carries no basis function"

println()
println("Hoppings S(R), T(R)   (Ecut = $ECUT, d = $D_SLAB Å, kgrid = $KGRID, norm-corrected)")
println("  n basis functions (GTO) = ", length(Wc.basis_functions))
println("  compression error (H¹)  = ", Wc.error)
println()
@printf("  %-20s %8s │ %12s %12s %9s │ %12s %12s %9s %9s\n",
        "R", "|R| (Å)", "S_ref", "S_gto", "relerr_S",
        "T_ref", "T_gto", "relerr_T", "err/T(0)")
for (label, R) in hoppings
    out = W2G.compare_hopping(Wc, w_fourier, basis_sc, R)
    @printf("  %-20s %8.3f │ %12.3e %12.3e %9.2e │ %12.3e %12.3e %9.2e %9.2e\n",
            label, ustrip(auconvert(u"Å", norm(R))),
            out.S_ref, out.S_gto, out.rel_err_S,
            out.T_ref, out.T_gto, out.rel_err_T, out.abs_err_T_per_T0)
end
println()
println("  T(0) reference (norm-corrected) = ",
        @sprintf("%.10f", W2G.compare_hopping(Wc, w_fourier, basis_sc, zeros(3)).T0_ref))
println("  On-site reference values (note 08): rel err T(0) = 0.41 % (Ecut 15), 1.89 % (Ecut 50)")
println("  H¹ compression errors           : 5.59 % (Ecut 15), 12.38 % (Ecut 50)")
println()
