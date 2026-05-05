#
# Phase A — diagnose the cond(S) blow-up identified in thesis 5.4.2.3.
#
# Loads a precomputed pz wannier (Ecut=15, smallest), reconstructs the matching
# DFTK plane-wave basis (no SCF needed), then runs the greedy compression with
# a custom callback that logs cond(S) of the accumulated basis at every
# iteration. Output: stdout stream + CSV file.
#

using Wannier2GTO
import Wannier2GTO as W2G
using LinearAlgebra
using DFTK
using DFTK.Unitful
using Printf
using Dates

const PROJECT_DIR = joinpath(splitpath(pathof(Wannier2GTO))[1:(end - 2)]...)
const WORKFLOW_DIR = joinpath(PROJECT_DIR, "workflow")
const OUTPUT_DIR = joinpath(WORKFLOW_DIR, "diag_outputs")
mkpath(OUTPUT_DIR)
const TIMESTAMP = Dates.format(now(), "yyyymmdd-HHMMSS")
const LOG_FILE = joinpath(OUTPUT_DIR, "diag_$(TIMESTAMP).csv")

include(joinpath(WORKFLOW_DIR, "monolayer_graphene.jl"))

# --- Setup ----------------------------------------------------------------
println("=" ^ 70)
println("Phase A diagnostic — cond(S) blow-up reproduction")
println("=" ^ 70)

# Match the params used to produce the JSON (found by `scan_basis_params.jl`):
# d=10 Å (vacuum half the default), Ecut=15, kgrid=[5,5,1].
low_params = (; d=10u"Å", kgrid=[5, 5, 1], Ecut=15)
println("\nBuilding PlaneWaveBasis with $low_params ...")
G = Graphene(; low_params...)
basis = G.basis()
basis_supercell = DFTK.cell_to_supercell(basis)
println("Plane-wave G-vectors at Γ: ",
        length(DFTK.G_vectors(basis_supercell, basis_supercell.kpoints[1])))

println("\nLoading precomputed wannier (Ecut=15) ...")
data = W2G.read_wannier_function(joinpath(WORKFLOW_DIR,
                                          "wannier_functions/wannier_pz_Ecut-15.json"))

println("Wannier center = ", data.center)
println("π_bond = (r=$(data.π_bond.r), θ=$(data.π_bond.θ))")
println("‖wannier‖ = ", round(norm(data.wannier); digits=10))

Wc = CompressedWannier(basis_supercell, normalize(data.wannier), data.center)
Wc.error_norm = 1   # H^1 norm

# --- Instrumented callback ------------------------------------------------
diag_log = NamedTuple[]

function diag_callback(info)
    n_iter = info.n_iter
    Wc = info.Wc
    Φs = Wc.basis_functions
    n_funcs = length(Φs)

    if n_funcs ≥ 1
        Φs_Four = [Φ(Wc.basis_supercell) for Φ in Φs]
        S = W2G.Hˢ_overlap(Wc.basis_supercell, Φs_Four; s=Wc.error_norm)
        σ = svdvals(Hermitian(Matrix(S)))
        σ_min = σ[end]
        σ_max = σ[1]
        cond_S = σ_max / σ_min
        n_SAGTOs = sum(length(Φ.SAGTOs) for Φ in Φs)

        # Direction 3: log the parameters of the function just added.
        Φ_new = Φs[end]
        ζ_new = Φ_new.SAGTOs[1].spread
        α_new = Φ_new.SAGTOs[1].center
        coeff_norm = norm(Φ_new.coeffs)
        @printf("    [new Φ] ζ=%.4e  α=[%+.3f, %+.3f, %+.3f]  ‖coeffs‖=%.3e\n",
                ζ_new, α_new[1], α_new[2], α_new[3], coeff_norm)
    else
        σ_min = NaN
        σ_max = NaN
        cond_S = NaN
        n_SAGTOs = 0
    end

    @printf("[iter=%2d]  n_funcs=%2d  n_SAGTOs=%3d  err=%.3e  cond(S)=%.3e  σ_min=%.3e  σ_max=%.3e\n",
            n_iter, n_funcs, n_SAGTOs, Wc.error, cond_S, σ_min, σ_max)
    flush(stdout)

    push!(diag_log,
          (; iter=n_iter, n_funcs=n_funcs, n_SAGTOs=n_SAGTOs,
           error=Wc.error, cond_S=cond_S, sigma_min=σ_min, sigma_max=σ_max))
    nothing
end

# --- Run compression with instrumentation ---------------------------------
# Tune `regularization` here. Available: :none, :tikhonov, :svd_truncation,
# :pivoted_cholesky. `regularization_param` is the relative scale (ε for
# :tikhonov, λ_min/λ_max threshold for the truncation methods).
const REGULARIZATION = :tikhonov
const REGULARIZATION_PARAM = 1e-8

println("\nLaunching compression with regularization=$(REGULARIZATION) param=$(REGULARIZATION_PARAM)")
println("Header: iter | n_funcs | n_SAGTOs | err | cond(S) | σ_min | σ_max\n")

res = try
    compress_graphene_pz_wannier(Wc, data.π_bond;
                                 max_iter=15,
                                 tol=1e-3,         # below thesis result, push beyond
                                 max_xy_order=3,
                                 max_z_order=3,
                                 regularization=REGULARIZATION,
                                 regularization_param=REGULARIZATION_PARAM,
                                 callback=diag_callback,
                                 file=joinpath(OUTPUT_DIR, "diag_compressed_wannier.json"))
catch e
    @warn "compress_graphene_pz_wannier raised: $e"
    nothing
end

# --- Persist log ----------------------------------------------------------
println("\nWriting CSV log to $LOG_FILE")
open(LOG_FILE, "w") do io
    println(io, "iter,n_funcs,n_SAGTOs,error,cond_S,sigma_min,sigma_max")
    for r in diag_log
        @printf(io, "%d,%d,%d,%.10e,%.10e,%.10e,%.10e\n",
                r.iter, r.n_funcs, r.n_SAGTOs, r.error, r.cond_S,
                r.sigma_min, r.sigma_max)
    end
end

println("\nDone. Iterations logged: ", length(diag_log))

# --- Direction 3: post-run inspection of pairwise H¹ overlaps -------------
println("\n" * "=" ^ 70)
println("Pairwise H¹ overlap |⟨Φ_i, Φ_j⟩| (all Φ are H¹-normalized via S\\Γ)")
println("=" ^ 70)

if length(Wc.basis_functions) ≥ 2
    Φs_Four = [Φ(basis_supercell) for Φ in Wc.basis_functions]
    n = length(Φs_Four)
    M = zeros(n, n)
    for i in 1:n, j in 1:n
        M[i, j] = abs(W2G.Hˢ_dot(basis_supercell, Φs_Four[i], Φs_Four[j];
                                  s=Wc.error_norm))
    end

    # Compact table with sorted spreads
    println("\nBasis function summary (sorted by spread ζ):")
    @printf("  %3s  %12s  %s\n", "i", "ζ", "α (primary center)")
    spreads = [Φ.SAGTOs[1].spread for Φ in Wc.basis_functions]
    centers = [Φ.SAGTOs[1].center for Φ in Wc.basis_functions]
    perm = sortperm(spreads)
    for k in perm
        @printf("  %3d  %12.4e  [%+.4f, %+.4f, %+.4f]\n",
                k, spreads[k], centers[k][1], centers[k][2], centers[k][3])
    end

    # Top redundant pairs
    println("\nTop 10 most-overlapping pairs (excluding diagonal):")
    pairs = [(i, j, M[i, j]) for i in 1:n for j in (i + 1):n]
    sort!(pairs; by=p -> -p[3])
    for k in 1:min(10, length(pairs))
        i, j, v = pairs[k]
        @printf("  Φ_%-2d ↔ Φ_%-2d : |⟨·,·⟩|_H¹ = %.4f   ζ_i=%.3e  ζ_j=%.3e\n",
                i, j, v, spreads[i], spreads[j])
    end

    # Save full overlap matrix
    overlap_file = joinpath(OUTPUT_DIR, "overlap_$(TIMESTAMP).csv")
    open(overlap_file, "w") do io
        for i in 1:n
            println(io, join((@sprintf("%.6e", M[i, j]) for j in 1:n), ","))
        end
    end
    println("\nFull overlap matrix saved to $overlap_file")
else
    println("(fewer than 2 functions; nothing to compare)")
end
