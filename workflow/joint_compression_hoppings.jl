#
# Joint compression with the hopping-targeted loss (design note 12).
#
# Runs a paired μ-scan: the μ = 0 twin (pure H¹) and each μ > 0 share the same
# master seed, hence identical per-restart inits — the controlled comparison
# required by design §5. After each run, rebuilds the best basis and logs the
# per-criterion evaluation (independent PASS/FAIL) on the validation targets.
#
# Usage:
#   W2G_ECUT=15 julia --project=. workflow/joint_compression_hoppings.jl   # smoke
#   julia --project=. workflow/joint_compression_hoppings.jl              # Ecut = 50
#
using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful
using Random
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "monolayer_graphene.jl"))

# --- Configuration --------------------------------------------------------
const ECUT        = parse(Int, get(ENV, "W2G_ECUT", "50"))
const KGRID       = [5, 5, 1]
const D           = 10.0u"Å"
const N_CENTERED  = 5
const N_PIBOND    = 5
const K_RESTARTS  = parse(Int, get(ENV, "W2G_K", "3"))   # W2G_K=1: single seed (cheap scan)
const MASTER_SEED = 42          # shared across the whole scan (paired runs)
const MAX_ITER    = parse(Int, get(ENV, "W2G_MAX_ITER", "100"))
const ζ_MIN       = 1e-2
const ζ_MAX       = ECUT / 4.0
const ε_TIKHONOV  = 1e-8
# μ scan; ν = μ (design §2). Override with W2G_MUS (comma-separated) to split
# the scan into one SLURM job per μ — e.g. W2G_MUS=10 for a single value. The
# paired summary (criterion 4) then needs the μ=0 job's H¹ from its own log.
const MUS         = [parse(Float64, m) for m in
                     split(get(ENV, "W2G_MUS", "0,1,10,100,1000"), ",")]
const H1_REGRESSION_FACTOR = 1.5                 # criterion 4 threshold

const ROOT = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]...)
wannier_file = joinpath(ROOT, "workflow", "wannier_functions",
                        "wannier_pz_Ecut-$(ECUT).json")
targets_file = joinpath(ROOT, "data", "hopping_targets_Ecut$(ECUT).json")
for f in (wannier_file, targets_file)
    isfile(f) || error("missing data file: $f — run precompute_hopping_targets.jl first")
end

# --- Setup ----------------------------------------------------------------
println("Building plane-wave basis (Ecut=$ECUT)...")
basis = Graphene(; d=D, kgrid=KGRID, Ecut=ECUT).basis()
basis_supercell = DFTK.cell_to_supercell(basis)

data = W2G.read_wannier_function(wannier_file)
w_z_fourier = normalize(data.wannier)
π_bond_unit = let θ = data.π_bond.θ; [cos(θ), sin(θ), 0.0] end
targets = W2G.HoppingTargets(targets_file)
layout = W2G.JointLayout(N_centered=N_CENTERED, N_pibond=N_PIBOND)
kw = (log_ζ_min=log(ζ_MIN), log_ζ_max=log(ζ_MAX), π_bond_unit=π_bond_unit)

rebuild(flat) = begin
    Φs = W2G.params_to_basis_functions(flat, layout; kw...)
    c = W2G.joint_inner_solve(Φs, w_z_fourier, basis_supercell; ε=ε_TIKHONOV)
    (Φs, c)
end

# --- Paired scan ----------------------------------------------------------
results = Dict{Float64, Any}()
for μ in MUS
    ν = μ
    println("\n", "=" ^ 70)
    println("RUN μ = $μ  (ν = $ν, master_seed = $MASTER_SEED)")
    println("=" ^ 70)
    res = W2G.run_joint_optim(layout, w_z_fourier, basis_supercell;
                              K=K_RESTARTS, master_seed=MASTER_SEED,
                              max_iter=MAX_ITER, ε=ε_TIKHONOV, kw...,
                              targets=(μ > 0 ? targets : nothing), μ, ν,
                              verbose=true,
                              # Per-μ output root: concurrent per-μ SLURM jobs
                              # must never share a run_<timestamp> directory.
                              output_root=joinpath(ROOT, "workflow", "joint_outputs",
                                                   "Ecut$(ECUT)_mu$(μ)"))
    best = res.restart_results[res.best_idx]
    Φs, c = rebuild(best.flat_final)
    # H¹ part alone (comparable across μ): recompute loss without penalty.
    h1² = W2G.joint_loss(best.flat_final, layout, w_z_fourier, basis_supercell;
                         kw..., ε=ε_TIKHONOV)
    rows = W2G.evaluate_hoppings(Φs, c, targets)
    crit = W2G.hopping_criteria(rows, targets)
    results[μ] = (; res, h1², rows, crit)

    @printf("H¹ relative error: %.4f%%   (loss with penalty: %.6e)\n",
            100 * sqrt(h1²), best.loss_final)
    @printf("%-18s %-11s │ %11s %11s │ %11s %11s %9s %9s %5s\n",
            "label", "set", "S_ref", "S_gto", "T_ref", "T_gto",
            "relerr_T", "err/T(0)", "sign")
    for r in rows
        @printf("%-18s %-11s │ %11.3e %11.3e │ %11.3e %11.3e %9.2e %9.2e %5s\n",
                r.label, r.set, r.S_ref, r.S_gto, r.T_ref, r.T_gto,
                r.rel_err_T, r.abs_err_T_per_T0, r.sign_ok ? "ok" : "BAD")
    end
    println("Criteria (independent):")
    for cr in crit
        @printf("  [%s] %-22s value = %.3e\n", cr.pass ? "PASS" : "FAIL",
                cr.name, float(cr.value))
    end
end

# --- Criterion 4: H¹ anti-regression vs the μ=0 twin ----------------------
println("\n", "=" ^ 70)
println("PAIRED SUMMARY (criterion 4: H¹ ≤ $(H1_REGRESSION_FACTOR)× the μ=0 twin)")
println("=" ^ 70)
if !haskey(results, 0.0)
    println("μ = 0 twin not in this job (W2G_MUS split): criterion 4 is computed")
    println("against the H¹ of the μ = 0 job's log when assembling the bilan.")
    for μ in MUS
        r = results[μ]
        c123 = join([cr.pass ? "PASS" : "FAIL" for cr in r.crit], " ")
        @printf("%-10.1e H¹ rel err = %.4f%% │ criteria 1-3: %s\n",
                μ, 100 * sqrt(r.h1²), c123)
    end
    exit(0)
end
h1_twin = results[0.0].h1²
@printf("%-10s %14s %10s %6s │ criteria 1-3\n", "μ", "H¹ rel err", "vs μ=0", "C4")
for μ in MUS
    r = results[μ]
    ratio = sqrt(r.h1² / h1_twin)
    c4 = ratio ≤ H1_REGRESSION_FACTOR
    c123 = join([cr.pass ? "PASS" : "FAIL" for cr in r.crit], " ")
    @printf("%-10.1e %13.4f%% %9.3fx %6s │ %s\n",
            μ, 100 * sqrt(r.h1²), ratio, c4 ? "PASS" : "FAIL", c123)
end
