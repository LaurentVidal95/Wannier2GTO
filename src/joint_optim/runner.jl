using Optim
using LineSearches
using Printf
using JSON3
using Dates
using Random

"""
Result of a single L-BFGS run.

- `flat_final`: parameter vector at convergence
- `loss_final`: loss value at convergence
- `loss_history`: vector of loss values per iteration
- `iterations`: number of iterations actually run
- `converged`: did Optim report convergence?
- `wallclock_seconds`: wall-clock time
"""
struct LBFGSResult
    flat_final::Vector{Float64}
    loss_final::Float64
    loss_history::Vector{Float64}
    iterations::Int
    converged::Bool
    wallclock_seconds::Float64
end

"""
Run L-BFGS on `joint_loss` from `flat_init`, with ForwardDiff-provided
gradients. ForwardDiff is the validated AD backend for phase B
(see notes/04-plan-phase-B.md, Task 8 v2). Switch to Enzyme can be done
in a future iteration if performance becomes a bottleneck.

Hyperparameters follow `notes/03-design-phase-B.md` §3 defaults.
"""
function run_lbfgs_once(flat_init::AbstractVector, layout::JointLayout,
                        w_z_fourier::AbstractVector,
                        basis_supercell::PlaneWaveBasis;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit::AbstractVector,
                        max_iter::Int = 200,
                        g_abstol::Real = 1e-5,
                        f_reltol::Real = 1e-8,
                        s::Int = 1,
                        ε::Real = 1e-8,
                        targets::Union{Nothing, HoppingTargets} = nothing,
                        μ::Real = 0.0,
                        ν::Real = 0.0,
                        verbose::Bool = false)
    f = θ -> joint_loss(θ, layout, w_z_fourier, basis_supercell;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit=π_bond_unit, s=s, ε=ε,
                        targets, μ, ν)
    g! = (G, θ) -> (G .= ForwardDiff.gradient(f, θ); nothing)

    history = Float64[]
    t_start = time()
    cb = function (state)
        push!(history, state.value)
        if verbose
            # Live progress for batch logs (flushed: SLURM buffers stdout).
            println("    iter ", lpad(state.iteration, 3), "  loss = ",
                    @sprintf("%.6e", state.value), "  elapsed = ",
                    round(time() - t_start; digits=1), "s")
            flush(stdout)
        end
        false   # keep going
    end

    options = Optim.Options(g_abstol=g_abstol, f_reltol=f_reltol,
                            iterations=max_iter, callback=cb,
                            show_trace=false)

    t0 = time()
    res = optimize(f, g!, collect(flat_init),
                   LBFGS(; linesearch=LineSearches.HagerZhang()),
                   options)
    elapsed = time() - t0

    LBFGSResult(Optim.minimizer(res), Optim.minimum(res),
                history, Optim.iterations(res),
                Optim.converged(res), elapsed)
end

"""
Aggregate result of a K-restart joint optimization run.
"""
struct JointOptimResult
    layout::JointLayout
    K::Int
    seeds::Vector{Int}
    restart_results::Vector{LBFGSResult}
    best_idx::Int
    output_dir::String
end

"""
Run K independent L-BFGS restarts, save artifacts, return the best result.

Output layout under `output_root/run_<timestamp>/`:
  - config.json     : hyperparameters and seeds
  - restart_<k>.json: per-restart trace (loss curve + final flat)
  - best.json       : best restart payload
  - summary.txt     : human-readable comparison
"""
function run_joint_optim(layout::JointLayout,
                         w_z_fourier::AbstractVector,
                         basis_supercell::PlaneWaveBasis;
                         K::Int = 5,
                         master_seed::Int = 42,
                         log_ζ_min, log_ζ_max,
                         π_bond_unit::AbstractVector,
                         output_root::String = joinpath(@__DIR__,
                                                        "..", "..",
                                                        "workflow", "joint_outputs"),
                         max_iter::Int = 200,
                         g_abstol::Real = 1e-5,
                         f_reltol::Real = 1e-8,
                         s::Int = 1,
                         ε::Real = 1e-8,
                         targets::Union{Nothing, HoppingTargets} = nothing,
                         μ::Real = 0.0,
                         ν::Real = 0.0,
                         verbose::Bool = false)
    timestamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    output_dir = joinpath(output_root, "run_" * timestamp)
    mkpath(output_dir)

    seeds = [master_seed + k for k in 0:(K - 1)]
    results = LBFGSResult[]

    for (k, seed) in enumerate(seeds)
        println("--- Restart $k/$K (seed=$seed) ---"); flush(stdout)
        flat0 = init_params(MersenneTwister(seed), layout)
        res = run_lbfgs_once(flat0, layout, w_z_fourier, basis_supercell;
                             log_ζ_min, log_ζ_max,
                             π_bond_unit=π_bond_unit,
                             max_iter, g_abstol, f_reltol, s, ε,
                             targets, μ, ν, verbose)
        push!(results, res)
        println("    restart done: loss_init=", first(res.loss_history),
                "  loss_final=", res.loss_final, "  iters=", res.iterations,
                "  wallclock=", round(res.wallclock_seconds; digits=1), "s")
        flush(stdout)
        @info "  loss_init=$(first(res.loss_history))  loss_final=$(res.loss_final)" *
              "  iters=$(res.iterations)  wallclock=$(round(res.wallclock_seconds; digits=1))s"
        # Per-restart persistence
        open(joinpath(output_dir, "restart_$(lpad(k, 2, '0')).json"), "w") do io
            JSON3.write(io,
                Dict("seed" => seed,
                     "loss_history" => res.loss_history,
                     "loss_final" => res.loss_final,
                     "iterations" => res.iterations,
                     "converged" => res.converged,
                     "wallclock_seconds" => res.wallclock_seconds,
                     "flat_final" => res.flat_final))
        end
    end

    best_idx = argmin([r.loss_final for r in results])

    # Aggregate config and best
    open(joinpath(output_dir, "config.json"), "w") do io
        JSON3.write(io,
            Dict("K" => K, "master_seed" => master_seed,
                 "seeds" => seeds, "max_iter" => max_iter,
                 "log_ζ_min" => log_ζ_min, "log_ζ_max" => log_ζ_max,
                 "ε" => ε, "s" => s,
                 "N_centered" => layout.N_centered,
                 "N_pibond" => layout.N_pibond,
                 "poly_dim_centered" => layout.poly_dim_centered,
                 "poly_dim_pibond" => layout.poly_dim_pibond,
                 "μ" => μ, "ν" => ν,
                 "targets_file" => "see caller"))
    end

    open(joinpath(output_dir, "best.json"), "w") do io
        JSON3.write(io,
            Dict("best_restart" => best_idx,
                 "best_seed" => seeds[best_idx],
                 "best_loss" => results[best_idx].loss_final,
                 "best_flat" => results[best_idx].flat_final))
    end

    open(joinpath(output_dir, "summary.txt"), "w") do io
        println(io, "Joint optim summary — $(timestamp)")
        println(io, "K = $K, layout = (N_centered=$(layout.N_centered), N_pibond=$(layout.N_pibond))")
        println(io, "")
        println(io, "Restart  Seed     Loss_init      Loss_final     Iters  Wall(s)")
        for (k, r) in enumerate(results)
            println(io, lpad(k, 6), "  ",
                    lpad(seeds[k], 5), "  ",
                    lpad(round(first(r.loss_history); digits=4), 12), "  ",
                    lpad(round(r.loss_final; digits=6), 14), "  ",
                    lpad(r.iterations, 5), "  ",
                    lpad(round(r.wallclock_seconds; digits=1), 7))
        end
        println(io, "")
        println(io, "Best: restart $best_idx with loss = $(results[best_idx].loss_final)")
    end

    JointOptimResult(layout, K, seeds, results, best_idx, output_dir)
end
