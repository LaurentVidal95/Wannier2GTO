using Optim
using LineSearches

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
                        ε::Real = 1e-8)
    f = θ -> joint_loss(θ, layout, w_z_fourier, basis_supercell;
                        log_ζ_min, log_ζ_max,
                        π_bond_unit=π_bond_unit, s=s, ε=ε)
    g! = (G, θ) -> (G .= ForwardDiff.gradient(f, θ); nothing)

    history = Float64[]
    cb = function (state)
        push!(history, state.value)
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
