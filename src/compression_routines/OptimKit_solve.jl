"""
Equivalent of the prompt routine but with OptimKit conventions
"""
function finalize!(ζ, E, ∇E, n_iter)
    if n_iter == 1
        println("ROHF direct energy minimization")
        println("Initial guess: $(ζ.guess)")
        header = ["Iter", "Energy","log10(ΔE)", "log10(||Π∇E||)"]
        println("-"^58)
        println(@sprintf("%-5s  %-16s  %-16s  %-16s", header...))
        println("-"^58)

        info_out = [n_iter-1, ζ.energy, " "^16, " "^16]
        println(@sprintf("%5i %16.12f %16s %16s", info_out...))
        flush(stdout)
    end
    # Print current iter infos
    log_ΔE = log(10, abs(E  - ζ.energy))
    residual = norm(∇E)
    info_out = [n_iter, E, log_ΔE, log10(residual)]
    println(@sprintf("%5i %16.12f %16.12f %16.12f", info_out...))
    flush(stdout)

    # Actualize energy
    energy!(ζ)
    # Update history
    ζ.history = vcat(ζ.history, reshape(info_out, 1, 4))

    # Return entry to match OptimKit.jl conventions
    ζ, E, ∇E
end
