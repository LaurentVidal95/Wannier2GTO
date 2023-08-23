function apply_term(term::Symbol, args...)
    is_atomic_term = term==:atomic
    if (length(args)<3) && is_atomic_term
        error("The position of the atom is required for atom - electron interaction")
    end

    W₁, W₂ = args[1:2]

    if !is_atomic_term
        return integral(W₁, W₂; type=term)
    else
        atomic_positions = args[3]
        Z_carbon = 6
        output = sum(map(atomic_positions) do R
                     integral(W₁, W₂, R; type=term)
                     end
                     )
        return -Z_carbon*output
    end
    error("Well, that was inexpected")
end
function hamiltonian_scalar_prod(TB::TightBindingModel, args...)
    sum(apply_term(term, args..., TB.atomic_positions)
        for term in TB.hamiltonian_terms)
end
