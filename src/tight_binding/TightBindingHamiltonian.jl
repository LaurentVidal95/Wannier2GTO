struct TightBindingHamiltonian{T<:Complex}
    terms::Vector{Symbol}
    potential::AbstractArray{T}
end

function hamiltonian(terms::Vector{Symbol}, potential::AbstractArray=Complex[])
    if isempty(potential) && (:PBE_potential ∈ terms)
        error("The KS potential must be provided to compute the PBE potential term")
    end
    TightBindingHamiltonian(terms, potential)
end

function apply_term(term::Symbol, args...)
    is_potential_term = (term==:PBE_potential)
    if (length(args)<3) && is_potential_term
        error("The KS potential is required to compute the :PBE_potential term contribution")
    end

    W₁, W₂ = args[1:2]

    if !is_atomic_term
        return integral(W₁, W₂; type=term)
    else
        
    end
    error("Well, that was inexpected")
end

function hamiltonian_scalar_prod(H, f1, f2)
    sum(apply_term(term, H.KS_potential) for term in H.terms)
end
