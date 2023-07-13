using Wannier2GTO.GaIn

symb_to_integral = Dict([:overlap => GaIn.overlap,
                         :laplacian => GaIn.laplacian,
                         :kinetic => GaIn.kinetic,
                         :ionic => GaIn.ionic,
                         :coulomb => GaIn.coulomb])

function analytic_norm(exps::Vector{Tuple{Int64, Int64, Int64}},
                       coeffs::Vector{T1}, center::Vector{T2},
                       spread::T3) where {T1, T2, T3 <: Real}
    output = zero(spread)
    TR = float(eltype(center))
    for ((nx₁,ny₁,nz₁), λ₁) in zip(exps, coeffs)
        for ((nx₂,ny₂,nz₂), λ₂) in zip(exps, coeffs)
            output += λ₁*λ₂*symb_to_integral[:overlap](spread, center, nx₁, ny₁, nz₁,
                                                       spread, center, nx₂, ny₂, nz₂)
        end
    end
    √(output)
end

function integral(X₁::GaussianPolynomial, X₂::GaussianPolynomial;
                  type=:overlap)
    output = zero(X₁.spread)
    TR = float(eltype(X₁.center))
    for ((nx₁,ny₁,nz₁), λ₁) in zip(pol_to_arrays(X₁.pol)...)
        for ((nx₂,ny₂,nz₂), λ₂) in zip(pol_to_arrays(X₂.pol)...)
            output += λ₁*λ₂*symb_to_integral[type](X₁.spread, TR.(X₁.center), nx₁, ny₁, nz₁,
                                                   X₂.spread, TR.(X₂.center), nx₂, ny₂, nz₂)
        end
    end
    output
end

function integral(Φ₁::BasisFunction, Φ₂::BasisFunction;
                  type=:overlap)
    TR = real(eltype(Φ₁.coeffs))
    output = zero(TR)
    for (λ₁, X₁) in zip(Φ₁.coeffs, Φ₁.SAGTOs)
        for (λ₂, X₂) in zip(Φ₂.coeffs, Φ₂.SAGTOs)
            output += λ₁*λ₂*integral(X₁, X₂; type)
        end
    end
    output       
end

function integral(Wc₁::CompressedWannier, Wc₂::CompressedWannier;
                  type=:overlap)
    TR = real(eltype(Wc₁.coefficients))
    output = zero(TR)
    for (λ₁, Φ₁) in zip(Wc₁.coefficients, Wc₁.basis_functions)
        for (λ₂, Φ₂) in zip(Wc₂.coefficients, Wc₂.basis_functions)
            output += λ₁*λ₂*integral(Φ₁, Φ₂; type)
        end
    end
    output
end
