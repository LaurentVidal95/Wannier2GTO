using Wannier2GTO.GaIn

function overlap(X₁::GaussianPolynomial, X₂::GaussianPolynomial)
    output = zero(X₁.spread)
    @show TR = float(eltype(X₁.center))
    for ((nx₁,ny₁,nz₁), λ₁) in zip(pol_to_arrays(X₁.pol)...)
        for ((nx₂,ny₂,nz₂), λ₂) in zip(pol_to_arrays(X₂.pol)...)
            output += λ₁*λ₂*GaIn.overlap(X₁.spread, TR.(X₁.center), nx₁, ny₁, nz₁,
                                         X₂.spread, TR.(X₂.center), nx₂, ny₂, nz₂)
        end
    end
    output
end

function overlap(Φ₁::BasisFunction, Φ₂::BasisFunction)
#     TR = real(eltype(Φ₁.coeffs))
#     output = zero(TR)
#     for (λ₁, X₁) in Φ₁.
    
end
