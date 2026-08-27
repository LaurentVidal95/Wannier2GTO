function apply_term(term::Symbol, W₁::CompressedWannier, W₂::CompressedWannier, scfres=nothing)
    KS_potential_term = (term==:PBE_local_potential)
    if isnothing(scfres) && KS_potential_term
        error("The KS potential is required to compute the :PBE_potential term contribution")
    end
    
    if !KS_potential_term
        return integral(f₁, f₂; type=term)
    else
        V = extract_local_KS_potential(scfres)
        return potential_scalar_prod(scfres.basis, V, W₁, W₂)
    end
    error("Well, that was inexpected")
end

function hamiltonian_scalar_prod(BG::BilayerGraphene, W₁::CompressedWannier,
                                 W₂::CompressedWannier; terms=[:kinetic, :PBE_local_potential])
    sum(apply_term(term, W₁, W₂, BG.scfres) for term in terms)
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{T},
                               X₁::GaussianPolynomial, X₂::GaussianPolynomial) where {T<:Complex}
    Ω = basis.model.unit_cell_volume
    X₃ = X₁*X₂
    X₃_four = ℱ.(Ref(X₃), G_vectors_cart(basis))
    (1/√Ω)*dot(X₃_four,V_KS_four)
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{TC},
                               Φ₁::BasisFunction, Φ₂::BasisFunction) where {TC<:Complex}
    output = zero(TC)
    for (λ₁, X₁) in zip(Φ₁.coeffs, Φ₁.SAGTOs)
        for (λ₂, X₂) in zip(Φ₂.coeffs, Φ₂.SAGTOs)
            @show output += λ₁*λ₂*potential_scalar_prod(basis, V_KS_four, X₁, X₂)
        end
    end
    output
end

function potential_scalar_prod(basis::PlaneWaveBasis, V_KS_four::AbstractArray{TC},
                               Wc₁::CompressedWannier, Wc₂::CompressedWannier) where {TC<:Complex}
    output = zero(TC)
    for (λ₁, Φ₁) in zip(Wc₁.coefficients, Wc₁.basis_functions)
        for (λ₂, Φ₂) in zip(Wc₂.coefficients, Wc₂.basis_functions)
            output += λ₁*λ₂*potential_scalar_prod(basis, V_KS_four, Φ₁, Φ₂)
        end
    end
    output
end

@doc raw"""
On-site one-body kinetic comparison, first link of the tight-binding validation
chain: contrast the kinetic integral obtained from the *true plane-wave*
Wannier against the one obtained from its *Gaussian-compressed* counterpart.

- reference (grid): ``\langle w, -\tfrac12\nabla^2 w\rangle
  = \tfrac12\langle w, |G|^2 w\rangle`` — on the plane-wave grid ``\nabla^2``
  is multiplication by ``-|G|^2``, so no new machinery is needed.
- Gaussian (analytic): `integral(Wc, Wc; type=:kinetic)`, native since
  `kinetic_julia`.

`w_fourier` and `Wc` must live on the same `basis_supercell` (Γ point of the
supercell, DFTK convention).

The compressed Wannier only approximates ``w``, so its norm is not exactly 1
while the reference is normalized. Raw values therefore conflate a shape error
with a norm error; the `_per_norm` fields divide by ``\|w\|^2`` and are the
physically meaningful on-site tight-binding element. Both are returned.
"""
function compare_onsite_kinetic(Wc::CompressedWannier, w_fourier::AbstractVector,
                                basis_supercell::PlaneWaveBasis)
    kpt = only(basis_supercell.kpoints)
    G² = [sum(abs2, Gpk) for Gpk in G_vectors_cart(basis_supercell, kpt)]

    T_ref = 0.5 * real(dot(w_fourier, G² .* w_fourier))
    T_gto = integral(Wc, Wc; type=:kinetic)

    norm²_ref = real(dot(w_fourier, w_fourier))
    norm²_gto = real(integral(Wc, Wc; type=:overlap))
    T_ref_per_norm = T_ref / norm²_ref
    T_gto_per_norm = T_gto / norm²_gto

    (; T_ref, T_gto,
       rel_err = abs(T_gto - T_ref) / abs(T_ref),
       T_ref_per_norm, T_gto_per_norm,
       rel_err_per_norm = abs(T_gto_per_norm - T_ref_per_norm) / abs(T_ref_per_norm),
       norm_ref = sqrt(norm²_ref), norm_gto = sqrt(norm²_gto))
end
