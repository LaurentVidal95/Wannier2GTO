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

@doc raw"""
Hopping comparison, second link of the tight-binding validation chain:
``S(R) = \langle w_0, w_R\rangle`` and
``T(R) = \langle w_0, -\tfrac12\nabla^2 w_R\rangle`` from the true plane-wave
Wannier vs its Gaussian-compressed counterpart. `R_cart` is Cartesian, Bohr.

- reference: exact periodic translation ``w(\cdot-R) \leftrightarrow
  w_G e^{-iG\cdot R}``, then grid dot products (no new machinery);
- Gaussian: `translate(Wc, R_cart)` then native `integral`.

The hoppings of a real Wannier are real; the imaginary parts of the reference
values are numerical leakage and asserted small, not returned.

Beware small denominators on distant hoppings: `rel_err_*` divides by the
reference value, `abs_err_*_per_T0` divides by the on-site kinetic ``T(0)`` —
report both (design note 09).
"""
function compare_hopping(Wc::CompressedWannier, w_fourier::AbstractVector,
                         basis_supercell::PlaneWaveBasis, R_cart::AbstractVector)
    @assert length(R_cart) == 3 "R_cart must be a 3-vector, got $(length(R_cart))"
    @assert all(isfinite, w_fourier) "non-finite coefficients in the reference Wannier"
    kpt = only(basis_supercell.kpoints)
    Gs = G_vectors_cart(basis_supercell, kpt)
    @assert length(w_fourier) == length(Gs) "Wannier/basis mismatch"

    phase = [cis(-dot(G, R_cart)) for G in Gs]
    wR = w_fourier .* phase
    G² = [sum(abs2, G) for G in Gs]

    S_ref_c = dot(w_fourier, wR)                 # ⟨w₀, w_R⟩
    T_ref_c = 0.5 * dot(w_fourier, G² .* wR)     # ⟨w₀, -½∇² w_R⟩
    norm²_ref = real(dot(w_fourier, w_fourier))
    T0_ref = 0.5 * real(dot(w_fourier, G² .* w_fourier))  # on-site scale

    IMAG_TOL = 1e-6  # real-Wannier sanity: imaginary leakage bound (relative)
    @assert abs(imag(S_ref_c)) ≤ IMAG_TOL * norm²_ref "Im S(R) leakage: $(imag(S_ref_c))"
    @assert abs(imag(T_ref_c)) ≤ IMAG_TOL * abs(T0_ref) "Im T(R) leakage: $(imag(T_ref_c))"
    S_ref, T_ref = real(S_ref_c), real(T_ref_c)

    WcR = translate(Wc, R_cart)
    S_gto = integral(Wc, WcR; type=:overlap)
    T_gto = integral(Wc, WcR; type=:kinetic)
    norm²_gto = real(integral(Wc, Wc; type=:overlap))

    # Norm-corrected values (translation is unitary: ‖w_R‖ = ‖w₀‖).
    S_ref_n, T_ref_n = S_ref / norm²_ref, T_ref / norm²_ref
    S_gto_n, T_gto_n = S_gto / norm²_gto, T_gto / norm²_gto

    (; S_ref = S_ref_n, S_gto = S_gto_n,
       T_ref = T_ref_n, T_gto = T_gto_n,
       T0_ref = T0_ref / norm²_ref,
       rel_err_S = abs(S_gto_n - S_ref_n) / abs(S_ref_n),
       rel_err_T = abs(T_gto_n - T_ref_n) / abs(T_ref_n),
       abs_err_S_per_T0 = abs(S_gto_n - S_ref_n) / abs(T0_ref / norm²_ref),
       abs_err_T_per_T0 = abs(T_gto_n - T_ref_n) / abs(T0_ref / norm²_ref),
       norm_ref = sqrt(norm²_ref), norm_gto = sqrt(norm²_gto))
end
