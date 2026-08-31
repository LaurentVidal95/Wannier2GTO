"""
Compute the optimal linear coefficients `c_i` via Tikhonov-regularized
projection onto the basis `Φs` in the H^s-inner-product space.

`s = 1` corresponds to H¹ (the regime used in compression).

This wraps `_make_inner_solver(:tikhonov, ε)` from BasisFunctions.jl. The
matrix `S` and the right-hand side `Γ` are computed in plane-wave Fourier.
Returned `c` has the same length as `Φs`.
"""
function joint_inner_solve(Φs::Vector{BasisFunction},
                           w_z_fourier::AbstractVector,
                           basis_supercell::PlaneWaveBasis;
                           s::Int = 1,
                           ε::Real = 1e-8)
    Φs_Four = [Φ(basis_supercell) for Φ in Φs]
    Γ = [Hˢ_dot(basis_supercell, w_z_fourier, Φ; s=s) for Φ in Φs_Four]
    S = Hˢ_overlap(basis_supercell, Φs_Four; s=s)
    solver = _make_inner_solver(:tikhonov, ε)
    c = solver(S, Γ)
    real(c)   # eigenvalues of S are real, c should be real up to noise
end

"""
Joint-fit loss: H¹ squared residual after Tikhonov-regularized
projection of `w_z` onto the basis described by `flat`.

Parameters
- `flat`: free parameters (length `flat_dim(layout)`)
- `layout`: parameter layout descriptor
- `w_z_fourier`: target wannier in plane-wave Fourier (supercell Γ point)
- `basis_supercell`: discretization basis
- `log_ζ_min`, `log_ζ_max`: spread box bounds (passed to encoding)
- `π_bond_unit`: unit vector along π-bond axis
- `s`: Sobolev exponent (default 1 for H¹)
- `ε`: Tikhonov regularization
- `max_xy_order`, `max_z_order`: SAGTO orders

Returns a positive real scalar.
"""
function joint_loss(flat::AbstractVector, layout::JointLayout,
                    w_z_fourier::AbstractVector,
                    basis_supercell::PlaneWaveBasis;
                    log_ζ_min, log_ζ_max,
                    π_bond_unit::AbstractVector,
                    s::Int = 1,
                    ε::Real = 1e-8,
                    max_xy_order::Int = 3,
                    max_z_order::Int = 3,
                    targets::Union{Nothing, HoppingTargets} = nothing,
                    μ::Real = 0.0,
                    ν::Real = 0.0)
    Φs = params_to_basis_functions(flat, layout;
                                    log_ζ_min, log_ζ_max,
                                    π_bond_unit=π_bond_unit,
                                    max_xy_order, max_z_order)
    Φs_Four = [Φ(basis_supercell) for Φ in Φs]
    Γ = [Hˢ_dot(basis_supercell, w_z_fourier, Φ; s=s) for Φ in Φs_Four]
    S = Hˢ_overlap(basis_supercell, Φs_Four; s=s)
    solver = _make_inner_solver(:tikhonov, ε)
    c = real(solver(S, Γ))
    residual = w_z_fourier - sum(c[i] .* Φs_Four[i] for i in eachindex(c))
    L = Hˢ_norm(basis_supercell, residual; s=s)^2
    if !isnothing(targets) && (μ > 0 || ν > 0)
        L += hopping_penalty(Φs, c, targets; μ, ν)
    end
    L
end

@doc raw"""
Raw (NOT norm-corrected) Gaussian-side hoppings of the fitted function
``w_{\mathrm{GTO}} = \sum_i c_i \Phi_i``:
``S = \langle w, w(\cdot - R)\rangle``, ``T = \langle w, -\tfrac12\nabla^2
w(\cdot - R)\rangle``, via analytic integrals. ForwardDiff-generic in the
parameters carried by `Φs` and `c`.
"""
function gto_hoppings(Φs::AbstractVector{<:BasisFunction}, c::AbstractVector,
                      R::AbstractVector)
    ΦsR = [translate(Φ, R) for Φ in Φs]
    S = zero(eltype(c)); T = zero(eltype(c))
    for i in eachindex(Φs), j in eachindex(ΦsR)
        S += c[i] * c[j] * integral(Φs[i], ΦsR[j]; type=:overlap)
        T += c[i] * c[j] * integral(Φs[i], ΦsR[j]; type=:kinetic)
    end
    (; S, T)
end

@doc raw"""
Hopping penalty of design note 12 §2:
``\mu \sum_{R\in\mathcal{T}} [\Delta T(R)^2 + \Delta S(R)^2]/T(0)^2
 + \nu \sum_{R_{\mathrm{intra}}} S_{\mathrm{GTO}}(R)^2``,
training entries only, norm-corrected GTO values.
"""
function hopping_penalty(Φs::AbstractVector{<:BasisFunction}, c::AbstractVector,
                         targets::HoppingTargets; μ::Real, ν::Real)
    norm² = gto_hoppings(Φs, c, zeros(3)).S
    pen = zero(norm²); ortho = zero(norm²)
    for i in eachindex(targets.labels)
        targets.sets[i] == :training || continue
        h = gto_hoppings(Φs, c, targets.Rs[i])
        Sn, Tn = h.S / norm², h.T / norm²
        pen += ((Tn - targets.T_ref[i])^2 + (Sn - targets.S_ref[i])^2) / targets.T0_ref^2
        (i in targets.ortho_idx) && (ortho += Sn^2)
    end
    μ * pen + ν * ortho
end
