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
