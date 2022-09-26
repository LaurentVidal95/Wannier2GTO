@doc raw"""
Construct one s orbital of center ``α`` and spread ``ζ``.
"""
s_orb(α, ζ) = GaussianPolynomial([(0,0,1)], [1.], α, ζ)

@doc raw"""
Matrix of the rotation of angle ``θ`` around ``[0,0]`` in the (x,y)-plane
"""
rot(θ) = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]

polar_to_cartesian_coords(origin, λ, θ) = origin .+ λ .*(rot(θ)*[1.; 1.; 0])

"""
Provides the angle between axis x and the first axis of the (x,y)-plane
D3 symmetry of pz-like wannier functions.
"""
function find_D3_sym_axis(basis_SC::PlaneWaveBasis, W_pz, α0)
    α(λ, θ) = polar_to_cartesian_coords(α0, λ, θ)
    res = optimize(X->norm(W_pz .- s_orb(α(X[1],X[2]), X[3])(basis_SC)),
                   [1., -1/2, 1/2], # Guess roughly close to wanted axis by experience.
                   ConjugateGradient(linesearch=BackTracking(order=3)),
                   Optim.Options(show_trace=true))
    res.minimizer[2]
end
