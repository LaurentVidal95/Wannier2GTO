@doc raw"""
Construct one s orbital of center ``α`` and spread ``ζ``.
"""
s_orb(α, ζ) = GaussianPolynomial([(0,0,1)], [1.], α, ζ)

@doc raw"""
Matrix of the rotation of angle ``θ`` around ``[0,0]`` in the (x,y)-plane
"""
rot(θ) = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]

polar_to_cartesian_coords(origin, λ, θ) = origin .+ λ .*(rot(θ)*[1.; 1.; 0])
