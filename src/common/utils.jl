@doc raw"""
Matrix of the rotation of angle ``θ`` around ``[0,0]`` in the (x,y)-plane
"""
@inline function rot(θ::T) where {T<:Real}
    [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
end

@inline function polar_to_cartesian_coords(origin, r::T, θ::T) where {T<:Real}
    origin + r*(rot(θ)*[1.; 1.; 0.])
end
