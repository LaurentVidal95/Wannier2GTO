@doc raw"""
Matrix of the rotation of angle ``θ`` around ``[0,0]`` in the (x,y)-plane
"""
@inline function rot(θ::T) where {T<:Real}
    [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
end

@inline function polar_to_cartesian_coords(origin, r::T, θ::T) where {T<:Real}
    origin + r*(rot(θ)*[1.; 1.; 0.])
end

"""
Two functions to manadge little complex values
"""
function safereal(tab; tol=1e-10)
    norm(imag(tab)) < tol && return real(tab)
    error("Non negligeable imaginary part")
end

function filter_small_coeffs(z::TC; tol=1e-5) where {TC<:Complex}
    a, b = reim(z) 
    a_new = norm(a)<tol ? zero(a) : a
    b_new = norm(b)<tol ? zero(b) : b
    a_new + im*b_new
end

function filter_dual(x::T) where T
    !(eltype(x) <: AbstractFloat) && (return x.value)
    x
end
