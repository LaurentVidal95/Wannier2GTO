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

@doc raw"""
For a given table `[u_G]` and linear operator ``L \colon G \to L(G)``, compute `[u_{L(G)}]`.
BEWARE: L has to act on reduced coordinates.
If L is defined in cartesian coordinate apply to ``L' = L·\textrm{cart_to_red}``.
"""
function apply_linear_operator(new_indexes, f_in::AbstractArray)
    f_out = zero(f_in)
    # Filter G_vectors that are maped outside the entry basis.    
    iG  = [index[1] for index in new_indexes if !isnothing(index[2])]
    iLG = [index[2] for index in new_indexes if !isnothing(index[2])]
    f_out[iLG] .= f_in[iG]
    f_out
end
function apply_linear_operator(basis::PlaneWaveBasis{T1}, kpt::Kpoint,
                               u::Vector{T2}, L::F) where {T1<:Real, T2<:Complex, F}
    new_indexes = [(iG, DFTK.index_G_vectors(basis, kpt, L(G))) for (iG, G) in
                   enumerate(G_vectors(basis, kpt))]
    apply_linear_operator(new_indexes, u)
end
function apply_linear_operator(basis::PlaneWaveBasis{T1}, G_vec::AbstractArray,
                               u::AbstractArray{T2}, L::F) where {T1<:Real, T2<:Complex, F}
    new_indexes = [(iG, DFTK.index_G_vectors(basis, L(G))) for (iG, G) in enumerate(G_vec)]
    apply_linear_operator(new_indexes, u)
end

# ParityOp(u, signs) = apply_linear_operator(basis, K, u, G -> signs .*G)
