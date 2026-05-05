#
# Debug routines that allow to store a wannier function
# to avoid redoing the wannierization each time
#

"""
Allow to store complex numbers with JSON3
"""
struct JCX
    real::Float64
    imaginary::Float64
end
JCX(x) = JCX(reim(x)...)

function store_wannier_function(W::Vector{ComplexF64}, center::Vector{T}, r::T, θ::T,
                                 prefix::String) where {T<:Real}
    @assert size(W,2)==1  "The wannier function is to be given as a single supercell vector"
    # Check if Ws are unitary
    !(norm(W) ≈ 1) && (@warn "Provided Wannier is not normalized")
    # Convert wanniers to complex JSON
    W_dict = Dict{String, Any}()
    W_dict["wannier"] = JCX.(W)
    W_dict["center"] = center
    W_dict["π_bond"] = [r, θ]
    open(io->JSON3.write(io, W_dict, allow_inf=true), prefix*".json", "w")
    nothing
end

function read_wannier_function(filename::String)
    data = open(JSON3.read, filename)
    JCX_to_complex(x) = x.real .+ im*x.imaginary
    wannier = JCX_to_complex.(data["wannier"])
    center = Float64.(data["center"])
    r, θ = Float64.(data["π_bond"])
    !(norm(wannier) ≈ 1) && @warn "The wannier function is not normalized"
    (;wannier, center, π_bond=(;r, θ))
end
