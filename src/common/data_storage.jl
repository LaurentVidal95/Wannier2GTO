# TODO: Update to new Wannier2GTO

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

# OLD

# function write_wannier_vtk(wn, basis::PlaneWaveBasis, prefix::String;
#                            basis_SC=cell_to_supercell(basis))
#     wn_SC = cell_to_supercell(wn, DFTK.unfold_bz(basis), basis_SC)
#     wn_real = real.(ifft(basis_SC, only(basis_SC.kpoints), sum(eachcol(wn_SC))))
#     r_cart = map(r->basis_SC.model.lattice*r, basis_SC.r_vectors)
#     x = [r[1] for r in r_cart]
#     y = [r[2] for r in r_cart]
#     z = [r[3] for r in r_cart]
#     vtk_grid(prefix, x, y, z) do vtk
#         vtk["wannier"] = wn_real
#     end
#     nothing
# end

# function read_compressed_wannier(basis_SC, filename::String)
#     data = open(JSON3.read, filename)
#     # Extract proj coeffs
#     proj_coeffs = Float64.(data[0])
#     N_SAGTOs = length(proj_coeffs)
#     p = Progress(N_SAGTOs, desc="Assembling compressed Wannier")
#     # Extract all MOs
#     Φs = map(1:N_SAGTOs) do i
#         Λs, αs, ζs, pol_orders = data[i]
#         next!(p) # progress bar
#         assemble_MO(basis_SC, Λs, αs, ζs, pol_orders...)
#     end
#     proj_coeffs'Φs
# end

# function assemble_MO(basis_SC, Λ, αs, ζs, xy_orders, z_orders)
#     # JSON type -> Float64
#     Λ = Float64.(Λ)
#     αs = [Float64.(α) for α in αs]
#     ζs = Float64.(ζs)

#     # Assemble MO
#     AOs = vcat(SAGTOs_basis.(αs, Ref(ζs), Ref(xy_orders), Ref(z_orders))...)
#     Λ'ThreadsX.map(X->X(basis_SC), AOs)
# end
