mutable struct CompressedWannier{TR<:Real, TC<:Complex}
    # Original Wannier function to compress
    basis :: PlaneWaveBasis{TR}
    basis_supercell :: PlaneWaveBasis{TR}
    wannier :: AbstractArray{TC}
    center :: Vector{TR}
    # SAGTO basis
    SAGTO_basis ::Vector{GaussianPolynomial{TR}}
    MO_coefficients :: Vector{AbstractArray{TR}}
    # Optim related
    residual :: AbstractArray{TC}
    residual_norm :: TR
    grad_norm :: TR
end

# function CompressedWannier(basis::PlaneWaveBasis, Wn::AbstractArray)
#     # Express Wnk as a single vector on the supercell    
#     basis_SC = cell_to_supercell(basis)
#     Wn_SC = cell_to_supercell(Wnk, basis, basis_SC)
#     # Initiate other args
#     TR = eltype(basis_SC)
#     SAGTO_basis = GaussianPolynomial{eltype(basis_SC)}[]
#     MO_coefficients = 
# end
