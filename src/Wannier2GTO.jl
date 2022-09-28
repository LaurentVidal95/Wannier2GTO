module Wannier2GTO

using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff
using DFTK

using Printf
using DelimitedFiles
using WriteVTK

# Utils
export Hs_scalar_prod
export Hs_overlap
export Hs_projection_on_AO_basis
export Hs_norm
export cell_to_supercell_vector
export extract_wannier_function
export store_wannier_function
export find_D3_sym_axis
include("common/Hs_scalar_prods.jl")
include("common/utils.jl")
include("common/data_storage.jl")
include("common/extract_wannier_functions.jl")
include("common/supercell.jl")
include("common/callback_info.jl")

# SAGTOs
export GaussianPolynomial
export SAGTOs_basis
include("SAGTOs/GaussianPolynomials.jl")
include("SAGTOs/symmetry_adapted_polynoms.jl")
include("SAGTOs/SAGTOs_basis.jl")

# Compression
export compress_graphene_pz_wannier
include("compression_routines/inner_optimization.jl")
include("compression_routines/compression.jl")
    
end # module
