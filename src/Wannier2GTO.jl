module Wannier2GTO

using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff
using DFTK

import DynamicPolynomials: @polyvar
using StaticPolynomials
using ThreadsX

# Visualisation
using Printf
using DelimitedFiles
using WriteVTK

# Utils
export Hs_dot
export Hs_norm
export Hs_overlap
export project_on_AO_basis
export extract_wannier_function
include("common/Hs_scalar_prods.jl")

export store_wannier_function
export read_wannier_function
export find_D3_sym_axis
include("common/utils.jl")
include("common/data_storage.jl")
include("common/extract_wannier_functions.jl")
include("common/callback_info.jl")

# SAGTOs basis
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
