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

# SAGTOs basis
export GaussianPolynomial
export SAGTOs_basis
include("SAGTOs/GaussianPolynomials.jl")
include("SAGTOs/symmetry_adapted_polynoms.jl")
include("SAGTOs/SAGTOs_basis.jl")

# Compression
export compress_graphene_pz_wannier
include("wannier_preprocessing.jl")
include("compression_routines/inner_optimization.jl")
include("compression_routines/compression.jl")

# Common features
export Hs_dot
export Hs_norm
export Hs_overlap
export project_on_AO_basis
include("common/Hs_scalar_prods.jl")
include("common/utils.jl")
include("common/callback_info.jl")

end # module
