module Wannier2GTO

# Core packages
using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff
using DFTK

# Fast evaluation of multivariate polynomials
import DynamicPolynomials: @polyvar
using StaticPolynomials
using ThreadsX

# Visualisation
using Printf
using DelimitedFiles
using WriteVTK

export CompressedWannier
include("CompressedWannier.jl")

# Common features
include("common/Hs_scalar_prods.jl")
include("common/utils.jl")
include("common/callback_info.jl")

# SAGTOs basis
export GaussianPolynomial
export SAGTO_basis
include("basis_functions/GaussianPolynomials.jl")
include("basis_functions/SAGTO.jl")
include("basis_functions/BasisFunctions.jl")

# Compression
export compress_graphene_pz_wannier
include("wannier_preprocessing.jl")
include("compression.jl")

end # module
