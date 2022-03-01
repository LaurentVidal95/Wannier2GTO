module GrapheneTB

using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff
using DFTK
using Printf
using DelimitedFiles

# Utils
export Hs_scalar_prod
export Hs_overlap
export Hs_projection_on_AO_basis
include("common/Hs_scalar_prods.jl")
include("common/utils.jl")
include("common/read_data.jl")

# SAGTOs
export GaussianPolynomial
export SAGTOs_basis
include("SAGTOs/GaussianPolynomials.jl")
include("SAGTOs/symmetry_adapted_polynoms.jl")
include("SAGTOs/SAGTOs_basis.jl")

# SCF and Wannierization
export scf_graphene
export run_wannier90_graphene
include("compute_graphene_wannier.jl")

# Compression
include("compression_routines/cost_function.jl")
# include("compression_routines/compression.jl")
    
end # module
