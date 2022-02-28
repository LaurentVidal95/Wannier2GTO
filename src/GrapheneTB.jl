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

# SAGTOs
export GaussianPolynomial
export construct_pz_SAGTO_basis
include("atomic_orbitals/GaussianPolynomials.jl")
include("atomic_orbitals/SAGTOs_parameters.jl")

# SCF and Wannierization
export scf_graphene
export run_wannier90_graphene
include("compute_graphene_wannier.jl")

# Compression
# include("compression_routines/cost_function.jl")
# include("compression_routines/compression.jl")
    
end # module
