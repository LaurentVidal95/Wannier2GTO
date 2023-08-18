module Wannier2GTO

# Core packages
using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff   # TODO Replace with ReverseDiff or equiv
using DFTK

# Fast evaluation of multivariate polynomials
import DynamicPolynomials: @polyvar
using StaticPolynomials
using ThreadsX

# I/O stuff and Visualisation
using Printf

# Basis functions are linear combinations of gaussian polynomials
export GaussianPolynomial
export SAGTO_basis
include("basis_functions/GaussianPolynomials.jl")
include("basis_functions/SAGTO.jl")
include("basis_functions/BasisFunctions.jl")

# Struct to store the compressed wannier
export CompressedWannier
export translate
export rotate
include("CompressedWannier.jl")

# Common stuff
include("common/Hs_scalar_prods.jl")
include("common/utils.jl")
include("common/callback_info.jl")

# Only there to gain time in debugging by skipping the
# SCF and wannierization.
include("common/data_storage.jl")

# Compression routine
export compress_graphene_pz_wannier
include("compression/wannier_preprocessing.jl")
include("compression/compression.jl")

# Gaussian polynomial analytic
export integral
include("integrals/GaIn.jl")
include("integrals/integrals.jl")

export real_hamiltonian_bloc
include("tight_binding/hamiltonian.jl")

end # module
