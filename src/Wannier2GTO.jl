module Wannier2GTO

# Core packages
using LinearAlgebra
using Optim
using LineSearches
using ForwardDiff   # TODO Replace with ReverseDiff or equiv
using Zygote
using Random
using DFTK
using DFTK.Unitful
using DFTK.UnitfulAtomic

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
export BilayerGraphene
include("CompressedWannier.jl")
include("BilayerGraphene.jl")

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
include("integrals/julia_integrals.jl")
include("integrals/integrals.jl")

# Phase B: joint optimization of SAGTOs
export JointLayout, flat_dim, get_centered, get_pibond, log_ζ_from_u, log_ζ_to_u, params_to_basis_functions, init_params, joint_loss, run_joint_optim
include("joint_optim/parametrization.jl")
include("joint_optim/init.jl")
include("joint_optim/loss.jl")
include("joint_optim/runner.jl")

export TightBindingModel
export R_vectors_cart
export real_hamiltonian
export bloch_transform
include("tight_binding/hamiltonian_scalar_product.jl")
include("tight_binding/TightBindingModel.jl")
include("tight_binding/tb_system.jl")


end # module
