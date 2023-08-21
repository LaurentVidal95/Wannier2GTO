using Wannier2GTO
import Wannier2GTO as W2G

dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]..., "workflow")
include(joinpath(dir, "monolayer_graphene.jl"))

# Initial data
low_params = (; kgrid=[5,5,1], Ecut=15)
med_params = (; kgrid=[5,5,1], Ecut=30)
high_params = (; kgrid=[5,5,1], Ecut=50)

BG = BilayerGraphene(zeros(2), zero(Float64))
W₁, W₂ = W2G.TightBindingModel(BG, 1e-8).basis_functions

# Read Wannier from file instead of re-doing wannierization
data = W2G.read_wannier_function(joinpath(dir, "wannier_functions/wannier_pz_Ecut-50.json"))
