using Wannier2GTO
import Wannier2GTO as W2G

dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]..., "workflow")
include(joinpath(dir,"produce_graphene_wanniers.jl"))

# Initial data
low_params = (; kgrid=[5,5,1], Ecut=15)
med_params = (; kgrid=[5,5,1], Ecut=30)
high_params = (; kgrid=[5,5,1], Ecut=50)

G = Graphene(; low_params...)
basis = G.basis()
basis_supercell = DFTK.cell_to_supercell(basis)

# Produce Wannier
n_bands = 65
scfres = G.scf(n_bands)
wann_model = G.wannierize(scfres; wannier_plot=false,
                          prefix=joinpath(dir, "wannier_output/graphene"))
# data = W2G.prepare_for_compression(wann_model, scfres)

# Read Wannier from file instead of re-doing wannierization
# data = W2G.read_wannier_function(joinpath(dir, "wannier_pz_Ecut-50.json"))

# # Assemble data in a CompressedWannier format
# Wc = CompressedWannier(basis_supercell, normalize(data.wannier), data.center)
# π_bond = data.π_bond

# # Choose H^1 norm for optimization
# s = 1
# Wc.error_norm = s

#Compress given Wannier
# res = compress_graphene_pz_wannier(Wc, π_bond; max_iter=10, file=dir*"outputs/compressed_wannier.json")
