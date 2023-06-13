using Wannier2GTO
dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]..., "tests")
include(joinpath(dir,"produce_graphene_wanniers.jl"))

# Initial data
G = Graphene();
basis = G.basis();
basis_supercell = DFTK.cell_to_supercell(basis)

# Read Wannier from file instead of re-doing wannierization
data = Wannier2GTO.read_wannier_function(joinpath(dir, "wannier_pz.json"))
Wc = CompressedWannier(basis_supercell, normalize(data.wannier), data.center)
π_bond = data.π_bond

# Choose H^1 norm for optimization
s = 1
Wc.error_norm = s

# Compress given Wannier
res = compress_graphene_pz_wannier(Wc, π_bond; max_iter=10, file=dir*"compressed_wannier.json")
