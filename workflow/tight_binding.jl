using Wannier2GTO
import Wannier2GTO as W2G

dir = joinpath(splitpath(pathof(Wannier2GTO))[1:end-2]..., "workflow")
include(joinpath(dir, "../tests/produce_graphene_wanniers.jl"))

# Initial data
low_params = (; kgrid=[5,5,1], Ecut=15)
med_params = (; kgrid=[5,5,1], Ecut=30)
high_params = (; kgrid=[5,5,1], Ecut=50)

G = Graphene(; high_params...)
basis = G.basis()
basis_supercell = DFTK.cell_to_supercell(basis)

# Read Wannier from file instead of re-doing wannierization
data = W2G.read_wannier_function(joinpath(dir, "wannier_pz_Ecut-50.json"))

# Import compressed Wannier for GaIn interface.
W₁ = CompressedWannier(basis_supercell, joinpath(dir,
      "../tests/compressed_wannier/compressed_wannier_55_H1_Ecut50.json"))

W₁.center = zeros(3)
positions_cart = map(R -> basis.model.lattice * R,  basis.model.positions)
W₂ = translate(rotate(W₁, π/3), positions_cart[2]);
