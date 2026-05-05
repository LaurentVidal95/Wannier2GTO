#
# Find which (d, Ecut, kgrid) tuple produces a supercell PlaneWaveBasis whose
# Γ G-vector count matches the precomputed wannier stored in the JSON file.
#

using Wannier2GTO
import Wannier2GTO as W2G
using DFTK
using DFTK.Unitful

include(joinpath(@__DIR__, "monolayer_graphene.jl"))

const JSON_FILE = joinpath(@__DIR__, "wannier_functions/wannier_pz_Ecut-15.json")
const N_TARGET = length(W2G.read_wannier_function(JSON_FILE).wannier)
println("target len(wannier) = ", N_TARGET)

const D_VALUES = (8.0, 10.0, 12.0, 15.0, 20.0)
const ECUT_VALUES = (10, 12, 15, 18, 20)
const KGRID_VALUES = ([3, 3, 1], [4, 4, 1], [5, 5, 1])

for d in D_VALUES, Ecut in ECUT_VALUES, kgrid in KGRID_VALUES
    try
        basis = Graphene(; d=d * u"Å", kgrid=kgrid, Ecut=Ecut).basis()
        sc = DFTK.cell_to_supercell(basis)
        n_G = length(DFTK.G_vectors(sc, sc.kpoints[1]))
        if abs(n_G - N_TARGET) < 200
            tag = (n_G == N_TARGET) ? " <<< EXACT MATCH" : " <<< close"
            println("d=", d, "  Ecut=", Ecut, "  kgrid=", kgrid,
                    "  len=", n_G, tag)
        end
    catch e
        # silently skip failed configs
    end
end

println("scan done")
