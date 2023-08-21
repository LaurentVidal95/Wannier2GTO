struct BilayerGraphene{T<:Real}
    lattice         ::AbstractArray
    atoms           ::AbstractVector{DFTK.Element}
    positions       ::AbstractVector
    # Angle of rotation and disregistry between layers A and B
    angle           ::T
    disregistry     ::AbstractVector{T}
end
function monolayer_graphene_model(; lattice_constant =  auconvert(2.641u"Å"),
                           z_axis_length = 40u"bohr")
    lz = z_axis_length / lattice_constant
    lattice = lattice_constant .* [1  -1/2    0;
                                   0  √3/2    0;
                                   0     0    lz]
    positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
    atoms = fill(ElementPsp(6, psp=load_psp("hgh/pbe/c-q4")), length(positions))
    DFTK.Model(lattice, atoms, positions)
end

function StackedBilayerGraphene(disregistry;
                                interlayer_distance=6.45u"bohr", kwargs...)
    # Compute monolayer geomtry
    mono = monolayer_graphene_model(; kwargs...)
    lattice = mono.lattice

    # Compute atomic positions for each layer
    d_red = interlayer_distance / lattice[end]
    # Remove unit for next computations if needed
    (typeof(d_red) <: Quantity) && (d_red = d_red.val)
    positions_A = mono.positions .- Ref([0., 0., d_red/2])
    positions_B = mono.positions .+ Ref([disregistry..., d_red/2])

    # Gather both layers
    atoms = mono.atoms
    append!(atoms, atoms)
    positions = [positions_A..., positions_B...]

    BilayerGraphene(lattice, atoms, positions, zero(eltype(disregistry)), disregistry)
end

function Base.show(io::IO, BG::BilayerGraphene)    
    name = iszero(BG.angle) ? "Stacked Bilayer Graphene" : "Twisted Bilayer Graphene"
    interlayer_distance = abs(BG.positions[1][3] - BG.positions[end][3]) * BG.lattice[end]
    
    println(io, name*"\n")
    println(io, "Interlayer distance : ", interlayer_distance*u"bohr")
    println(io, "Disregistry : $(BG.disregistry)")
    println(io, "Angle of rotation between layers : $(BG.angle)")
    println(io, "Periodic system : ", iszero(BG.angle))
end

function TwistedBilayerGraphene(disregistry, angle;
                                interlayer_distance=6.45u"bohr", kwargs...)
    # positions_B = map(R->rot(angle)*R + [disregistry..., d_red/2],
    #                   mono.positions)
    nothing
end

function BilayerGraphene(disregistry, angle;
                         interlayer_distance=6.45u"bohr", kwargs...)
    iszero(angle) && return StackedBilayerGraphene(disregistry; interlayer_distance, kwargs...)
    TwistedBilayerGraphene(disregistry, angle; interlayer_distance, kwargs...)
end
