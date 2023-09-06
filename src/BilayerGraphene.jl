struct BilayerGraphene{T<:Real}
    lattice         ::AbstractArray
    atoms           ::AbstractVector{DFTK.Element}
    positions       ::AbstractVector
    # Angle of rotation and disregistry between layers A and B
    angle           ::T
    disregistry     ::AbstractVector{T}
end
geometry(BG::BilayerGraphene) = (BG.lattice, BG.atoms, BG.positions)

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

    println(io, name)
    println(io, "- interlayer distance : ", interlayer_distance*u"bohr")
    println(io, "- disregistry : $(BG.disregistry) (reduced coordinates)")
    println(io, "- angle of rotation between layers : $(BG.angle)°")
    println(io, "- periodic system : ", iszero(BG.angle))
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

const default_convergence_parameters = (;Ecut=30, kgrid=[5,5,1], scf_tol=1e-5)
function compute_scfres(BG::BilayerGraphene;
                        terms=nothing,
                        convergence_parameters = default_convergence_parameters,
                        # Temperature
                        temperature=1e-3,
                        smearing=DFTK.Smearing.Gaussian(),
                        # To be able to ask for a specific fft grid
                        basis_kwargs...
                        )
    # Build model (PBE unless terms are given)
    model = nothing
    if isnothing(terms)
        model = model_PBE(geometry(BG)...; temperature, smearing)
    else
        model = Model(geometry(BG)...; temperature, smearing, terms)
    end

    # Build PW basis and launch SCF routine
    conv = convergence_parameters
    basis = PlaneWaveBasis(model; conv.Ecut, conv.kgrid, basis_kwargs...)
    self_consistent_field(basis; tol=conv.scf_tol)
end

function extract_local_KS_potential(scfres)
    basis = scfres.basis
    fft_size = basis.fft_size
    Nz = fft_size[3]

    # Apply substraction by far values
    V = DFTK.total_local_potential(scfres.ham)[:,:,:,1]

    # The potential V is not 0at infinity. We substract its value at infinity to correct
    # that. TODO: Is it a bug ? Maybe not.
    slab_average = sum(V; dims=(1,2))/prod(fft_size[1:2])  # Average of V over x and y
    Ns_mid = Nz/2; Ns1 = Ns_mid - Nz/10; Ns2 = Ns_mid + Nz/10
    Ns_mid, Ns1, Ns2 = floor.(Int, (Ns_mid, Ns1, Ns2))
    mid_slab = slab_average[Ns1:Ns2]
    m = sum(mid_slab)/length(mid_slab)

    V .-= m
    fft(basis, V)
end

function local_KS_potential(BG::BilayerGraphene;
                      terms=nothing,
                      convergence_parameters = default_convergence_parameters,
                      # Temperature
                      temperature=1e-3,
                      smearing=DFTK.Smearing.Gaussian(),
                      # To be able to ask for a specific fft grid
                      basis_kwargs...
                      )
    scfres = compute_scfres(BG; terms, convergence_parameters, temperature, smearing, basis_kwargs...)
    extract_local_KS_potential(scfres)
end
