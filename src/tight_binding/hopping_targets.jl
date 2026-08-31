using Random

@doc raw"""
Hopping-loss targets (design note 12 §4): labeled displacement set split into
training/validation, with reference values ``S_{\mathrm{ref}}(R),
T_{\mathrm{ref}}(R)`` (norm-corrected) and the on-site scale ``T(0)``.
`ortho_idx` indexes the training intralayer entries whose exact reference is
``S = 0`` (MLWF orthonormality) — they receive the extra ν penalty.
All lengths in Bohr, all Rs Cartesian.
"""
struct HoppingTargets{T<:Real}
    labels    :: Vector{String}
    sets      :: Vector{Symbol}       # :training | :validation
    Rs        :: Vector{Vector{T}}
    S_ref     :: Vector{T}
    T_ref     :: Vector{T}
    T0_ref    :: T
    ortho_idx :: Vector{Int}
end

@doc raw"""
Displacement sets of design note 12 §4. Training: note-09 set minus 2a₁
(periodic-image contamination). Validation: deterministic core (δ=a₁/4 shift,
two alternative z distances, a₂ for the C₃ check) plus `n_random` seeded
random draws in the overlap shell `r_min ≤ |R| ≤ r_max`.
`d_inter` and `z_val` are z distances in Bohr; the caller converts from Å.
"""
function hopping_R_sets(a₁::AbstractVector, a₂::AbstractVector;
                        d_inter::Real, z_val::Tuple, seed::Int,
                        n_random::Int, r_min::Real, r_max::Real)
    @assert 0 < r_min < r_max "invalid shell: [$r_min, $r_max]"
    ẑ(d) = [0.0, 0.0, float(d)]
    entries = Tuple{String, Vector{Float64}, Symbol}[
        ("intra a1",        Vector{Float64}(a₁),                :training),
        ("intra a1+a2",     Vector{Float64}(a₁ + a₂),           :training),
        ("inter AA",        ẑ(d_inter),                          :training),
        ("inter AB-like",   ẑ(d_inter) + (a₁ + a₂) / 3,          :training),
        ("inter mid",       ẑ(d_inter) + a₁ / 2,                 :training),
        ("val inter d=a1/4", ẑ(d_inter) + a₁ / 4,                :validation),
        ("val inter z1",     ẑ(z_val[1]),                        :validation),
        ("val inter z2",     ẑ(z_val[2]),                        :validation),
        ("val intra a2",     Vector{Float64}(a₂),                :validation),
    ]
    rng = MersenneTwister(seed)
    for k in 1:n_random
        u = randn(rng, 3); u /= norm(u)
        r = r_min + (r_max - r_min) * rand(rng)
        push!(entries, ("rand $k", r * u, :validation))
    end
    labels = [e[1] for e in entries]
    Rs     = [e[2] for e in entries]
    sets   = [e[3] for e in entries]
    ortho_idx = findall(l -> startswith(l, "intra"), labels)
    (; labels, Rs, sets, ortho_idx)
end

function build_hopping_targets(w_fourier::AbstractVector,
                               basis_supercell::PlaneWaveBasis, rsets)
    refs = [reference_hopping(w_fourier, basis_supercell, R) for R in rsets.Rs]
    T0 = reference_hopping(w_fourier, basis_supercell, zeros(3)).T0_ref
    HoppingTargets(rsets.labels, rsets.sets, rsets.Rs,
                   [r.S_ref for r in refs], [r.T_ref for r in refs],
                   T0, rsets.ortho_idx)
end

function store(targets::HoppingTargets; file::String)
    data = Dict("labels" => targets.labels,
                "sets" => String.(targets.sets),
                "Rs" => targets.Rs,
                "S_ref" => targets.S_ref, "T_ref" => targets.T_ref,
                "T0_ref" => targets.T0_ref, "ortho_idx" => targets.ortho_idx)
    open(io -> JSON3.write(io, data), file, "w")
    nothing
end

function HoppingTargets(file::String)
    d = open(JSON3.read, file)
    HoppingTargets(String.(d.labels), Symbol.(d.sets),
                   [Vector{Float64}(R) for R in d.Rs],
                   Vector{Float64}(d.S_ref), Vector{Float64}(d.T_ref),
                   Float64(d.T0_ref), Vector{Int}(d.ortho_idx))
end
