module GaIn

# TODO: replace with gain_jll binary
const GaIn_lib = "/home/lvidal/programs/GaIn/lib/libGaIn.so"

# struct Ccplx
#     r::Cdouble
#     i::Cdouble
# end

export overlap, kinetic, coulomb, laplacian, atomic

function overlap(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                 ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int) where {T1, T2 <:Real}
    @ccall GaIn_lib.c_overlap_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
end

function overlap_upper_bound(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                             ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int) where {T1, T2 <:Real}
    l1 = max(nx1, ny1, nz1)
    l2 = max(nx2, ny2, nz2)
    @ccall GaIn_lib.overlap_upper_bound_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(l1))::Ptr{Cint},
                                        Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(l2))::Ptr{Cint},
                                        )::Cdouble
end

function coulomb(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                 ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int) where {T1, T2 <:Real}
    @ccall GaIn_lib.c_coulomb_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
end


function laplacian(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                   ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int) where {T1, T2 <:Real}
    res = @ccall GaIn_lib.c_laplacian_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                   Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                         Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
    -res
end

function kinetic(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                 ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int) where {T1, T2 <:Real}
    (1/2) * laplacian(ζ1, R1, nx1, ny1, nz1, ζ2, R2, nx2, ny2, nz2)
end

function atomic(ζ1::T1, R1::AbstractVector{T2}, nx1::Int, ny1::Int, nz1::Int,
                ζ2::T1, R2::AbstractVector{T2}, nx2::Int, ny2::Int, nz2::Int,
                R_at::AbstractVector{T2}) where {T1, T2 <:Real}
    @ccall GaIn_lib.cc_coulomb_ion_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                    Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                    Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint}, R_at::Ptr{Cdouble})::Cdouble
end

end
