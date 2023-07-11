module GaIn

# TODO: replace with gain_jll binary
const GaIn_lib = "/home/lvidal/programs/GaIn/lib/libGaIn.so"

# struct Ccplx
#     r::Cdouble
#     i::Cdouble
# end

export overlap, kinetic, ionic, coulomb, laplacian

function overlap(ζ1::Cdouble, R1::Vector{Cdouble}, nx1::Int, ny1::Int, nz1::Int,
                 ζ2::Cdouble, R2::Vector{Cdouble}, nx2::Int, ny2::Int, nz2::Int)
    @ccall GaIn_lib.c_overlap_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
end

function kinetic(ζ1::Cdouble, R1::Vector{Cdouble}, nx1::Int, ny1::Int, nz1::Int,
                 ζ2::Cdouble, R2::Vector{Cdouble}, nx2::Int, ny2::Int, nz2::Int)
    res = @ccall GaIn_lib.c_laplacian_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
    -res/2
end

function ionic(ζ1::Cdouble, R1::Vector{Cdouble}, nx1::Int, ny1::Int, nz1::Int,
               ζ2::Cdouble, R2::Vector{Cdouble}, nx2::Int, ny2::Int, nz2::Int)
    res = @ccall GaIn_lib.c_laplacian_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
    -res
end


function coulomb(ζ1::Cdouble, R1::Vector{Cdouble}, nx1::Int, ny1::Int, nz1::Int,
                     ζ2::Cdouble, R2::Vector{Cdouble}, nx2::Int, ny2::Int, nz2::Int)
    @ccall GaIn_lib.c_coulomb_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                 Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                 Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
end


function laplacian(ζ1::Cdouble, R1::Vector{Cdouble}, nx1::Int, ny1::Int, nz1::Int,
                   ζ2::Cdouble, R2::Vector{Cdouble}, nx2::Int, ny2::Int, nz2::Int)
    @ccall GaIn_lib.c_laplacian_c_(Ref(ζ1)::Ptr{Cdouble}, R1::Ptr{Cdouble}, Ref(Cint(nx1))::Ptr{Cint}, Ref(Cint(ny1))::Ptr{Cint},
                                   Ref(Cint(nz1))::Ptr{Cint}, Ref(ζ2)::Ptr{Cdouble}, R2::Ptr{Cdouble}, Ref(Cint(nx2))::Ptr{Cint},
                                   Ref(Cint(ny2))::Ptr{Cint}, Ref(Cint(nz2))::Ptr{Cint})::Cdouble
end

end
