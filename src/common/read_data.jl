# Extract stored data
function read_wannier_dir(dir)
    readdlm_complex(file) = readdlm(file, Float64) |> x -> Complex.(x[:,1], x[:,2])
    wn = [readdlm_complex(joinpath(dir,"wn_$(k).dat")) for k in 1:64];
    αn = vec(readdlm(joinpath(dir,"wn_center.dat")))
    θn = only(readdlm(joinpath(dir,"wn_blob_angle.dat")))
    wn, αn, θn
end
