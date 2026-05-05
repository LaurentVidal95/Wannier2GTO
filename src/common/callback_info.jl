# TODO add number of MOs and AOs

function default_callback()
    function callback(info)
        if info.n_iter == 0
            header = ["Iter", "Compression error (H$(info.Wc.error_norm) norm)",
                      "number of Gaussians"]
            println("-"^60); @printf "%-7s %-30s %-20s \n" header...; println("-"^60)
            flush(stdout)
        else
            # Compute number of Gaussians
            n_GTOs = sum(length(Φ.coeffs) for Φ in info.Wc.basis_functions)
            println(@sprintf("%5i %30.12f %20i", info.n_iter, info.Wc.error, n_GTOs))
            flush(stdout)
        end
        nothing
    end
end
