function default_callback()
    function callback(info)
        if info.n_iter == 0
            header = ["Iter", "Residual H$(info.Wc.error_norm) norm"]
            println("-"^40); @printf "%5s %-16s \n" header...; println("-"^40)
            flush(stdout)
        else
            println(@sprintf("%5i %16.12f", info.n_iter, info.Wc.error))
            flush(stdout)
        end
        nothing
    end
end
