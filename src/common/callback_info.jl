function default_callback()
    function callback(info)
        if info.n_iter == 0
            header = ["Iter", "J", "Residual H$(info.s) norm"]
            println("-"^40); @printf "%5s %-16s %16s \n" header...; println("-"^40)
            flush(stdout)
        else
            output = [info.n_iter, info.J, info.res_norm]
            println(@sprintf("%5i %16.12f %16.12f", output...))
            flush(stdout)
        end
        nothing
    end
end
