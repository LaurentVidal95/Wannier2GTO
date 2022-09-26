function default_callback()
    function callback(info)
        if info.n_iter == 0
            header = ["Iter", "J", "||Residual||"]
            println("-"^40); @printf "%-10s %-10s %10s \n" header...; println("-"^40)
            flush(stdout)
        end
    end
    callback
end
