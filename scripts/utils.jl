@views av(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

function preprocess_flowline(filename, outname)
    data = readdlm(filename; skipstart=1)

    x, y, h, b, w = data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]

    l = [0; cumsum(sqrt.(diff(x) .^ 2 .+ diff(y) .^ 2))]
    A = diff(l) .* av(w)

    open(outname, "w") do io
        println(io, "length(m),elevation(m),bed(m),area(m2)")
        writedlm(io, [av(l) av(h) av(b) A], ',')
    end

    return
end

function read_flowline(filename)
    data = readdlm(filename, ','; skipstart=1)
    l, h, b, A = data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    return l, h, b, A
end
