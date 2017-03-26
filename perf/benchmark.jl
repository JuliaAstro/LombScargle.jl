#!/usr/bin/env julia

using LombScargle, BenchmarkTools

N  = round.(Int, logspace(1, 5, 15))
for nthreads in (1, 4)
    FFTW.set_num_threads(nthreads)
    info(string(nthreads) * " FFTW thread(s)")
    open(joinpath(@__DIR__, "julia_times-" * string(nthreads) * ".dat"), "w") do file
        for (j, n) in enumerate(N)
            println("Iteration ", j, " (", n, " datapoints)")
            t = collect(linspace(0.01, 10, n))
            s = sin.(t) .+ 1.5*cospi.(4t) .+ 3
            plan = LombScargle.plan(t, s, fast = true)
            time = minimum(@benchmark(lombscargle($plan)).times)/1e9
            println(" Julia:  ", time, " seconds")
            println(file, n, "\t", time)
        end
    end
end
