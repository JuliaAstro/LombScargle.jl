#!/usr/bin/env julia

using LombScargle, BenchmarkTools, FFTW

N  = round.(Int, exp10.(range(1, stop=5, length=15)))
for nthreads in (1, 4)
    FFTW.set_num_threads(nthreads)
    @info(string(nthreads) * " FFTW thread(s)")
    open(joinpath(@__DIR__, "julia_times-" * string(nthreads) * ".dat"), "w") do file
        for (j, n) in enumerate(N)
            println("Iteration ", j, " (", n, " datapoints)")
            t = collect(range(0.01, stop=10, length=n))
            s = sin.(t) .+ 1.5*cospi.(4t) .+ 3
            plan = LombScargle.plan(t, s, fast = true, flags = FFTW.MEASURE)
            time = minimum(@benchmark(lombscargle($plan)).times)/1e9
            println(" Julia:  ", time, " seconds")
            println(file, n, "\t", time)
        end
    end
end
