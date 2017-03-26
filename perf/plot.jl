#!/usr/bin/env julia

using Plots
pyplot()

julia1 = readdlm("julia_times-1.dat")
julia4 = readdlm("julia_times-4.dat")
python = readdlm("python_times.dat")

plot(xaxis = (:log,), yaxis = (:log,),
     xlabel = "Datapoints", ylabel = "Time (seconds)",
     size=(900, 600))
plot!(julia1[:,1], julia1[:,2], linewidth = 2, marker = (:auto,), lab = "LombScargle.jl - single thread")
plot!(julia4[:,1], julia4[:,2], linewidth = 2, marker = (:auto,), lab = "LombScargle.jl - 4 threads")
plot!(python[:,1], python[:,2], linewidth = 2, marker = (:auto,), lab = "AstroPy")
savefig("benchmarks.png")
