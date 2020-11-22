using Documenter, LombScargle

# # Generate all images
# include(joinpath(@__DIR__, "src", "images.jl"))
cp(joinpath(@__DIR__, "..", "perf", "benchmarks.png"),
   joinpath(@__DIR__, "src", "benchmarks.png"))

makedocs(
    modules = [LombScargle],
    sitename = "LombScargle",
)

deploydocs(
    repo = "github.com/JuliaAstro/LombScargle.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
