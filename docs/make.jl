using Documenter
using LombScargle

# # Generate all images
# include(joinpath(@__DIR__, "src", "images.jl"))
cp(joinpath(@__DIR__, "..", "perf", "benchmarks.png"),
   joinpath(@__DIR__, "src", "benchmarks.png"))

# gives `pages` and `bib`
include("pages.jl")

makedocs(
    modules = [LombScargle],
    sitename = "LombScargle",
    pages = pages,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaAstro/LombScargle.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
