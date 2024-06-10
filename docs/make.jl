using Documenter
using DocumenterCitations
using LombScargle

# # Generate all images
# include(joinpath(@__DIR__, "src", "images.jl"))
cp(joinpath(@__DIR__, "..", "perf", "benchmarks.png"),
   joinpath(@__DIR__, "src", "benchmarks.png"))

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style=:authoryear,
)

makedocs(
    modules = [LombScargle],
    sitename = "LombScargle",
    plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaAstro/LombScargle.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
