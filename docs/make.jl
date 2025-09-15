using Documenter
using LombScargle

# # Generate all images
# include(joinpath(@__DIR__, "src", "images.jl"))
cp(joinpath(@__DIR__, "..", "perf", "benchmarks.png"),
   joinpath(@__DIR__, "src", "benchmarks.png"))

# gives `pages` and `bib`
include("pages.jl")

makedocs(;
    modules = [LombScargle],
    sitename = "LombScargle",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.org/LombScargle/stable/",
    ),
    pages = pages,
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/JuliaAstro/LombScargle.jl.git",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
