using Documenter, LombScargle

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
