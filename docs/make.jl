using Documenter, LombScargle

makedocs(
    modules = [LombScargle],
    format = :html,
    sitename = "LombScargle",
)

deploydocs(
    repo = "github.com/JuliaAstro/LombScargle.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
