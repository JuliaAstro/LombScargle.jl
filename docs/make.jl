using Documenter, LombScargle

makedocs(
    modules = [LombScargle],
    format = :html,
    sitename = "LombScargle",
)

deploydocs(
    repo = "github.com/giordano/LombScargle.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
