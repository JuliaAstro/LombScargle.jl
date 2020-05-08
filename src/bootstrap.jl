### bootstrap.jl
#
# Copyright (C) 2016, 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle, bootstrapping
#
# This file is a part of LombScargle.jl.
#
# License is BSD 3-clause "New" or "Revised".
#
### Commentary
#
# This file contains facilities to perform bootstrapping and to calculate
# false-alarm probability and its inverse.
#
### Code:

using Random

struct Bootstrap{T<:AbstractFloat}
    p::Vector{T} # Vector of highest peaks
end

const _default_rng = if VERSION < v"1.3"
    () -> Random.GLOBAL_RNG
else
    Random.default_rng
end

bootstrap(rng::AbstractRNG, N::Integer, p::PeriodogramPlan) =
    Bootstrap(sort!([maximum(normalize!(_periodogram!(shuffle(rng, p.times), p), p)) for _ in 1:N], rev = true))

bootstrap(rng::AbstractRNG, N::Integer, t::AbstractVector{<:Real}, rest...; kwargs...) =
    bootstrap(rng, N, plan(t, rest...; kwargs...))

bootstrap(N::Integer, t::AbstractVector{<:Real}, rest...; kwargs...) =
    bootstrap(_default_rng(), N, plan(t, rest...; kwargs...))

"""
    LombScargle.bootstrap(N::Integer,
                          times::AbstractVector{Real},
                          signal::AbstractVector{Real},
                          errors::AbstractVector{Real}=ones(signal); ...)

Create `N` bootstrap samples, perform the Lomb–Scargle analysis on them, and
store all the highest peaks for each one in a `LombScargle.Bootstrap` object.
All the arguments after `N` are passed around to [`lombscargle`](@ref).
"""
bootstrap(::Integer, ::AbstractVector{<:Real})

"""
    LombScargle.bootstrap(N::Integer, plan::PeriodogramPlan)

Create `N` bootstrap samples, perform the Lomb–Scargle analysis on them for the given
`plan`, and store all the highest peaks for each one in a `LombScargle.Bootstrap` object.

See documentation of [`LombScargle.plan`](@ref) for how to plan a Lomb–Scargle periodogram.
"""
bootstrap(::Integer, ::PeriodogramPlan)

"""
    fap(b::Bootstrap, power::Real)

Return the false-alarm probability for `power` in the bootstrap sample `b`.

Its inverse is the [`fapinv`](@ref) function.
"""
fap(b::Bootstrap{<:AbstractFloat}, power::Real) =
    length(findall(x -> x >= power, b.p))/length(b.p)

"""
    fapinv(b::Bootstrap, prob::Real)

Return the power value whose false-alarm probability is `prob` in the bootstrap
sample `b`.

It returns `NaN` if the requested probability is too low and the power cannot be
determined with the bootstrap sample `b`.  In this case, you should enlarge your
bootstrap sample so that `N*fap` can be rounded to an integer larger than or
equal to 1.

This is the inverse of [`fap`](@ref) function.
"""
fapinv(b::Bootstrap{<:AbstractFloat}, prob::Real) =
    get(b.p, round(Int, length(b.p)*prob), NaN)
