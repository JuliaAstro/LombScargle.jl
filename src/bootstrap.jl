### bootstrap.jl
#
# Copyright (C) 2016, 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle, bootstrapping
#
# This file is a part of LombScargle.jl.
#
# License is MIT "Expat".
#
### Commentary
#
# This file contains facilities to perform bootstrapping and to calculate
# false-alarm probability and its inverse.
#
### Code:

immutable Bootstrap{T<:AbstractFloat}
    p::Vector{T} # Vector of highest peaks
end

bootstrap(N::Integer, p::PeriodogramPlan) =
    Bootstrap(sort!([maximum(normalize!(_periodogram!(shuffle(p.times), p), p)) for _ in 1:N], rev = true))

bootstrap(N::Integer, t::AbstractVector{<:Real}, rest...; kwargs...) =
    bootstrap(N, plan(t, rest...; kwargs...))

"""
    LombScargle.bootstrap(N::Integer,
                          times::AbstractVector{Real},
                          signal::AbstractVector{Real},
                          errors::AbstractVector{Real}=ones(signal); ...)

Create `N` bootstrap samples, perform the Lomb–Scargle analysis on them, and
store all the highest peaks for each one in a `LombScargle.Bootstrap` object.
All the arguments after `N` are passed around to `lombscargle`, which see.
"""
bootstrap(::Integer, ::AbstractVector{<:Real})

"""
    LombScargle.bootstrap(N::Integer, plan::PeriodogramPlan)

Create `N` bootstrap samples, perform the Lomb–Scargle analysis on them for the given
`plan`, and store all the highest peaks for each one in a `LombScargle.Bootstrap` object.

See documentation of `LombScargle.plan` for how to plan a Lomb–Scargle periodogram.
"""
bootstrap(::Integer, ::PeriodogramPlan)

"""
    fap(b::Bootstrap, power::Real)

Return the false-alarm probability for `power` in the bootstrap sample `b`.

Its inverse is the `fapinv` function.
"""
fap{T<:AbstractFloat}(b::Bootstrap{T}, power::Real) =
    length(find(x -> x >= power, b.p))/length(b.p)

"""
    fapinv(b::Bootstrap, prob::Real)

Return the power value whose false-alarm probability is `prob` in the bootstrap
sample `b`.

It returns `NaN` if the requested probability is too low and the power cannot be
determined with the bootstrap sample `b`.  In this case, you should enlarge your
bootstrap sample so that `N*fap` can be rounded to an integer larger than or
equal to 1.

This is the inverse of `fap` function.
"""
fapinv{T<:AbstractFloat}(b::Bootstrap{T}, prob::Real) =
    get(b.p, round(Int, length(b.p)*prob), NaN)
