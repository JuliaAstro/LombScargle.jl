### LombScargle.jl ---  Perform Lomb–Scargle periodogram
#
# Copyright (C) 2016, 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle
#
# This file is a part of LombScargle.jl.
#
# License is BSD 3-clause "New" or "Revised".
#
### Code:

__precompile__()

module LombScargle

using FFTW, Measurements, LinearAlgebra, Statistics, SpecialFunctions

export lombscargle

# This is similar to Periodogram type of DSP.Periodograms module, but for
# unevenly spaced frequencies.
struct Periodogram{P<:AbstractFloat, F<:AbstractVector, T<:AbstractVector}
    power::Vector{P}
    freq::F
    # XXX: the `times' vector is only in the `M' function (see utils.jl), but
    # only maximum(times) and minimum(times) are used.  We could consider the
    # possibility to keep in this type only the extrema of t, instead of the
    # whole array.
    times::T
    norm::Symbol
end

abstract type PeriodogramPlan end

include("utils.jl")
include("bootstrap.jl")

include("townsend.jl")
include("gls.jl")
include("press-rybicki.jl")

include("planning.jl")

function normalize!(P::AbstractVector{<:AbstractFloat},
                    signal::AbstractVector{<:Real},
                    psdfactor::Real,
                    N::Integer,
                    noise_level::Real,
                    normalization::Symbol)
    if normalization == :standard
        return P
    elseif normalization == :model
        return P ./= (1 .- P)
    elseif normalization == :log
        return P .= -log.(1 .- P)
    elseif normalization == :psd
        return P .*= psdfactor ./ 2
    elseif normalization == :Scargle
        return P ./= noise_level
    elseif normalization == :HorneBaliunas
        return P .*= (N .- 1) ./ 2
    elseif normalization == :Cumming
        M = maximum(P)
        return P .*= (N .- 3) ./ (1 .- M) ./ 2
    else
        error("normalization \"", string(normalization), "\" not supported")
    end
end

normalize!(P::AbstractVector{<:AbstractFloat}, p::PeriodogramPlan) =
    normalize!(P, p.signal, p.YY * p.sumw, length(p.signal), p.noise, p.norm)

lombscargle(p::PeriodogramPlan) =
    Periodogram(normalize!(_periodogram!(p), p), p.freq, p.times, p.norm)

lombscargle(args...; kwargs...) = lombscargle(plan(args...; kwargs...))

"""
    lombscargle(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},
                [errors::AbstractVector{<:Real}]; keywords...)

Compute the Lomb–Scargle periodogram of the `signal` vector, observed at
`times`.  You can also specify the uncertainties for each signal point with
`errors` argument.  All these vectors must have the same length.

All optional keywords are described in the docstring of
[`LombScargle.plan`](@ref).

If the signal has uncertainties, the `signal` vector can also be a vector of
`Measurement` objects (from
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package), in
which case you don’t need to pass a separate `errors` vector for the
uncertainties of the signal.
"""
lombscargle(::AbstractVector{<:Real}, rest...)

"""
    lombscargle(plan::PeriodogramPlan)

Compute the Lomb–Scargle periodogram for the given `plan`.  This method has no
other arguments.  See documentation of [`LombScargle.plan`](@ref) for how to
plan a Lomb–Scargle periodogram.
"""
lombscargle(::PeriodogramPlan)

_periodogram!(p::PeriodogramPlan) = _periodogram!(p.P, p.times, p)

end # module
