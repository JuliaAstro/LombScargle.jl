### LombScargle.jl ---  Perform Lomb–Scargle periodogram
#
# Copyright (C) 2016 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle
#
# This file is a part of LombScargle.jl.
#
# License is MIT "Expat".
#
### Code:

__precompile__()

module LombScargle

using Measurements

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
                    w::AbstractVector{<:Real},
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
        YY = w ⋅ (signal .^ 2)
        return P .*= N .* YY ./ 2
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
    normalize!(P, p.signal, p.w, length(p.signal), p.noise, p.norm)

lombscargle(p::PeriodogramPlan) =
    Periodogram(normalize!(_periodogram!(p), p), p.freq, p.times, p.norm)

lombscargle(args...; kwargs...) = lombscargle(plan(args...; kwargs...))

"""
    lombscargle(times::AbstractVector{Real}, signal::AbstractVector{Real},
                errors::AbstractVector{Real}=ones(signal);
                normalization::Symbol=:standard,
                noise_level::Real=1,
                center_data::Bool=true,
                fit_mean::Bool=true,
                fast::Bool=true,
                oversampling::Integer=5,
                Mfft::Integer=4,
                samples_per_peak::Integer=5,
                nyquist_factor::Integer=5,
                minimum_frequency::Real=NaN,
                maximum_frequency::Real=NaN,
                frequencies::AbstractVector{Real}=
                autofrequency(times,
                              samples_per_peak=samples_per_peak,
                              nyquist_factor=nyquist_factor,
                              minimum_frequency=minimum_frequency,
                              maximum_frequency=maximum_frequency))

Compute the Lomb–Scargle periodogram of the `signal` vector, observed at
`times`.  You can also specify the uncertainties for each signal point with in
`errors` argument.  All these vectors must have the same length.

Optional keywords arguments are:

* `normalization`: how to normalize the periodogram.  Valid choices are:
  `:standard`, `:model`, `:log`, `:psd`, `:Scargle`, `:HorneBaliunas`,
  `:Cumming`
* `noise_level`: the noise level used to normalize the periodogram when
  `normalization` is set to `:Scargle`
* `fit_mean`: if `true`, fit for the mean of the signal using the Generalised
  Lomb–Scargle algorithm (see Zechmeister & Kürster paper below).  If this is
  `false` and no uncertainty on the signal is provided, the original algorithm
  by Lomb and Scargle will be employed (see Townsend paper below)
* `center_data`: if `true`, subtract the weighted mean of `signal` from `signal`
  itself before performing the periodogram.  This is especially important if
  `fit_mean` is `false`
* `frequencies`: the frequecy grid (not angular frequencies) at which the
  periodogram will be computed, as a vector.  If not provided, it is an evenly
  spaced grid of type `Range`, automatically determined with
  `LombScargle.autofrequency` function, which see.  See below for other
  available keywords that can be used to affect the frequency grid without
  directly setting `frequencies`
* `fast`: whether to use the fast method by Press & Rybicki, overriding the
  default choice.  In any case, `frequencies` must be a `Range` object in order
  to use this method (this is the default)
* `oversampling`: oversampling the frequency factor for the approximation;
  roughly the number of time samples across the highest-frequency sinusoid.
  This parameter contains the tradeoff between accuracy and speed.  Used only
  when the fast method is employed
* `Mfft`: the number of adjacent points to use in the FFT approximation.  Used
  only when the fast method is employed

In addition, you can use all optional keyword arguments of
`LombScargle.autofrequency` function in order to tune the `frequencies` vector
without calling the function:

* `samples_per_peak`: the approximate number of desired samples across the
  typical peak
* `nyquist_factor`: the multiple of the average Nyquist frequency used to choose
  the maximum frequency if `maximum_frequency` is not provided
* `minimum_frequency`: if specified, then use this minimum frequency rather than
  one chosen based on the size of the baseline
* `maximum_frequency`: if specified, then use this maximum frequency rather than
  one chosen based on the average Nyquist frequency

If the signal has uncertainties, the `signal` vector can also be a vector of
`Measurement` objects (from
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package), in
which case you don’t need to pass a separate `errors` vector for the
uncertainties of the signal.  See `Measurements.jl` manual at
http://measurementsjl.readthedocs.io/ for details on how to create a vector of
`Measurement` objects.

The algorithm used here are reported in the following papers:

* Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
  http://dx.doi.org/10.1086/167197, Bibcode:
  http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
* Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
  http://dx.doi.org/10.1088/0067-0049/191/2/247,
  Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
* Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
  http://dx.doi.org/10.1051/0004-6361:200811296,
  Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
"""
lombscargle(::AbstractVector{<:Real}, rest...)

"""
    lombscargle(plan::PeriodogramPlan)

Compute the Lomb–Scargle periodogram for the given `plan`.  This method has no other
argument.  See documentation of `LombScargle.plan` for how to plan a Lomb–Scargle
periodogram.
"""
lombscargle(::PeriodogramPlan)

end # module
