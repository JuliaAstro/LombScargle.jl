### utils.jl
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
### Commentary
#
# This file contains some utilities for the LombScargle.jl package.
#
### Code:

"""
    power(p::Periodogram)

Return the power vector of Lomb-Scargle periodogram `p`.
"""
power(p::Periodogram) = p.power

"""
    freq(p::Periodogram)

Return the frequency vector of Lomb-Scargle periodogram `p`.
"""
freq(p::Periodogram) = p.freq

"""
    freqpower(p::Periodogram)

Return the 2-tuple `(freq(p), power(p))`, where `freq(p)` and `power(p)` are the
frequency vector and the power vector of Lomb-Scargle periodogram `p`
respectively.
"""
freqpower(p::Periodogram) = (freq(p), power(p))

"""
    findmaxpower(p::Periodogram)

Return the highest power of the periodogram `p`.
"""
findmaxpower(p::Periodogram) = maximum(power(p))

"""
    findmaxfreq(p::Periodogram, threshold::Real=findmaxpower(p))

Return the array of frequencies with the highest power in the periodogram `p`.
If a second argument `threshold` is provided, return the frequencies with power
larger than or equal to `threshold`.
"""
findmaxfreq(p::Periodogram, threshold::Real=findmaxpower(p)) =
    freq(p)[find(x -> x >= threshold, power(p))]

"""
    autofrequency(times::AbstractVector{Real};
                  samples_per_peak::Integer=5,
                  nyquist_factor::Integer=5,
                  minimum_frequency::Real=NaN,
                  maximum_frequency::Real=NaN)

Determine a suitable frequency grid for the given vector of `times`.

Optional keyword arguments are:

* `samples_per_peak`: the approximate number of desired samples across the
  typical peak
* `nyquist_factor`: the multiple of the average Nyquist frequency used to choose
  the maximum frequency if `maximum_frequency` is not provided
* `minimum_frequency`: if specified, then use this minimum frequency rather than
  one chosen based on the size of the baseline
* `maximum_frequency`: if specified, then use this maximum frequency rather than
  one chosen based on the average Nyquist frequency

This is based on prescription given at
https://jakevdp.github.io/blog/2015/06/13/lomb-scargle-in-python/ and uses the
same keywords names adopted in Astropy.
"""
function autofrequency{R<:Real}(times::AbstractVector{R};
                                samples_per_peak::Integer=5,
                                nyquist_factor::Integer=5,
                                minimum_frequency::Real=NaN,
                                maximum_frequency::Real=NaN)
    T = maximum(times) - minimum(times)
    δf = inv(samples_per_peak*T)
    f_min = isfinite(minimum_frequency) ? minimum_frequency : 0.5*δf
    if isfinite(maximum_frequency)
        return f_min:δf:maximum_frequency
    else
        return f_min:δf:0.5*nyquist_factor*length(times)/T
    end
end

"""
    prob(P::Periodogram, p_0::Real)

Return the probability that the periodogram power can exceed the value `p_{0}`.

Its inverse is the `probinv` function.
"""
function prob(P::Periodogram, p_0::Real)
    N = length(P.times)
    normalization = P.norm
    if normalization == "standard"
        return (1.0 - p_0)^((N - 3.0)*0.5)
    elseif normalization == "Scargle"
        return exp(-p_0)
    elseif normalization == "HorneBaliunas"
        return (1.0 - 2p_0/(N - 1))^((N - 3.0)*0.5)
    elseif normalization == "Cumming"
        return (1.0 + 2p_0/(N - 3.0))^((3.0 - N)*0.5)
    else
        error("normalization \"", normalization, "\" not supported")
    end
end

"""
    probinv(P::Periodogram, prob::Real)

Return the `p_0` value of the periodogram power whose probability is `prob`.

This is the inverse of `prob` function.
"""
function probinv(P::Periodogram, prob::Real)
    N = length(P.times)
    normalization = P.norm
    if normalization == "standard"
        return 1.0 - prob^(2.0/(N - 3.0))
    elseif normalization == "Scargle"
        return -log(prob)
    elseif normalization == "HorneBaliunas"
        return 0.5*(N - 1.0)*(1.0 - prob^(2.0/(N - 3.0)))
    elseif normalization == "Cumming"
        return 0.5*(N - 3.0)*(prob^(2.0/(3.0 - N)) - 1.0)
    else
        error("normalization \"", normalization, "\" not supported")
    end
end

"""
    M(P::Periodogram)

Estimates the number of independent frequencies in the periodogram `P`.
"""
function M(P::Periodogram)
    t = P.times
    f = P.freq
    return (maximum(t) - minimum(t))*(maximum(f) - minimum(f))
end

"""
    fap(P::Periodogram, p_0::Real)

Return the false-alarm probability for periodogram `P` and power value `p_0`.

Its inverse is the `fapinv` function.
"""
fap(P::Periodogram, p_0::Real) = 1.0 - (1.0 - prob(P, p_0))^M(P)

"""
    fapinv(P::Periodogram, fap::Real)

Return the `p_0` value of the periodogram power whose false-alarm probability is
`fap`.

This is the inverse of `fap` function.
"""
fapinv(P::Periodogram, fap::Real) = probinv(P, 1.0 - (1.0 - fap)^(inv(M(P))))
