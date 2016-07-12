### LombScargle.jl ---  Perform Lomb-Scargle periodogram
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

export lombscargle, power, freq, freqpower, findmaxfreq, autofrequency

# This is similar to Periodogram type of DSP.Periodograms module, but for
# unevenly spaced frequencies.
immutable Periodogram{T<:AbstractFloat}
    power::AbstractVector{T}
    freq::AbstractVector{T}
end

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
    findmaxfreq(p::Periodogram)
    findmaxfreq(p::Periodogram, threshold::Real)

Return the frequencies with the highest power in the periodogram `p`.  If a
second argument `threshold` is provided, return the frequencies with power
larger than or equal to `threshold`.
"""
findmaxfreq(p::Periodogram, threshold::Real=maximum(power(p))) =
    freq(p)[find(x -> x >= threshold, power(p))]

"""
    autofrequency(times::AbstractVector{Real};
                  samples_per_peak::Int=5,
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
function autofrequency{R<:Real}(times::AbstractVector{R}; samples_per_peak::Int=5,
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

# Original algorithm that doesn't take into account uncertainties and doesn't
# fit the mean of the signal.  This is implemented following the recipe by
# * Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
#   http://dx.doi.org/10.1088/0067-0049/191/2/247,
#   Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
function _lombscargle_orig{T<:Real}(times::AbstractVector{T}, signal::AbstractVector{T},
                                    freqs::AbstractVector{T}, center_data::Bool)
    P = Vector{T}(freqs)
    # If "center_data" keyword is true, subtract the mean from each point.
    signal_mean = center_data ? mean(signal) : zero(T)
    @inbounds for n in eachindex(freqs)
        ω = 2pi*freqs[n]
        XX = XC = XS = CC = SS = CS = zero(float(T))
        for j in eachindex(times)
            ωt = ω*times[j]
            C = cos(ωt)
            S = sin(ωt)
            X = signal[j] - signal_mean
            XX += X*X
            XC += X*C
            XS += X*S
            CC += C*C
            SS += S*S
            CS += C*S
        end
        τ       = 0.5*atan2(2CS, CC - SS)/ω
        c_τ     = cos(ω*τ)
        s_τ     = sin(ω*τ)
        c_τ2    = c_τ*c_τ
        s_τ2    = s_τ*s_τ
        cs_τ_CS = 2c_τ*s_τ*CS
        P[n] = (abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))/XX
    end
    return Periodogram(P, freqs)
end

# Generalised Lomb-Scargle algorithm: this takes into account uncertainties and
# fit the mean of the signal.  This is implemented following the recipe by
# * Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
#   http://dx.doi.org/10.1051/0004-6361:200811296,
#   Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
function _generalised_lombscargle{T<:Real}(times::AbstractVector{T},
                                           signal::AbstractVector{T},
                                           freqs::AbstractVector{T},
                                           errors::AbstractVector{T},
                                           center_data::Bool, fit_mean::Bool)
    P = Vector{T}(freqs)
    nil = zero(float(T))
    # If "center_data" keyword is true, subtract the mean from each point.
    signal_mean = center_data ? mean(signal) : nil
    w = 1.0./errors.^2
    w /= sum(w)
    @inbounds for n in eachindex(freqs)
        ω = 2pi*freqs[n]

        # Find τ for current angular frequency
        C = S = CS = CC = SS = nil
        for i in eachindex(times)
            ωt  = ω*times[i]
            W   = w[i]
            c   = cos(ωt)
            s   = sin(ωt)
            C  += W*c
            S  += W*s
            CS += W*c*s
            CC += W*c*c
            SS += W*s*s
        end
        CS -= C*S
        CC -= C*C
        SS -= S*S
        τ   = 0.5*atan2(2CS, CC - SS)/ω

        # Now we can compute the power
        Y = YY = C_τ = S_τ = YC_τ = YS_τ = CC_τ = SS_τ = nil
        for i in eachindex(times)
            y     = signal[i] - signal_mean
            ωt    = ω*(times[i] - τ)
            W     = w[i]
            c     = cos(ωt)
            s     = sin(ωt)
            Y    += W*y
            YY   += W*y*y
            C_τ  += W*c
            S_τ  += W*s
            YC_τ += W*y*c
            YS_τ += W*y*s
            CC_τ += W*c*c
            SS_τ += W*s*s
        end
        P[n] = (abs2(YC_τ - Y*C_τ)/(CC_τ - C_τ*C_τ) +
                abs2(YS_τ - Y*S_τ)/(SS_τ - S_τ*S_τ))/(YY - Y*Y)
    end
    return Periodogram(P, freqs)
end

"""
    lombscargle(times::AbstractVector{Real}, signal::AbstractVector{Real},
                freqs::AbstractVector{Real}=autofrequency(times),
                errors::AbstractVector{Real}=ones(signal);
                center_data::Bool=true, fit_mean::Bool=true)

Compute the Lomb-Scargle periodogram of the `signal` vector, observed at
`times`.  The frequecy grid `freqs` at which the periodogram will be computed is
automatically determined with `autofrequency` function, which see, if not
specified.  You can also specify the uncertainties for each signal point with in
`errors` argument.

Optional keywords arguments are:

* `fit_mean`: if `true`, fit for the mean of the signal using the Generalised
  Lomb-Scargle algorithm.  If this is `false`, the original algorithm by Lomb
  and Scargle will be employed, which does not take into account a non-null mean
  and uncertainties for observations
* `center_data`: if `true`, subtract the mean of `signal` from `signal` itself
  before performing the periodogram.  This is especially important if `fit_mean`
  is `false`

The algorithm used here are reported in the following papers:

* Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
  http://dx.doi.org/10.1088/0067-0049/191/2/247,
  Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
* Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
  http://dx.doi.org/10.1051/0004-6361:200811296,
  Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
"""
function lombscargle{T<:Real}(times::AbstractVector{T}, signal::AbstractVector{T},
                              freqs::AbstractVector{T}=autofrequency(times),
                              errors::AbstractVector{T}=ones(signal);
                              center_data::Bool=true, fit_mean::Bool=true)
    @assert length(times) == length(signal) == length(errors)
    if fit_mean
        return _generalised_lombscargle(times, signal, freqs, errors,
                                        center_data, fit_mean)
    else
        return _lombscargle_orig(times, signal, freqs, center_data)
    end
end

end # module
