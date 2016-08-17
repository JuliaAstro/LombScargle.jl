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

using Measurements

export lombscargle, power, freq, freqpower, findmaxpower, findmaxfreq,
prob, probinv, fap, fapinv

include("extirpolation.jl")

# This is similar to Periodogram type of DSP.Periodograms module, but for
# unevenly spaced frequencies.
immutable Periodogram{T<:AbstractFloat}
    power::AbstractVector{T}
    freq::AbstractVector{T}
    times::AbstractVector{T}
    norm::AbstractString
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

# Original algorithm that doesn't take into account uncertainties and doesn't
# fit the mean of the signal.  This is implemented following the recipe by
# * Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
#   http://dx.doi.org/10.1088/0067-0049/191/2/247,
#   Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
function _lombscargle_orig{R1<:Real,R2<:Real,R3<:Real}(times::AbstractVector{R1},
                                                       signal::AbstractVector{R2},
                                                       freqs::AbstractVector{R3},
                                                       center_data::Bool)
    P_type = promote_type(float(R1), float(R2))
    P = Vector{P_type}(freqs)
    N = length(signal)
    nil = zero(P_type)
    # If "center_data" keyword is true, subtract the mean from each point.
    if center_data
        X = signal - mean(signal)
    else
        X = signal
    end
    XX = X⋅X

    @inbounds for n in eachindex(freqs)
        ω = 2pi*freqs[n]
        XC = XS = CC = CS = nil
        for j in eachindex(times)
            ωt = ω*times[j]
            C = cos(ωt)
            S = sin(ωt)
            XC += X[j]*C
            XS += X[j]*S
            CC += C*C
            CS += C*S
        end
        SS      = N - CC
        τ       = 0.5*atan2(2CS, CC - SS)/ω
        c_τ     = cos(ω*τ)
        s_τ     = sin(ω*τ)
        c_τ2    = c_τ*c_τ
        s_τ2    = s_τ*s_τ
        cs_τ_CS = 2c_τ*s_τ*CS
        P[n] = (abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))/XX
    end
    return P, freqs
end

# Generalised Lomb-Scargle algorithm: this takes into account uncertainties and
# fit the mean of the signal.  This is implemented following the recipe by
# * Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
#   http://dx.doi.org/10.1051/0004-6361:200811296,
#   Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
# In addition, some tricks suggested by
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# to make computation faster are adopted.
function _generalised_lombscargle{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                                       signal::AbstractVector{R2},
                                                                       w::AbstractVector{R3},
                                                                       freqs::AbstractVector{R4},
                                                                       center_data::Bool, fit_mean::Bool)
    P_type = promote_type(float(R1), float(R2))
    P = Vector{P_type}(freqs)
    nil = zero(P_type)
    # If "center_data" keyword is true, subtract the mean from each point.
    if center_data
        y = signal - mean(signal)
    else
        y = signal
    end
    YY = w⋅(y.^2)
    if fit_mean
        Y  = w⋅y
    end

    @inbounds for n in eachindex(freqs)
        ω = 2pi*freqs[n]

        # Find τ for current angular frequency
        C = S = CS = CC = nil
        for i in eachindex(times)
            ωt  = ω*times[i]
            W   = w[i]
            c   = cos(ωt)
            s   = sin(ωt)
            CS += W*c*s
            CC += W*c*c
            if fit_mean
                C  += W*c
                S  += W*s
            end
        end
        if fit_mean
            CS -= C*S
            SS  = 1.0 - CC - S*S
            CC -= C*C
        else
            SS  = 1.0 - CC
        end
        τ       = 0.5*atan2(2CS, CC - SS)/ω
        # These quantities will be used below.  See Press & Rybicki paper.
        # NOTE: They can be computed in a faster way, for the time being we keep
        # this simple and straightforward implementation.
        cos_ωτ  = cos(ω*τ)
        sin_ωτ  = sin(ω*τ)

        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = nil
        for i in eachindex(times)
            ωt    = ω*(times[i] - τ)
            W     = w[i]
            c     = cos(ωt)
            s     = sin(ωt)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
        end
        if fit_mean
            # "C_τ" and "S_τ" are computed following equation (7) of Press &
            # Rybicki, this formula simply comes from angle difference
            # trigonometric identities.
            C_τ   = C*cos_ωτ + S*sin_ωτ
            S_τ   = S*cos_ωτ - C*sin_ωτ
            YC_τ -= Y*C_τ
            YS_τ -= Y*S_τ
            SS_τ  = 1.0 - CC_τ - S_τ*S_τ
            CC_τ -= C_τ*C_τ
            YY   -= Y*Y
        else
            SS_τ  = 1.0 - CC_τ
        end
        P[n] = (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY
    end
    return P, freqs
end

# This is the switch to select the appropriate function to run
function _lombscargle{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                           signal::AbstractVector{R2},
                                                           with_errors::Bool,
                                                           w::AbstractVector{R3}=ones(signal)/length(signal);
                                                           normalization::AbstractString="standard",
                                                           noise_level::Real=1.0,
                                                           center_data::Bool=true, fit_mean::Bool=true,
                                                           samples_per_peak::Integer=5,
                                                           nyquist_factor::Integer=5,
                                                           minimum_frequency::Real=NaN,
                                                           maximum_frequency::Real=NaN,
                                                           frequencies::AbstractVector{R4}=
                                                           autofrequency(times,
                                                                         samples_per_peak=samples_per_peak,
                                                                         nyquist_factor=nyquist_factor,
                                                                         minimum_frequency=minimum_frequency,
                                                                         maximum_frequency=maximum_frequency))
    @assert length(times) == length(signal) == length(w)
    if fit_mean || with_errors
        P, f = _generalised_lombscargle(times, signal, w, frequencies,
                                        center_data, fit_mean)
    else
        P, f = _lombscargle_orig(times, signal, frequencies, center_data)
    end

    N = length(signal)
    # Normalize periodogram
    if normalization == "standard"
    elseif normalization == "model"
        P = P./(1.0 - P)
    elseif normalization == "log"
        P = -log(1.0 - P)
    elseif normalization == "psd"
        P *= 0.5*N*(w⋅signal.^2)
    elseif normalization == "Scargle"
        P /= noise_level
    elseif normalization == "HorneBaliunas"
        P *= 0.5*(N - 1.0)
    elseif normalization == "Cumming"
        P *= 0.5*(N - 3.0)/(1.0 - maximum(P))
    else
        error("normalization \"", normalization, "\" not supported")
    end

    return Periodogram(P, f, times, normalization)
end

# No uncertainties
lombscargle{R1<:Real,R2<:Real}(times::AbstractVector{R1},
                               signal::AbstractVector{R2};
                               kwargs...) =
                                   _lombscargle(times,
                                                signal,
                                                false;
                                                kwargs...)

# Uncertainties provided
function lombscargle{R1<:Real,R2<:Real,R3<:Real}(times::AbstractVector{R1},
                                                 signal::AbstractVector{R2},
                                                 errors::AbstractVector{R3};
                                                 kwargs...)
    # Compute weights vector
    w = 1.0./errors.^2
    w /= sum(w)
    return _lombscargle(times, signal, true, w; kwargs...)
end

# Uncertainties provided via Measurement type
lombscargle{T<:Real,F<:AbstractFloat}(times::AbstractVector{T},
                                      signal::AbstractVector{Measurement{F}};
                                      kwargs...) =
                                          lombscargle(times,
                                                      value(signal),
                                                      uncertainty(signal);
                                                      kwargs...)

"""
    lombscargle(times::AbstractVector{Real}, signal::AbstractVector{Real},
                errors::AbstractVector{Real}=ones(signal);
                normalization::AbstractString="standard",
                noise_level::Real=1.0,
                center_data::Bool=true, fit_mean::Bool=true,
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

Compute the Lomb-Scargle periodogram of the `signal` vector, observed at
`times`.  You can also specify the uncertainties for each signal point with in
`errors` argument.  All these vectors must have the same length.

Optional keywords arguments are:

* `normalization`: how to normalize the periodogram.  Valid choices are:
  `"standard"`, `"model"`, `"log"`, `"psd"`, `"Scargle"`, `"HorneBaliunas"`,
  `"Cumming"`
* `noise_level`: the noise level used to normalize the periodogram when
  `normalization` is set to `"Scargle"`
* `fit_mean`: if `true`, fit for the mean of the signal using the Generalised
  Lomb-Scargle algorithm (see Zechmeister & Kürster paper below).  If this is
  `false` and no uncertainty on the signal is provided, the original algorithm
  by Lomb and Scargle will be employed (see Townsend paper below)
* `center_data`: if `true`, subtract the mean of `signal` from `signal` itself
  before performing the periodogram.  This is especially important if `fit_mean`
  is `false`
* `frequencies`: the frequecy grid (not angular frequencies) at which the
  periodogram will be computed, as a vector.  If not provided, it is
  automatically determined with `LombScargle.autofrequency` function, which see.
  See below for other available keywords that can be used to affect the
  frequency grid without directly setting `frequencies`

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

* Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
  http://dx.doi.org/10.1088/0067-0049/191/2/247,
  Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
* Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
  http://dx.doi.org/10.1051/0004-6361:200811296,
  Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
"""
lombscargle

end # module
