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

export lombscargle

# This is similar to Periodogram type of DSP.Periodograms module, but for
# unevenly spaced frequencies.
immutable Periodogram{T<:AbstractFloat}
    power::AbstractVector{T}
    freq::AbstractVector{T}
    times::AbstractVector{T}
    norm::AbstractString
end

include("extirpolation.jl")
include("utils.jl")

# Original algorithm that doesn't take into account uncertainties and doesn't
# fit the mean of the signal.  This is implemented following the recipe by
# * Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
#   http://dx.doi.org/10.1088/0067-0049/191/2/247,
#   Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
function _lombscargle_orig{R1<:Real,R2<:Real,R3<:Real}(times::AbstractVector{R1},
                                                       signal::AbstractVector{R2},
                                                       freqs::AbstractVector{R3},
                                                       center_data::Bool)
    P_type = promote_type(float(R1), float(R2), float(R3))
    P = Vector{P_type}(length(freqs))
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
        ωτ      = 0.5*atan2(2CS, CC - SS)
        c_τ     = cos(ωτ)
        s_τ     = sin(ωτ)
        c_τ2    = c_τ*c_τ
        s_τ2    = s_τ*s_τ
        cs_τ_CS = 2c_τ*s_τ*CS
        P[n] = (abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))/XX
    end
    return P
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
                                                                       center_data::Bool,
                                                                       fit_mean::Bool)
    P_type = promote_type(float(R1), float(R2), float(R3), float(R4))
    P = Vector{P_type}(length(freqs))
    nil = zero(P_type)
    # If "center_data" keyword is true, subtract the mean from each point.
    if center_data || fit_mean
        y = signal - mean(signal)
    else
        y = signal
    end
    YY = w⋅(y.^2)
    if fit_mean
        Y   = w⋅y
        YY -= Y*Y
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
        ωτ   = 0.5*atan2(2CS, CC - SS)

        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = SS_τ = nil
        for i in eachindex(times)
            ωt    = ω*times[i] - ωτ
            W     = w[i]
            c     = cos(ωt)
            s     = sin(ωt)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
            SS_τ += W*s*s
        end
        if fit_mean
            # "C_τ" and "S_τ" are computed following equation (7) of Press &
            # Rybicki, this formula simply comes from angle difference
            # trigonometric identities.
            cos_ωτ = cos(ωτ)
            sin_ωτ = sin(ωτ)
            C_τ    = C*cos_ωτ + S*sin_ωτ
            S_τ    = S*cos_ωτ - C*sin_ωτ
            YC_τ  -= Y*C_τ
            YS_τ  -= Y*S_τ
            CC_τ  -= C_τ*C_τ
            # Note: it is possible to calculate SS_τ non-iteratively using the
            # formula "1.0 - CC_τ - S_τ*S_τ" (or "1.0 - CC_τ" if `fit_mean' is
            # false), but this seems to lead to numerical issues, like SS_τ ==
            # 0.0 and so P[n] == Inf.
            SS_τ  -= S_τ*S_τ
        end
        P[n] = (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY
    end
    return P
end

# Fast, but approximate, method to compute the Lomb-Scargle periodogram for
# evenly spaced data.  See
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# This is adapted from Astropy implementation of the method.  See
# `lombscargle_fast' function.
function _press_rybicki{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::Range{R1},
                                                             signal::AbstractVector{R2},
                                                             w::AbstractVector{R3},
                                                             freqs::Range{R4},
                                                             center_data::Bool,
                                                             fit_mean::Bool,
                                                             oversampling::Int,
                                                             Mfft::Int)
    P_type = promote_type(float(R1), float(R2))
    P = Vector{P_type}(length(freqs))
    nil = zero(P_type)
    # If "center_data" keyword is true, subtract the mean from each point.
    if center_data || fit_mean
        y = signal - w⋅signal
    else
        y = signal
    end

    df = step(freqs)
    N  = length(freqs)
    f0 = minimum(freqs)
    #---------------------------------------------------------------------------
    # 1. compute functions of the time-shift tau at each frequency
    Ch, Sh = trig_sum(times, w .* y, df, N, f0, 1, oversampling, Mfft)
    C2, S2 = trig_sum(times, w,      df, N, f0, 2, oversampling, Mfft)
    if fit_mean
        C, S = trig_sum(times, w, df, N, f0, 1, oversampling, Mfft)
        tan_2ωτ = (S2 .- 2.0 * S .* C) ./ (C2 .- (C .* C .- S .* S))
    else
        tan_2ωτ = S2 ./ C2
    end
    # This is what we're computing below; the straightforward way is slower and
    # less stable, so we use trig identities instead
    #
    # ωτ = 0.5 * atan(tan_2ωτ)
    # S2w, C2w = sin(2 * ωτ), cos(2 * ωτ)
    # Sw, Cw = sin(ωτ), cos(ωτ)
    C2w = 1./(sqrt(1.0 + tan_2ωτ .* tan_2ωτ))
    S2w = tan_2ωτ .* C2w
    Cw = sqrt(0.5 * (1.0 + C2w))
    Sw = sqrt(0.5) .* sign(S2w) .* sqrt(1.0 - C2w)
    #---------------------------------------------------------------------------
    # 2. Compute the periodogram, following Zechmeister & Kurster
    #    and using tricks from Press & Rybicki.
    YY = w⋅(y.^2)
    YC = Ch .* Cw .+ Sh .* Sw
    YS = Sh .* Cw .- Ch .* Sw
    CC = 0.5 * (1.0 + C2 .* C2w .+ S2 .* S2w)
    SS = 0.5 * (1.0 - C2 .* C2w .- S2 .* S2w)
    if fit_mean
        CC -= (C .* Cw .+ S .* Sw) .^ 2
        SS -= (S .* Cw .- C .* Sw) .^ 2
    end
    P = (YC .* YC ./ CC .+ YS .* YS ./ SS)/YY
    return P
end

# This is the switch to select the appropriate function to run
function _lombscargle{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                           signal::AbstractVector{R2},
                                                           floatrange::Bool,
                                                           with_errors::Bool,
                                                           w::AbstractVector{R3}=ones(signal)/length(signal);
                                                           normalization::AbstractString="standard",
                                                           noise_level::Real=1.0,
                                                           center_data::Bool=true,
                                                           fit_mean::Bool=true,
                                                           oversampling::Integer=5,
                                                           Mfft::Integer=4,
                                                           samples_per_peak::Integer=5,
                                                           nyquist_factor::Integer=5,
                                                           minimum_frequency::Real=NaN,
                                                           maximum_frequency::Real=NaN,
                                                           frequencies::AbstractVector{R4}=
                                                           autofrequency(times,
                                                                         samples_per_peak=samples_per_peak,
                                                                         nyquist_factor=nyquist_factor,
                                                                         minimum_frequency=minimum_frequency,
                                                                         maximum_frequency=maximum_frequency),
                                                           fast::Bool=(length(frequencies) > 200))
    @assert length(times) == length(signal) == length(w)
    # By default, we will use the fast algorithm only if there are more than 200
    # frequencies.  In any case, times must have been passed as a Range.
    if floatrange && fast
        P = _press_rybicki(times, signal, w, frequencies, center_data,
                              fit_mean, oversampling, Mfft)
    else
        if fit_mean || with_errors
            P = _generalised_lombscargle(times, signal, w, frequencies,
                                            center_data, fit_mean)
        else
            P = _lombscargle_orig(times, signal, frequencies, center_data)
        end
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

    return Periodogram(P, frequencies, times, normalization)
end

# No uncertainties
lombscargle{R1<:Real,R2<:Real}(times::Range{R1},
                               signal::AbstractVector{R2};
                               kwargs...) =
                                   _lombscargle(times,
                                                signal,
                                                true,
                                                false;
                                                kwargs...)

lombscargle{R1<:Real,R2<:Real}(times::AbstractVector{R1},
                               signal::AbstractVector{R2};
                               kwargs...) =
                                   _lombscargle(times,
                                                signal,
                                                false,
                                                false;
                                                kwargs...)

# Uncertainties provided
function _lombscargle_with_errors{R1<:Real,R2<:Real,R3<:Real}(times::AbstractVector{R1},
                                                              signal::AbstractVector{R2},
                                                              errors::AbstractVector{R3},
                                                              floatrange::Bool;
                                                              kwargs...)
    # Compute weights vector
    w = 1.0./errors.^2
    w /= sum(w)
    return _lombscargle(times, signal, floatrange, true, w; kwargs...)
end

function lombscargle{R1<:Real,R2<:Real,R3<:Real}(times::Range{R1},
                                                 signal::AbstractVector{R2},
                                                 errors::AbstractVector{R3};
                                                 kwargs...)
    return _lombscargle_with_errors(times, signal, errors, true; kwargs...)
end

function lombscargle{R1<:Real,R2<:Real,R3<:Real}(times::AbstractVector{R1},
                                                 signal::AbstractVector{R2},
                                                 errors::AbstractVector{R3};
                                                 kwargs...)
    return _lombscargle_with_errors(times, signal, errors, false; kwargs...)
end

# Uncertainties provided via Measurement type
function lombscargle{T<:Real,F<:AbstractFloat}(times::Range{T},
                                               signal::AbstractVector{Measurement{F}};
                                               kwargs...)
    return _lombscargle_with_errors(times,
                                    value(signal),
                                    uncertainty(signal),
                                    true;
                                    kwargs...)
end

function lombscargle{T<:Real,F<:AbstractFloat}(times::AbstractVector{T},
                                               signal::AbstractVector{Measurement{F}};
                                               kwargs...)
    return _lombscargle_with_errors(times,
                                    value(signal),
                                    uncertainty(signal),
                                    false;
                                    kwargs...)
end

"""
    lombscargle(times::AbstractVector{Real}, signal::AbstractVector{Real},
                errors::AbstractVector{Real}=ones(signal);
                normalization::AbstractString="standard",
                noise_level::Real=1.0,
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
* `fast`: whether to use the fast method by Press & Rybicki, overriding the
  default choice.  In any case, `times` must be a `Range` object in order to use
  this method
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
lombscargle

end # module
