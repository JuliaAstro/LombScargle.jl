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
immutable Periodogram{P<:AbstractFloat, F<:AbstractVector,
                      T<:AbstractVector, S<:AbstractString}
    power::Vector{P}
    freq::F
    # XXX: the `times' vector is only in the `M' function (see utils.jl), but
    # only maximum(times) and minimum(times) are used.  We could consider the
    # possibility to keep in this type only the extrema of t, instead of the
    # whole array.
    times::T
    norm::S
end

include("extirpolation.jl")
include("utils.jl")
include("bootstrap.jl")

# XXX: Until Julia's bug #15276
# (https://github.com/JuliaLang/julia/issues/15276) will be fixed, we have to
# move the body loop to an external function in order to effectively gain
# scaling when multiple threads are used.

function _lombscargle_orig_loop!(P, freqs, times, X, XX, nil, N)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = 2pi*freqs[n]
        XC = XS = CC = CS = nil
        @inbounds for j in eachindex(times)
            ωt = ω*times[j]
            C = cos(ωt)
            S = sin(ωt)
            XC += X[j]*C
            XS += X[j]*S
            CC += C*C
            CS += C*S
        end
        SS      = N - CC
        ωτ      = atan2(CS, CC - N / 2) / 2
        c_τ     = cos(ωτ)
        s_τ     = sin(ωτ)
        c_τ2    = c_τ*c_τ
        s_τ2    = s_τ*s_τ
        cs_τ_CS = 2c_τ*s_τ*CS
        P[n] = (abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))/XX
    end
end

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
    # If "center_data" keyword is true, subtract the mean from each point.
    if center_data
        X = signal - mean(signal)
    else
        X = signal
    end
    XX = X⋅X
    _lombscargle_orig_loop!(P, freqs, times, X, XX, zero(P_type), length(signal))
    return P
end

function _generalised_lombscargle_loop!(P, freqs, times, y, w, Y, YY, nil)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = 2pi*freqs[n]
        # Find τ for current angular frequency
        C = S = CS = CC = nil
        @inbounds for i in eachindex(times)
            ωt  = ω*times[i]
            W   = w[i]
            c   = cos(ωt)
            s   = sin(ωt)
            CS += W*c*s
            CC += W*c*c
            C  += W*c
            S  += W*s
        end
        CS -= C*S
        SS  = 1 - CC - S*S
        CC -= C*C
        ωτ   = atan2(2CS, CC - SS) / 2
        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = nil
        @inbounds for i in eachindex(times)
            ωt    = ω*times[i] - ωτ
            W     = w[i]
            c     = cos(ωt)
            s     = sin(ωt)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
        end
        # "C_τ" and "S_τ" are computed following equation (7) of Press &
        # Rybicki, this formula simply comes from angle difference trigonometric
        # identities.
        cos_ωτ = cos(ωτ)
        sin_ωτ = sin(ωτ)
        C_τ    = C*cos_ωτ + S*sin_ωτ
        S_τ    = S*cos_ωτ - C*sin_ωτ
        YC_τ  -= Y*C_τ
        YS_τ  -= Y*S_τ
        SS_τ   = 1 - CC_τ - S_τ*S_τ
        CC_τ  -= C_τ*C_τ
        P[n] = (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY
    end
end

function _generalised_lombscargle_loop!(P, freqs, times, y, w, YY, nil)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = 2pi*freqs[n]
        # Find τ for current angular frequency
        C = S = CS = CC = nil
        @inbounds for i in eachindex(times)
            ωt  = ω*times[i]
            W   = w[i]
            c   = cos(ωt)
            s   = sin(ωt)
            CS += W*c*s
            CC += W*c*c
        end
        SS  = 1 - CC
        ωτ   = atan2(2CS, CC - SS) / 2
        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = nil
        @inbounds for i in eachindex(times)
            ωt    = ω*times[i] - ωτ
            W     = w[i]
            c     = cos(ωt)
            s     = sin(ωt)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
        end
        SS_τ  = 1 - CC_τ
        P[n] = (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY
    end
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
    # If "center_data" or "fit_mean" keywords are true,
    # subtract the weighted mean from each point.
    if center_data || fit_mean
        y = signal - (w ⋅ signal)/sum(w)
    else
        y = signal
    end
    YY = w⋅(y.^2)
    if fit_mean
        Y   = w⋅y
        YY -= Y*Y
        _generalised_lombscargle_loop!(P, freqs, times, y, w, Y, YY, nil)
    else
        _generalised_lombscargle_loop!(P, freqs, times, y, w, YY, nil)
    end
    return P
end

# Fast, but approximate, method to compute the Lomb-Scargle periodogram for
# evenly spaced frequency grid.  See
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# This is adapted from Astropy implementation of the method.  See
# `lombscargle_fast' function.
function _press_rybicki{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                             signal::AbstractVector{R2},
                                                             w::AbstractVector{R3},
                                                             freqs::Range{R4},
                                                             center_data::Bool,
                                                             fit_mean::Bool,
                                                             oversampling::Int,
                                                             Mfft::Int)
    @assert Mfft > 0
    nil = zero(promote_type(float(R1), float(R2)))
    # If "center_data" keyword is true, subtract the weighted mean from each
    # point.
    if center_data || fit_mean
        y = signal - w⋅signal
    else
        y = signal
    end
    df = step(freqs)
    N  = length(freqs)
    f0 = minimum(freqs)
    #---------------------------------------------------------------------------
    # 0. prepare for the computation of IFFT.
    # `ifft' can take arrays of any length in input, but
    # it's faster when the length is exactly a power of 2.
    Nfft = nextpow2(N * oversampling)
    ifft_vec = Vector{Complex{promote_type(float(R2), float(R3))}}(Nfft)
    grid = similar(ifft_vec)
    plan = plan_ifft(ifft_vec, flags = FFTW.MEASURE)
    #---------------------------------------------------------------------------
    # 1. compute functions of the time-shift tau at each frequency
    Ch, Sh = trig_sum!(grid, ifft_vec, plan, times, w .* y, df, N, Nfft, f0,
                       1, Mfft)
    C2, S2 = trig_sum!(grid, ifft_vec, plan, times, w,      df, N, Nfft, f0,
                       2, Mfft)
    if fit_mean
        C, S = trig_sum!(grid, ifft_vec, plan, times, w, df,    N, Nfft, f0,
                         1, Mfft)
        tan_2ωτ = (S2 .- 2 .* S .* C) ./ (C2 .- (C .* C .- S .* S))
    else
        tan_2ωτ = S2 ./ C2
    end
    # This is what we're computing below; the straightforward way is slower and
    # less stable, so we use trig identities instead
    #   ωτ = 0.5 * atan(tan_2ωτ)
    #   S2w, C2w = sin(2 * ωτ), cos(2 * ωτ)
    #   Sw, Cw = sin(ωτ), cos(ωτ)
    C2w = 1 ./ (sqrt.(1 .+ tan_2ωτ .* tan_2ωτ))
    S2w = tan_2ωτ .* C2w
    Cw = sqrt.((1 .+ C2w) ./ 2)
    Sw = sign(S2w) .* sqrt.((1 .- C2w) ./ 2)
    #---------------------------------------------------------------------------
    # 2. Compute the periodogram, following Zechmeister & Kurster
    #    and using tricks from Press & Rybicki.
    YY = w⋅(y.^2)
    YC = Ch .* Cw .+ Sh .* Sw
    YS = Sh .* Cw .- Ch .* Sw
    CC = (1 .+ C2 .* C2w .+ S2 .* S2w) ./ 2
    SS = (1 .- C2 .* C2w .- S2 .* S2w) ./ 2 # = 1 .- CC
    if fit_mean
        CC .-= (C .* Cw .+ S .* Sw) .^ 2
        SS .-= (S .* Cw .- C .* Sw) .^ 2
    end
    return (YC .* YC ./ CC .+ YS .* YS ./ SS) ./ YY
end

# Switches to select the appropriate algorithm to compute the periodogram.
function periodogram_no_fast{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                                  signal::AbstractVector{R2},
                                                                  w::AbstractVector{R3},
                                                                  frequencies::AbstractVector{R4},
                                                                  with_errors::Bool,
                                                                  center_data::Bool,
                                                                  fit_mean::Bool)
    if fit_mean || with_errors
        return _generalised_lombscargle(times, signal, w, frequencies,
                                        center_data, fit_mean)
    else
        return _lombscargle_orig(times, signal, frequencies, center_data)
    end
end

# There are two "periodogram" methods, the only different argument being
# "frequencies".  When it is a `Range' object (first method) it could be
# possible to use the fast method, provided that fast == true, otherwise we can
# only use the non fast methods.
function periodogram{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                          signal::AbstractVector{R2},
                                                          w::AbstractVector{R3},
                                                          frequencies::Range{R4},
                                                          fast::Bool,
                                                          with_errors::Bool,
                                                          center_data::Bool,
                                                          fit_mean::Bool,
                                                          oversampling::Integer,
                                                          Mfft::Integer)
    if fast
        return _press_rybicki(times, signal, w, frequencies, center_data,
                              fit_mean, oversampling, Mfft)
    else
        return periodogram_no_fast(times, signal, w, frequencies,
                                   with_errors, center_data, fit_mean)
    end
end

function periodogram{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                          signal::AbstractVector{R2},
                                                          w::AbstractVector{R3},
                                                          frequencies::AbstractVector{R4},
                                                          fast::Bool,
                                                          with_errors::Bool,
                                                          center_data::Bool,
                                                          fit_mean::Bool,
                                                          oversampling::Integer,
                                                          Mfft::Integer)
    return periodogram_no_fast(times, signal, w, frequencies,
                               with_errors, center_data, fit_mean)
end

# The main purpose of this function is to compute normalization of the
# periodogram computed with one of the functions above.
function _lombscargle{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                           signal::AbstractVector{R2},
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
    P = periodogram(times, signal, w, frequencies, fast, with_errors,
                    center_data, fit_mean, oversampling, Mfft)
    N = length(signal)
    # Normalize periodogram
    if normalization == "standard"
    elseif normalization == "model"
        P = P./(1.0 - P)
    elseif normalization == "log"
        P = -log.(1.0 - P)
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

### Main interface functions

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
function lombscargle{T<:Real,F<:AbstractFloat}(times::AbstractVector{T},
                                               signal::AbstractVector{Measurement{F}};
                                               kwargs...)
    return lombscargle(times, value.(signal), uncertainty.(signal); kwargs...)
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
lombscargle

end # module
