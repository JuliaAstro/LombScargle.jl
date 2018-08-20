### press-rybicki.jl ---  Plan the Lomb–Scargle periodogram
#
# Copyright (C) 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle
#
# This file is a part of LombScargle.jl.
#
# License is BSD 3-clause "New" or "Revised".
#
### Code:

# Switches to select the appropriate algorithm to compute the periodogram.
function _plan_no_fast(times::AbstractVector{R1}, signal::AbstractVector{R2}, sumw::Real,
                       w::AbstractVector{R3}, frequencies::AbstractVector{R4},
                       with_errors::Bool, center_data::Bool,
                       fit_mean::Bool, noise::Real,
                       normalization::Symbol) where {R1<:Real,R2<:Real,R3<:Real,R4<:Real}
    P_type = promote_type(float(R1), float(R2), float(R3), float(R4))
    P = Vector{P_type}(undef, length(frequencies))
    if fit_mean || with_errors
        # If "center_data" or "fit_mean" keywords are true,
        # subtract the weighted mean from each point.
        if center_data || fit_mean
            y = signal .- (w ⋅ signal)
        else
            y = signal
        end
        YY = w ⋅ (y .^ 2)
        if fit_mean
            Y   = w ⋅ y
            YY -= Y * Y
            return GLSPlan_fit_mean(times, signal, frequencies, sumw, w, y, Y,
                                    YY, noise, normalization, P)
        else
            return GLSPlan(times, signal, frequencies, sumw, w, y,
                           YY, noise, normalization, P)
        end
    else
        if center_data
            X = signal .- mean(signal)
        else
            X = signal
        end
        return LSPlan(times, signal, frequencies, sumw, w, X, X ⋅ X, noise, normalization, P)
    end
end

# There are two "_plan" methods, the only different argument being "frequencies".  When it
# is a `Range' object (first method) it could be possible to use the fast method, provided
# that fast == true, otherwise we can only use the non fast methods.
function _plan(times::AbstractVector{<:Real}, signal::AbstractVector{R1}, sumw::Real,
               w::AbstractVector{R2}, frequencies::AbstractRange{<:Real}, fast::Bool,
               with_errors::Bool, center_data::Bool, fit_mean::Bool, flags::Integer,
               timelimit::Real, oversampling::Integer, Mfft::Integer, noise::Real,
               normalization::Symbol) where {R1<:Real,R2<:Real}
    if fast
        @assert Mfft > 0
        @assert step(frequencies) > 0
        if center_data || fit_mean
            y = signal .- w⋅signal
        else
            y = signal
        end
        YY = w ⋅ (y .^ 2)
        N = length(frequencies)
        Nfft = nextpow(2, N * oversampling)
        T = promote_type(float(R1), float(R2))
        bfft_vect = Vector{Complex{T}}(undef, Nfft)
        bfft_grid = Vector{Complex{T}}(undef, Nfft)
        fftgrid   = Vector{Complex{T}}(undef, N)
        bfft_plan = FFTW.plan_bfft(bfft_vect, flags = flags, timelimit = timelimit)
        if fit_mean
            return FastGLSPlan_fit_mean(times, signal, frequencies, sumw, w, y, YY,
                                        fftgrid, bfft_vect, bfft_grid, bfft_plan, Mfft,
                                        noise, normalization, Vector{T}(undef, N))
        else
            return FastGLSPlan(times, signal, frequencies, sumw, w, y, YY,
                               fftgrid, bfft_vect, bfft_grid, bfft_plan, Mfft,
                               noise, normalization, Vector{T}(undef, N))
        end
    else
        return _plan_no_fast(times, signal, sumw, w, frequencies, with_errors,
                             center_data, fit_mean, noise, normalization)
    end
end

_plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real}, sumw::Real,
      w::AbstractVector{<:Real}, frequencies::AbstractVector{<:Real},
      fast::Bool, with_errors::Bool, center_data::Bool, fit_mean::Bool, flags::Integer,
      timelimit::Real, oversampling::Integer, Mfft::Integer,
      noise::Real, normalization::Symbol) =
          _plan_no_fast(times, signal, sumw, w, frequencies,
                        with_errors, center_data, fit_mean, noise, normalization)

function _plan(times::AbstractVector{<:Real},
               signal::AbstractVector{<:Real},
               with_errors::Bool,
               sumw::Real=length(signal),
               w::AbstractVector{<:Real}=fill!(similar(signal), one(eltype(signal)))/sumw;
               normalization::Symbol=:standard,
               noise_level::Real=1,
               center_data::Bool=true,
               fit_mean::Bool=true,
               flags::Integer=FFTW.ESTIMATE,
               timelimit::Real=Inf,
               oversampling::Integer=5,
               Mfft::Integer=4,
               samples_per_peak::Integer=5,
               nyquist_factor::Integer=5,
               minimum_frequency::Real=NaN,
               maximum_frequency::Real=NaN,
               frequencies::AbstractVector{<:Real}=
               autofrequency(times,
                             samples_per_peak=samples_per_peak,
                             nyquist_factor=nyquist_factor,
                             minimum_frequency=minimum_frequency,
                             maximum_frequency=maximum_frequency),
               fast::Bool=(length(frequencies) > 200))
    @assert length(times) == length(signal) == length(w)
    return _plan(times, signal, sumw, w, frequencies, fast, with_errors, center_data,
                 fit_mean, flags, timelimit, oversampling, Mfft, noise_level, normalization)
end

### Main interface functions

# No uncertainties
plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real}; kwargs...) =
    _plan(times, signal, false; kwargs...)

# Uncertainties provided
function plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},
              errors::AbstractVector{<:Real}; kwargs...)
    # Compute weights vector
    w = 1 ./ errors .^ 2
    sumw = sum(w)
    w ./= sumw
    return _plan(times, signal, true, sumw, w; kwargs...)
end

# Uncertainties provided via Measurement type
plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Measurement}; kwargs...) =
         plan(times, Measurements.value.(signal), Measurements.uncertainty.(signal); kwargs...)

"""
    LombScargle.plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},
                     [errors::AbstractVector{<:Real}];
                     normalization::Symbol=:standard,
                     noise_level::Real=1,
                     center_data::Bool=true,
                     fit_mean::Bool=true,
                     fast::Bool=true,
                     flags::Integer=FFTW.ESTIMATE,
                     timelimit::Real=Inf,
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

Pre-plan the Lomb–Scargle periodogram of the `signal` vector, observed at
`times`.  The periodogram can then be computed by passing the result of this
function to `lombscargle`.

You can also specify the uncertainties for each signal point with `errors`
argument.  All these vectors must have the same length.

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

You can explicitely require to use or not the fast method by Press & Rybicki,
overriding the default choice, by setting the `fast` keyword.  In any case,
`frequencies` must be a `Range` object (this is the default) in order to
actually use this method.  A few other keywords are available to adjust the
settings of the periodogram when the fast method is used (otherwise they are
ignored):

* `fast`: whether to use the fast method.
* `flags`: this integer keyword is a bitwise-or of FFTW planner flags,
  defaulting to `FFTW.ESTIMATE`.  Passing `FFTW.MEASURE` or `FFTW.PATIENT` will
  instead spend several seconds (or more) benchmarking different possible FFT
  algorithms and picking the fastest one; see the FFTW manual for more
  information on planner flags.
* `timelimit`: specifies a rough upper bound on the allowed planning time, in seconds.
* `oversampling`: oversampling the frequency factor for the approximation;
  roughly the number of time samples across the highest-frequency sinusoid.
  This parameter contains the tradeoff between accuracy and speed.
* `Mfft`: the number of adjacent points to use in the FFT approximation.

In addition, you can use all optional keyword arguments of
[`LombScargle.autofrequency`](@ref) function in order to tune the `frequencies`.

If the signal has uncertainties, the `signal` vector can also be a vector of
`Measurement` objects (from
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package), in
which case you don’t need to pass a separate `errors` vector for the
uncertainties of the signal.
"""
plan
