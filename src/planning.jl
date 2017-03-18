### press-rybicki.jl ---  Planthe Lomb-Scargle periodogram
#
# Copyright (C) 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle
#
# This file is a part of LombScargle.jl.
#
# License is MIT "Expat".
#
### Code:

# Switches to select the appropriate algorithm to compute the periodogram.
function _plan_no_fast{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(times::AbstractVector{R1},
                                                            signal::AbstractVector{R2},
                                                            w::AbstractVector{R3},
                                                            frequencies::AbstractVector{R4},
                                                            with_errors::Bool,
                                                            center_data::Bool,
                                                            fit_mean::Bool,
                                                            noise::Real,
                                                            normalization::Symbol)
    P_type = promote_type(float(R1), float(R2), float(R3), float(R4))
    P = Vector{P_type}(length(frequencies))
    if fit_mean || with_errors
        # If "center_data" or "fit_mean" keywords are true,
        # subtract the weighted mean from each point.
        if center_data || fit_mean
            y = signal .- (w ⋅ signal) ./ sum(w)
        else
            y = signal
        end
        YY = w ⋅ (y .^ 2)
        if fit_mean
            Y   = w ⋅ y
            YY -= Y * Y
            return GLSPlan_fit_mean(times, signal, frequencies, w, y, Y,
                                    YY, noise, normalization, P)
        else
            return GLSPlan(times, signal, frequencies, w, y,
                           YY, noise, normalization, P)
        end
    else
        if center_data
            X = signal .- mean(signal)
        else
            X = signal
        end
        return LSPlan(times, signal, frequencies, w, X, X ⋅ X, noise, normalization, P)
    end
end

# There are two "_plan" methods, the only different argument being "frequencies".  When it
# is a `Range' object (first method) it could be possible to use the fast method, provided
# that fast == true, otherwise we can only use the non fast methods.
function _plan{R1<:Real,R2<:Real}(times::AbstractVector{<:Real}, signal::AbstractVector{R1},
                                  w::AbstractVector{R2}, frequencies::Range{<:Real},
                                  fast::Bool, with_errors::Bool, center_data::Bool,
                                  fit_mean::Bool, oversampling::Integer, Mfft::Integer,
                                  noise::Real, normalization::Symbol)
    if fast
        @assert Mfft > 0
        @assert step(frequencies) > 0
        if center_data || fit_mean
            y = signal .- w⋅signal
        else
            y = signal
        end
        YY = w ⋅ (y .^ 2)
        Nfft = nextpow2(length(frequencies) * oversampling)
        T = promote_type(float(R1), float(R2))
        bfft_vect = Vector{Complex{T}}(Nfft)
        bfft_grid = Vector{Complex{T}}(Nfft)
        bfft_plan = plan_bfft(bfft_vect, flags = FFTW.MEASURE)
        if fit_mean
            return FastGLSPlan_fit_mean(times, signal, frequencies, w, y, YY,
                                        bfft_vect, bfft_grid, bfft_plan, Mfft,
                                        noise, normalization, Vector{T}(length(frequencies)))
        else
            return FastGLSPlan(times, signal, frequencies, w, y, YY,
                               bfft_vect, bfft_grid, bfft_plan, Mfft,
                               noise, normalization, Vector{T}(length(frequencies)))
        end
    else
        return _plan_no_fast(times, signal, w, frequencies, with_errors,
                             center_data, fit_mean, noise, normalization)
    end
end

_plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},
      w::AbstractVector{<:Real}, frequencies::AbstractVector{<:Real},
      fast::Bool, with_errors::Bool, center_data::Bool, fit_mean::Bool,
      oversampling::Integer, Mfft::Integer, noise::Real, normalization::Symbol) =
          _plan_no_fast(times, signal, w, frequencies,
                        with_errors, center_data, fit_mean, noise, normalization)

# The main purpose of this function is to compute normalization of the
# periodogram computed with one of the functions above.
function _plan(times::AbstractVector{<:Real},
               signal::AbstractVector{<:Real},
               with_errors::Bool,
               w::AbstractVector{<:Real}=ones(signal)/length(signal);
               normalization::Symbol=:standard,
               noise_level::Real=1,
               center_data::Bool=true,
               fit_mean::Bool=true,
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
    return _plan(times, signal, w, frequencies, fast, with_errors, center_data,
                 fit_mean, oversampling, Mfft, noise_level, normalization)
end

### Main interface functions

# No uncertainties
plan{R1<:Real,R2<:Real}(times::AbstractVector{R1}, signal::AbstractVector{R2}; kwargs...) =
    _plan(times, signal, false; kwargs...)

# Uncertainties provided
function plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},
              errors::AbstractVector{<:Real}; kwargs...)
    # Compute weights vector
    w = 1 ./ errors .^ 2
    w ./= sum(w)
    return _plan(times, signal, true, w; kwargs...)
end

# Uncertainties provided via Measurement type
plan{T<:Real,F<:AbstractFloat}(times::AbstractVector{T},
                               signal::AbstractVector{Measurements.Measurement{F}};
                               kwargs...) =
                                   plan(times, Measurements.value.(signal),
                                        Measurements.uncertainty.(signal); kwargs...)
