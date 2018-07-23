### utils.jl
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
### Commentary
#
# This file contains some utilities for the LombScargle.jl package.
#
### Code:

export power, freq, freqpower, findmaxpower, findmaxfreq, period, periodpower,
    findmaxperiod, prob, probinv, fap, fapinv

"""
    power(p::Periodogram)

Return the power vector of Lomb–Scargle periodogram `p`.
"""
power(p::Periodogram) = p.power

"""
    freq(p::Periodogram)

Return the frequency vector of Lomb–Scargle periodogram `p`.
"""
freq(p::Periodogram) = p.freq

"""
    freqpower(p::Periodogram)

Return the 2-tuple `(freq(p), power(p))`, where `freq(p)` and `power(p)` are the
frequency vector and the power vector of Lomb–Scargle periodogram `p`
respectively.
"""
freqpower(p::Periodogram) = (freq(p), power(p))

"""
    findmaxpower(p::Periodogram)

Return the highest power of the periodogram `p`.
"""
findmaxpower(p::Periodogram) = maximum(power(p))


_findmaxfreq(freq::AbstractVector{<:Real}, power::AbstractVector{<:Real}, threshold::Real) =
    freq[findall(x -> x >= threshold, power)]

"""
    findmaxfreq(p::Periodogram, [interval::AbstractVector{Real}], threshold::Real=findmaxpower(p))

Return the array of frequencies with the highest power in the periodogram `p`.
If a scalar real argument `threshold` is provided, return the frequencies with
power larger than or equal to `threshold`.  If you want to limit the search to a
narrower frequency range, pass as second argument a vector with the extrema of
the interval.
"""
findmaxfreq(p::Periodogram, threshold::Real=findmaxpower(p)) =
    _findmaxfreq(freqpower(p)..., threshold)

function findmaxfreq(p::Periodogram, interval::AbstractVector{<:Real}, threshold::Real=NaN)
    f = freq(p)
    lo, hi = extrema(interval)
    indices = findall(x -> lo <= x <= hi, f)
    pow = power(p)[indices]
    if isnan(threshold)
        threshold = maximum(pow)
    end
    return _findmaxfreq(f[indices], pow, threshold)
end

"""
    power(p::Periodogram)

Return the period vector of Lomb–Scargle periodogram `p`.  It is equal to `1 ./ freq(p)`.
"""
period(p::Periodogram) = 1 ./ freq(p)

"""
    periodpower(p::Periodogram)

Return the 2-tuple `(period(p), power(p))`, where `period(p)` and `power(p)` are
the period vector and the power vector of Lomb–Scargle periodogram `p`
respectively.
"""
periodpower(p::Periodogram) = (period(p), power(p))


"""
    findmaxperiod(p::Periodogram, [interval::AbstractVector{Real}], threshold::Real=findmaxpower(p))

Return the array of periods with the highest power in the periodogram `p`.  If a
scalar real argument `threshold` is provided, return the period with power
larger than or equal to `threshold`.  If you want to limit the search to a
narrower period range, pass as second argument a vector with the extrema of the
interval.
"""
findmaxperiod(p::Periodogram, threshold::Real=findmaxpower(p)) =
    1 ./ findmaxfreq(p, threshold)
findmaxperiod(p::Periodogram, interval::AbstractVector{<:Real},
              threshold::Real=NaN) =
                  1 ./ findmaxfreq(p, 1 ./ interval, threshold)

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
<https://jakevdp.github.io/blog/2015/06/13/lomb-scargle-in-python/> and uses the
same keywords names adopted in Astropy.
"""
function autofrequency(times::AbstractVector{<:Real};
                       samples_per_peak::Integer=5,
                       nyquist_factor::Integer=5,
                       minimum_frequency::Real=NaN,
                       maximum_frequency::Real=NaN)
    T = maximum(times) - minimum(times)
    δf = inv(samples_per_peak * T)
    f_min = isfinite(minimum_frequency) ? minimum_frequency : (δf / 2)
    if isfinite(maximum_frequency)
        return f_min:δf:maximum_frequency
    else
        return f_min:δf:(nyquist_factor * length(times) / 2T)
    end
end

"""
    prob(P::Periodogram, pow::Real)

Return the probability that the periodogram power can exceed the value `pow`.

Its inverse is the `probinv` function.
"""
function prob(P::Periodogram, pow::Real)
    N = length(P.times)
    normalization = P.norm
    if normalization == :standard
        return (1 - pow)^((N - 3) / 2)
    elseif normalization == :Scargle
        return exp(-pow)
    elseif normalization == :HorneBaliunas
        return (1 - 2*pow/(N - 1))^((N - 3) / 2)
    elseif normalization == :Cumming
        return (1 + 2*pow/(N - 3))^((3 - N) / 2)
    else
        error("normalization \"", string(normalization), "\" not supported")
    end
end

"""
    probinv(P::Periodogram, prob::Real)

Return the power value of the periodogram power whose probability is `prob`.

This is the inverse of `prob` function.
"""
function probinv(P::Periodogram, prob::Real)
    N = length(P.times)
    normalization = P.norm
    if normalization == :standard
        return 1 - prob^(2/(N - 3))
    elseif normalization == :Scargle
        return -log(prob)
    elseif normalization == :HorneBaliunas
        return (N - 1) * (1 - prob ^ (2 / (N - 3))) / 2
    elseif normalization == :Cumming
        return (N - 3) * (prob ^ (2 / (3 - N)) - 1) / 2
    else
        error("normalization \"", string(normalization), "\" not supported")
    end
end

"""
    LombScargle.M(P::Periodogram)

Estimates the number of independent frequencies in the periodogram `P`.
"""
function M(P::Periodogram)
    tmin, tmax = extrema(P.times)
    fmin, fmax = extrema(P.freq)
    return (tmax - tmin)*(fmax - fmin)
end

"""
    fap(P::Periodogram, pow::Real)

Return the false-alarm probability for periodogram `P` and power value `pow`.

Its inverse is the [`fapinv`](@ref) function.
"""
fap(P::Periodogram, pow::Real) = 1 - (1 - prob(P, pow))^M(P)

"""
    fapinv(P::Periodogram, prob::Real)

Return the power value of the periodogram whose false-alarm probability is
`prob`.

This is the inverse of [`fap`](@ref) function.
"""
fapinv(P::Periodogram, prob::Real) = probinv(P, 1 - (1 - prob)^(inv(M(P))))

"""
    LombScargle.model(times::AbstractVector{Real},
                      signal::AbstractVector{R2},
                      [errors::AbstractVector{R3},]
                      frequency::Real,
                      [times_fit::AbstractVector{R4}];
                      center_data::Bool=true,
                      fit_mean::Bool=true)

Return the best fitting Lomb–Scargle model for the given signal at the given
frequency.

Mandatory arguments are:

* `times`: the observation times
* `signal`: the signal, sampled at `times` (must have the same length as
  `times`)
* `frequency`: the frequency at which to calculate the model

Optional arguments are:

* `errors`: the vector of uncertainties of the signal.  If provided, it must
  have the same length as `signal` and `times`, and be the third argument.  Like
  for [`lombscargle`](@ref), if the signal has uncertainties, the `signal`
  vector can also be a vector of `Measurement` objects, and this argument should
  be omitted

* `times_fit`: the vector of times at which the model will be calculated.  It
  defaults to `times`.  If provided, it must come after `frequency`

Optional keyword arguments `center_data` and `fit_mean` have the same meaning as
in [`lombscargle`](@ref):

- `fit_mean`: whether to fit for the mean. If this is `false`, like in the
  original Lomb--Scargle periodogram, ``\\mathbf{A}`` does not have the third
  column of ones, ``c_f`` is set to ``0`` and the unknown vector to be determined
  becomes ``x = [a_f, b_f]^\\text{T}``
- `center_data`: whether the data should be pre-centered before solving the
  linear system. This is particularly important if `fit_mean=false`
"""
function model(t::AbstractVector{<:Real}, s::AbstractVector{T},
               errors::AbstractVector{<:Real}, f::Real,
               t_fit::AbstractVector{<:Real}=t;
               center_data::Bool=true, fit_mean::Bool=true) where {T<:Real}
    @assert length(t) == length(s) == length(errors)
    if center_data
        # Compute weights vector
        w = 1 ./ (errors .^ 2)
        m = (w⋅s)/sum(w)
        y = (s .- m) ./ errors
    else
        # We don't want to center the data: the mean is 0 and the signal is left
        # unchanged
        m = zero(T)
        y = s./errors
    end
    # The Lomb–Scargle periodogram looks for the best fitting sinusoidal
    # function
    #     a·cos(ωt) + b·sin(ωt) + c
    # In order to find the coefficients a, b, c for the given frequency we can
    # solve the linear system A·x = y, where A is the matrix with rows:
    # [cos(ωt) sin(ωt) 1]; x is the column vector [a, b, c], and y is the
    # column vector of the signal
    ω = 2 * f
    if fit_mean
        a, b, c = [cospi.(ω .* t) sinpi.(ω .* t)  fill(1, size(t))] ./ errors \ y
        return a .* cospi.(ω .* t_fit) .+ b .* sinpi.(ω .* t_fit) .+ (c + m)
    else
        # If fit_mean==false, the model to be fitted is a·cos(ωt) + b·sin(ωt)
        a, b = [cospi.(ω .* t) sinpi.(ω .* t)] ./ errors \ y
        return a .* cospi.(ω .* t_fit) .+ b .* sinpi.(ω .* t_fit) .+ m
    end
end

# No uncertainties: errors=ones(s)
model(t::AbstractVector{<:Real}, s::AbstractVector{<:Real},
      f::Real, t_fit::AbstractVector{<:Real}=t; kwargs...) =
          model(t, s, fill(1, size(s)), f, t_fit; kwargs...)

# Uncertainties provided via Measurement type
model(t::AbstractVector{<:Real}, s::AbstractVector{<:Measurement},
      f::Real, t_fit::AbstractVector{<:Real}=t; kwargs...) =
          model(t, Measurements.value.(s), Measurements.uncertainty.(s),
                f, t_fit; kwargs...)
