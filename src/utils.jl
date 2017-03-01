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

export power, freq, freqpower, findmaxpower, findmaxfreq, period, periodpower,
    findmaxperiod, prob, probinv, fap, fapinv

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


_findmaxfreq{R1<:Real,R2<:Real}(freq::AbstractVector{R1},
                                power::AbstractVector{R2},
                                threshold::Real) =
                                    freq[find(x -> x >= threshold, power)]

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

function findmaxfreq{R<:Real}(p::Periodogram, interval::AbstractVector{R},
                              threshold::Real=NaN)
    f = freq(p)
    lo, hi = extrema(interval)
    indices = find(x -> lo <= x <= hi, f)
    pow = power(p)[indices]
    if isnan(threshold)
        threshold = maximum(pow)
    end
    return _findmaxfreq(f[indices], pow, threshold)
end

"""
    power(p::Periodogram)

Return the period vector of Lomb-Scargle periodogram `p`.  It is equal to `1./freq(p)`.
"""
period(p::Periodogram) = 1./freq(p)

"""
    periodpower(p::Periodogram)

Return the 2-tuple `(period(p), power(p))`, where `period(p)` and `power(p)` are
the period vector and the power vector of Lomb-Scargle periodogram `p`
respectively.
"""
periodpower(p::Periodogram) = (period(p), power(p))


"""
    findmaxperiod(p::Periodogram, [interval::AbstractVector{Real}], threshold::Real=findmaxpower(p))

Return the array of period with the highest power in the periodogram `p`.  If a
scalar real argument `threshold` is provided, return the period with power
larger than or equal to `threshold`.  If you want to limit the search to a
narrower period range, pass as second argument a vector with the extrema of the
interval.
"""
findmaxperiod(p::Periodogram, threshold::Real=findmaxpower(p)) =
    1./findmaxfreq(p, threshold)
findmaxperiod{R<:Real}(p::Periodogram, interval::AbstractVector{R},
                       threshold::Real=NaN) =
                           1./findmaxfreq(p, 1./interval, threshold)

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
    prob(P::Periodogram, pow::Real)

Return the probability that the periodogram power can exceed the value `pow`.

Its inverse is the `probinv` function.
"""
function prob(P::Periodogram, pow::Real)
    N = length(P.times)
    normalization = P.norm
    if normalization == :standard
        return (1.0 - pow)^((N - 3.0)*0.5)
    elseif normalization == :Scargle
        return exp(-pow)
    elseif normalization == :HorneBaliunas
        return (1.0 - 2*pow/(N - 1))^((N - 3.0)*0.5)
    elseif normalization == :Cumming
        return (1.0 + 2*pow/(N - 3.0))^((3.0 - N)*0.5)
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
        return 1.0 - prob^(2.0/(N - 3.0))
    elseif normalization == :Scargle
        return -log(prob)
    elseif normalization == :HorneBaliunas
        return 0.5*(N - 1.0)*(1.0 - prob^(2.0/(N - 3.0)))
    elseif normalization == :Cumming
        return 0.5*(N - 3.0)*(prob^(2.0/(3.0 - N)) - 1.0)
    else
        error("normalization \"", string(normalization), "\" not supported")
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
    fap(P::Periodogram, pow::Real)

Return the false-alarm probability for periodogram `P` and power value `pow`.

Its inverse is the `fapinv` function.
"""
fap(P::Periodogram, pow::Real) = 1.0 - (1.0 - prob(P, pow))^M(P)

"""
    fapinv(P::Periodogram, prob::Real)

Return the power value of the periodogram whose false-alarm probability is
`prob`.

This is the inverse of `fap` function.
"""
fapinv(P::Periodogram, prob::Real) = probinv(P, 1.0 - (1.0 - prob)^(inv(M(P))))

"""
    LombScargle.model(times::AbstractVector{Real},
                      signal::AbstractVector{R2},
                      [errors::AbstractVector{R3},]
                      frequency::Real,
                      [times_fit::AbstractVector{R4}];
                      center_data::Bool=true,
                      fit_mean::Bool=true)

Return the best fitting Lomb-Scargle model for the given signal at the given
frequency.

Mandatory arguments are:

* `times`: the observation times
* `signal`: the signal, sampled at `times` (must have the same length as
  `times`)
* `frequency`: the frequency at which to calculate the model

Optional arguments are:

* `errors`: the vector of uncertainties of the signal.  If provided, it must
  have the same length as `signal` and `times`, and be the third argument
* `times_fit`: the vector of times at which the model will be calculated.  It
  defaults to `times`.  If provided, it must come after `frequency`

Optional keyword arguments `center_data` and `fit_mean` have the same meaning as
in `lombscargle` function, which see.
"""
function model{R1<:Real,R2<:Real,R3<:Real,R4<:Real}(t::AbstractVector{R1},
                                                    s::AbstractVector{R2},
                                                    errors::AbstractVector{R3},
                                                    f::Real,
                                                    t_fit::AbstractVector{R4}=t;
                                                    center_data::Bool=true,
                                                    fit_mean::Bool=true)
    @assert length(t) == length(s) == length(errors)
    if center_data
        # Compute weights vector
        w = 1 ./ (errors .^ 2)
        m = (w⋅s)/sum(w)
        y = (s .- m) ./ errors
    else
        # We don't want to center the data: the mean is 0 and the signal is left
        # unchanged
        m = zero(R2)
        y = s./errors
    end
    # The Lomb-Scargle periodogram looks for the best fitting sinusoidal
    # function
    #     a·cos(ωt) + b·sin(ωt) + c
    # In order to find the coefficients a, b, c for the given frequency we can
    # solve the linear system A·x = y, where A is the matrix with rows:
    # [cos(ωt) sin(ωt) 1]; x is the column vector [a, b, c], and y is the
    # column vector of the signal
    ω = 2pi*f
    if fit_mean
        a, b, c = [cos.(ω .* t) sin.(ω .* t)  ones(t)] ./ errors \ y
        return a .* cos.(ω .* t_fit) .+ b .* sin.(ω .* t_fit) .+ (c + m)
    else
        # If fit_mean==false, the model to be fitted is a·cos(ωt) + b·sin(ωt)
        a, b = [cos.(ω .* t) sin.(ω .* t)] ./ errors \ y
        return a .* cos.(ω .* t_fit) .+ b .* sin.(ω .* t_fit) .+ m
    end
end

# No uncertainties: errors=ones(s)
model{R1<:Real,R2<:Real,R3<:Real}(t::AbstractVector{R1},
                                  s::AbstractVector{R2},
                                  f::Real,
                                  t_fit::AbstractVector{R3}=t;
                                  kwargs...) =
                                      model(t, s, ones(s), f, t_fit; kwargs...)

# Uncertainties provided via Measurement type
model{R1<:Real,R2<:AbstractFloat,R3<:Real}(t::AbstractVector{R1},
                                           s::AbstractVector{Measurement{R2}},
                                           f::Real,
                                           t_fit::AbstractVector{R3}=t;
                                           kwargs...) =
                                               model(t,
                                                     Measurements.value.(s),
                                                     Measurements.uncertainty.(s),
                                                     f,
                                                     t_fit;
                                                     kwargs...)
