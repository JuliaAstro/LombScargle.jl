### extirpolation.jl
#
# Copyright (C) 2016, 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle, extirpolation
#
# This file is a part of LombScargle.jl.
#
# License is MIT "Expat".
#
### Commentary:
#
# Utility functions used for implementing the fast algorithm proposed here:
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# The actual implementation is adapted from Astropy project.  See
# stats/lombscargle/implementations/utils.py in Astropy source code.
#
### Code:

function add_at!(arr::AbstractVector{T1}, ind::AbstractVector{<:Integer},
                 vals::AbstractVector{T2}) where {T1,T2}
    @inbounds for i in eachindex(ind)
        arr[ind[i]] += vals[i]
    end
end

function extirpolate!(result, x::AbstractVector{<:Real}, y::AbstractVector{NU},
                      N::Integer, M::Integer=4) where {NU<:Number}
    # Now use legendre polynomial weights to populate the results array; This is
    # an efficient recursive implementation (See Press et al. 1989)
    fill!(result, zero(NU))
    # first take care of the easy cases where x is an integer
    integers = find(isinteger, x)
    add_at!(result, mod.(trunc.(Int, x[integers]), N) .+ 1, y[integers])
    deleteat!(x, integers)
    deleteat!(y, integers)
    # For each remaining x, find the index describing the extirpolation range.
    # i.e. ilo[i] < x[i] < ilo[i] + M with x[i] in the center, adjusted so that
    # the limits are within the range 0...N
    ilo = clamp!(trunc.(Int, x .- div(M, 2)), 0, N - M)
    v = collect(0:(M - 1))
    numerator = [y[j] * prod(x[j] - ilo[j] - v) for j in eachindex(x)]
    denominator = float(factorial(M - 1))
    ilo .+= M .+ 1
    @inbounds for j in v
        if j > 0
            denominator *= j/(j - M)
        end
        ilo .-= 1
        add_at!(result, ilo, numerator ./ (denominator .* (x .- ilo .+ 1)))
    end
    return result
end

function trig_sum!(grid, bfft_vec, bfft_plan, t::AbstractVector{<:Real},
                   h::AbstractVector{<:Real}, df::Real, N::Integer,
                   Nfft::Integer, t0::Real, f0::Real=0.0,
                   freq_factor::Integer=1, Mfft::Integer=4)
    df *= freq_factor
    f0 *= freq_factor
    if f0 > 0
        H = h .* cis.(f0 .* (t .- t0) .* 2 .* pi)
    else
        H = complex(h)
    end
    tnorm = mod.(((t .- t0) .* Nfft .* df), Nfft)
    extirpolate!(grid, tnorm, H, Nfft, Mfft)
    A_mul_B!(bfft_vec, bfft_plan, grid)
    fftgrid = @view(bfft_vec[1:N])
    if t0 != 0
        fftgrid .*= cis.(t0 .* (f0 .+ df .* (0:(N - 1))) .* 2 .* pi)
    end
    return real(fftgrid), imag(fftgrid)
end
