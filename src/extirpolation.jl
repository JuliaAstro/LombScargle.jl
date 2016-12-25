### extirpolation.jl
#
# Copyright (C) 2016 Mosè Giordano.
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

function add_at!{T1,I<:Integer,T2}(arr::AbstractVector{T1},
                                   ind::AbstractVector{I},
                                   vals::AbstractVector{T2})
    @inbounds for i in eachindex(ind)
        arr[ind[i]] += vals[i]
    end
end

function extirpolate!{RE<:Real,NU<:Number}(result,
                                           x::AbstractVector{RE},
                                           y::AbstractVector{NU},
                                           N::Integer, M::Integer=4)
    # Now use legendre polynomial weights to populate the results array; This is
    # an efficient recursive implementation (See Press et al. 1989)
    fill!(result, zero(NU))
    # first take care of the easy cases where x is an integer
    integers = find(isinteger, x)
    add_at!(result, mod(trunc(Int, x[integers]), N) .+ 1, y[integers])
    deleteat!(x, integers)
    deleteat!(y, integers)
    # For each remaining x, find the index describing the extirpolation range.
    # i.e. ilo[i] < x[i] < ilo[i] + M with x[i] in the center, adjusted so that
    # the limits are within the range 0...N
    ilo = clamp(trunc(Int, x .- div(M, 2)), 0, N - M)
    v = collect(0:(M - 1))
    numerator = y .* [prod(x[j] - ilo[j] - v) for j in eachindex(x)]
    denominator = float(factorial(M - 1))
    @inbounds for j in v
        if j > 0
            denominator *= j/(j - M)
        end
        ind = ilo + (M - j)
        add_at!(result, ind, numerator ./ (denominator * (x .- ind .+ 1)))
    end
    return result
end

function trig_sum!{R1<:Real,R2<:Real}(grid, ifft_vec, ifft_plan,
                                      t::AbstractVector{R1},
                                      h::AbstractVector{R2}, df::Real,
                                      N::Integer, Nfft::Integer,
                                      f0::Real=0.0,
                                      freq_factor::Integer=1,
                                      Mfft::Integer=4)
    df *= freq_factor
    f0 *= freq_factor
    @assert df > 0
    t0 = minimum(t)
    if f0 > 0
        H = h .* exp.(2im .* pi .* f0 .* (t .- t0))
    else
        H = complex(h)
    end
    tnorm = mod(((t .- t0) .* Nfft .* df), Nfft)
    extirpolate!(grid, tnorm, H, Nfft, Mfft)
    A_mul_B!(ifft_vec, ifft_plan, grid)
    fftgrid = Nfft .* ifft_vec[1:N]
    if t0 != 0
        fftgrid .*= exp.(2im .* pi .* t0 .* (f0 .+ df .* (0:N - 1)))
    end
    C = real(fftgrid)
    S = imag(fftgrid)
    return C, S
end
