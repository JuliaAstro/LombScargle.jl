### press-rybicki.jl ---  Perform fast but approximate Lomb–Scargle periodogram
#
# Copyright (C) 2016 Astropy Developers
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

struct FastGLSPlan{T,A,B<:AbstractVector{T},C,D,E,F,G,H,I,J,K,L} <: PeriodogramPlan
    times::A
    signal::B
    freq::C
    sumw::D
    w::E
    y::B
    YY::T
    fftgrid::F
    bfft_vect::G
    bfft_grid::H
    bfft_plan::I
    Mfft::J
    noise::K
    norm::Symbol
    P::L
end

struct FastGLSPlan_fit_mean{T,A,B<:AbstractVector{T},C,D,E,F,G,H,I,J,K,L} <: PeriodogramPlan
    times::A
    signal::B
    freq::C
    sumw::D
    w::E
    y::B
    YY::T
    fftgrid::F
    bfft_vect::G
    bfft_grid::H
    bfft_plan::I
    Mfft::J
    noise::K
    norm::Symbol
    P::L
end

include("extirpolation.jl")

# Fast, but approximate, method to compute the Lomb–Scargle periodogram for
# evenly spaced frequency grid.  See
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# This is adapted from Astropy implementation of the method.  See
# `lombscargle_fast' function.

function _press_rybicki!(P, times::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                         w::AbstractVector{<:Real}, t0, df, N , f0, YY::Real,
                         fftgrid, bfft_vec::AbstractVector{Complex{T}},
                         grid, plan, Nfft::Integer, Mfft::Integer) where {T<:AbstractFloat}
    #---------------------------------------------------------------------------
    # 1. compute functions of the time-shift tau at each frequency
    Ch, Sh = trig_sum!(grid, fftgrid, bfft_vec, plan, times, w .* y, df, N, Nfft, t0, f0, 1, Mfft)
    C2, S2 = trig_sum!(grid, fftgrid, bfft_vec, plan, times, w,      df, N, Nfft, t0, f0, 2, Mfft)
    #---------------------------------------------------------------------------
    # 2. Compute the periodogram, following Zechmeister & Kurster
    #    and using tricks from Press & Rybicki.
    tan_2ωτ = S2 ./ C2
    C2w = 1 ./ (sqrt.(1 .+ tan_2ωτ .* tan_2ωτ)) # = cos(2 * ωτ)
    S2w = tan_2ωτ .* C2w # = sin(2 * ωτ)
    Cw  = sqrt.((1 .+ C2w) ./ 2) # = cos(ωτ)
    Sw  = sign.(S2w) .* sqrt.((1 .- C2w) ./ 2) # = sin(ωτ)
    return P .= ((Ch .* Cw .+ Sh .* Sw) .^ 2 ./ ((1 .+ C2 .* C2w .+ S2 .* S2w) ./ 2) .+
                 (Sh .* Cw .- Ch .* Sw) .^ 2 ./ ((1 .- C2 .* C2w .- S2 .* S2w) ./ 2)) ./ YY
end

function _press_rybicki_fit_mean!(P, times::AbstractVector{<:Real},
                                  y::AbstractVector{<:Real}, w::AbstractVector{<:Real},
                                  t0, df, N , f0, YY::Real, fftgrid,
                                  bfft_vec::AbstractVector{Complex{T}}, grid, plan,
                                  Nfft::Integer, Mfft::Integer) where {T<:AbstractFloat}
    #---------------------------------------------------------------------------
    # 1. compute functions of the time-shift tau at each frequency
    Ch, Sh = trig_sum!(grid, fftgrid, bfft_vec, plan, times, w .* y, df, N, Nfft, t0, f0, 1, Mfft)
    C2, S2 = trig_sum!(grid, fftgrid, bfft_vec, plan, times, w,      df, N, Nfft, t0, f0, 2, Mfft)
    C, S   = trig_sum!(grid, fftgrid, bfft_vec, plan, times, w,      df, N, Nfft, t0, f0, 1, Mfft)
    #---------------------------------------------------------------------------
    # 2. Compute the periodogram, following Zechmeister & Kurster
    #    and using tricks from Press & Rybicki.
    tan_2ωτ = (S2 .- 2 .* S .* C) ./ (C2 .- (C .* C .- S .* S))
    C2w = 1 ./ (sqrt.(1 .+ tan_2ωτ .* tan_2ωτ)) # = cos(2 * ωτ)
    S2w = tan_2ωτ .* C2w # = sin(2 * ωτ)
    Cw  = sqrt.((1 .+ C2w) ./ 2) # = cos(ωτ)
    Sw  = sign.(S2w) .* sqrt.((1 .- C2w) ./ 2) # = sin(ωτ)
    return P .= ((Ch .* Cw .+ Sh .* Sw) .^ 2 ./
                 ((1 .+ C2 .* C2w .+ S2 .* S2w) ./ 2 .- (C .* Cw .+ S .* Sw) .^ 2) .+
                 (Sh .* Cw .- Ch .* Sw) .^ 2 ./
                 ((1 .- C2 .* C2w .- S2 .* S2w) ./ 2 .- (S .* Cw .- C .* Sw) .^ 2)) ./ YY
end

_periodogram!(P::AbstractVector, times, p::FastGLSPlan) =
    _press_rybicki!(P, times, p.y, p.w, minimum(p.times), step(p.freq),
                    length(p.freq), minimum(p.freq), p.YY, p.fftgrid, p.bfft_vect,
                    p.bfft_grid, p.bfft_plan, length(p.bfft_vect), p.Mfft)
_periodogram!(p::FastGLSPlan) = _periodogram!(p.P, p.times, p)
_periodogram!(P::AbstractVector, times, p::FastGLSPlan_fit_mean) =
    _press_rybicki_fit_mean!(P, times, p.y, p.w, minimum(p.times), step(p.freq),
                             length(p.freq), minimum(p.freq), p.YY, p.fftgrid, p.bfft_vect,
                             p.bfft_grid, p.bfft_plan, length(p.bfft_vect), p.Mfft)
_periodogram!(p::FastGLSPlan_fit_mean) = _periodogram!(p.P, p.times, p)
