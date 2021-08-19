### gls.jl ---  Perform Generalised Lomb–Scargle periodogram
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

struct GLSPlan{T,A,B<:AbstractVector{T},C,D,E,F,G} <: PeriodogramPlan
    times::A
    signal::B
    freq::C
    sumw::D
    w::E
    y::B
    YY::T
    noise::F
    norm::Symbol
    P::G
end

struct GLSPlan_fit_mean{T,A,B<:AbstractVector{T},C,D,E,F,G,H} <: PeriodogramPlan
    times::A
    signal::B
    freq::C
    sumw::D
    w::E
    y::B
    Y::F
    YY::T
    noise::G
    norm::Symbol
    P::H
end

# Generalised Lomb–Scargle algorithm: this takes into account uncertainties and
# fit the mean of the signal.  This is implemented following the recipe by
# * Zechmeister, M., Kürster, M. 2009, A&A, 496, 577  (URL:
#   http://dx.doi.org/10.1051/0004-6361:200811296,
#   Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)
# In addition, some tricks suggested by
# * Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
#   http://dx.doi.org/10.1086/167197,
#   Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
# to make computation faster are adopted.

function _generalised_lombscargle!(P, freqs, times, y, w, Y, YY, nil)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = freqs[n] * 2 * pi
        # Find τ for current angular frequency
        C = S = CS = CC = nil
        @inbounds for i in eachindex(times)
            W    = w[i]
            s, c = sincos(ω * times[i])
            CS  += W*c*s
            CC  += W*c*c
            C   += W*c
            S   += W*s
        end
        CS -= C*S
        SS  = 1 - CC - S*S
        CC -= C*C
        ωτ   = atan(2CS, CC - SS) / 2
        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = nil
        @inbounds for i in eachindex(times)
            W     = w[i]
            s, c  = sincos(ω*times[i] - ωτ)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
        end
        # "C_τ" and "S_τ" are computed following equation (7) of Press &
        # Rybicki, this formula simply comes from angle difference trigonometric
        # identities.
        sin_ωτ, cos_ωτ = sincos(ωτ)
        C_τ    = C*cos_ωτ + S*sin_ωτ
        S_τ    = S*cos_ωτ - C*sin_ωτ
        YC_τ  -= Y*C_τ
        YS_τ  -= Y*S_τ
        SS_τ   = 1 - CC_τ - S_τ*S_τ
        CC_τ  -= C_τ*C_τ
        P[n] = (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY
    end
    return P
end

function _generalised_lombscargle!(P, freqs, times, y, w, YY, nil)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = freqs[n] * 2 * pi
        # Find τ for current angular frequency
        C = S = CS = CC = nil
        @inbounds for i in eachindex(times)
            W    = w[i]
            s, c = sincos(ω * times[i])
            CS  += W*c*s
            CC  += W*c*c
        end
        ωτ   = atan(2CS, 2CC - 1) / 2
        # Now we can compute the power
        YC_τ = YS_τ = CC_τ = nil
        @inbounds for i in eachindex(times)
            W     = w[i]
            s, c  = sincos(ω*times[i] - ωτ)
            YC_τ += W*y[i]*c
            YS_τ += W*y[i]*s
            CC_τ += W*c*c
        end
        # P[n] should be (abs2(YC_τ)/CC_τ + abs2(YS_τ)/SS_τ)/YY, but we guard the values of
        # CC_τ and SS_τ in case they're zeros (note that SS_τ = 1 - CC_τ).
        frac_C = ifelse(iszero(CC_τ), nil, abs2(YC_τ) / CC_τ)
        frac_S = ifelse(CC_τ == 1,    nil, abs2(YS_τ) / (1 - CC_τ))
        P[n] = (frac_C + frac_S)/YY
    end
    return P
end

_periodogram!(P::AbstractVector, times, p::GLSPlan_fit_mean) =
    _generalised_lombscargle!(P, p.freq, times, p.y, p.w, p.Y, p.YY, zero(eltype(p.P)))
_periodogram!(p::GLSPlan_fit_mean) = _periodogram!(p.P, p.times, p)
_periodogram!(P::AbstractVector, times, p::GLSPlan) =
    _generalised_lombscargle!(P, p.freq, times, p.y, p.w, p.YY, zero(eltype(p.P)))
_periodogram!(p::GLSPlan) = _periodogram!(p.P, p.times, p)
