### townsend.jl ---  Perform original Lomb–Scargle periodogram
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

struct LSPlan{T,A,B<:AbstractVector{T},C,D,E,F,G} <: PeriodogramPlan
    times::A
    signal::B
    freq::C
    sumw::D
    w::E
    X::B
    YY::T
    noise::F
    norm::Symbol
    P::G
end

# Original algorithm that doesn't take into account uncertainties and doesn't
# fit the mean of the signal.  This is implemented following the recipe by
# * Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
#   http://dx.doi.org/10.1088/0067-0049/191/2/247,
#   Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
function _lombscargle_orig!(P::AbstractVector{T}, times::AbstractVector{<:Real},
                            X::AbstractVector{<:Real}, freqs::AbstractVector{<:Real},
                            XX::Real) where {T<:Real}
    N = length(X)
    @inbounds Threads.@threads for n in eachindex(freqs)
        ω = freqs[n] * 2 * pi
        XC = XS = CC = CS = zero(T)
        @inbounds for j in eachindex(times)
            S, C = sincos(ω * times[j])
            XC  += X[j]*C
            XS  += X[j]*S
            CC  += C*C
            CS  += C*S
        end
        SS       = N - CC
        s_τ, c_τ = sincos(atan2(CS, CC - N / 2) / 2)
        c_τ2     = c_τ*c_τ
        s_τ2     = s_τ*s_τ
        cs_τ_CS  = 2c_τ*s_τ*CS
        P[n]     = (abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                    abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))/XX
    end
    return P
end

_periodogram!(p::LSPlan) = _lombscargle_orig!(p.P, p.times, p.X, p.freq, p.YY)
_periodogram!(times, p::LSPlan) = _lombscargle_orig!(p.P, times, p.X, p.freq, p.YY)
