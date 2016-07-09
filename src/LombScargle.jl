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

export lombscargle

# Original algorithm that doesn't take into account uncertainties and doesn't
# fit the mean of the signal.  This is implemented following the recipe by
# * Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
#   http://dx.doi.org/10.1088/0067-0049/191/2/247)
function _lombscargle_orig{T<:Real}(times::AbstractVector{T}, signal::AbstractVector{T},
                                    angfreqs::AbstractVector{T}, center_data::Bool)
    P = Vector{T}(angfreqs)
    # If "center_data" keyword is true, subtract the mean from each point.
    signal_mean = center_data ? mean(signal) : zero(T)
    @inbounds for n in eachindex(angfreqs)
        ω = angfreqs[n]
        XC = XS = CC = SS = CS = zero(float(T))
        for j in eachindex(times)
            ωt = ω*times[j]
            C = cos(ωt)
            S = sin(ωt)
            X = signal[j] - signal_mean
            XC += X*C
            XS += X*S
            CC += C*C
            SS += S*S
            CS += C*S
        end
        τ       = 0.5*atan2(2CS, CC - SS)/ω
        c_τ     = cos(ω*τ)
        s_τ     = sin(ω*τ)
        c_τ2    = c_τ*c_τ
        s_τ2    = s_τ*s_τ
        cs_τ_CS = 2c_τ*s_τ*CS
        P[n] = 0.5*(abs2(c_τ*XC + s_τ*XS)/(c_τ2*CC + cs_τ_CS + s_τ2*SS) +
                    abs2(c_τ*XS - s_τ*XC)/(c_τ2*SS - cs_τ_CS + s_τ2*CC))
    end
    return angfreqs, P
end

function lombscargle{T<:Real}(times::AbstractVector{T}, signal::AbstractVector{T},
                              angfreqs::AbstractVector{T};
                              center_data::Bool=true)
    @assert length(times) == length(signal)
    return _lombscargle_orig(times, signal, angfreqs, center_data)
end

end # module
