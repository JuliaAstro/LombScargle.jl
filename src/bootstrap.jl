### bootstrap.jl
#
# Copyright (C) 2016 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle, bootstrapping
#
# This file is a part of LombScargle.jl.
#
# License is MIT "Expat".
#
### Commentary
#
# This file contains facilities to perform bootstrapping and to calculate
# false-alarm probability and its inverse.
#
### Code:

immutable Bootstrap{T<:AbstractFloat}
    p::AbstractVector{T} # Vector of highest peaks
end

# XXX: possible improvement: compute IFFT plan for the fast method only once.
# However, this would require some work, instead current implementation is
# extremely simple.
function bootstrap{R1<:Real,R2<:Real}(N::Integer,
                                      t::AbstractVector{R1},
                                      s::AbstractVector{R2},
                                      rest...;
                                      args...)
    # Allocate vector
    high_peaks = Vector{promote_type(R1, R2)}(N)
    # Run N simulations.
    @inbounds for i in eachindex(high_peaks)
        # Store the highest peaks
        high_peaks[i] = findmaxpower(lombscargle(shuffle(t), s, rest...; args...))
    end
    # Create a `Bootstrap' object with the vector sorted in descending order.
    return Bootstrap(sort(high_peaks, rev = true))
end

function bootstrap{R<:Real,F<:AbstractFloat}(N::Integer,
                                             t::AbstractVector{R},
                                             s::AbstractVector{Measurement{F}};
                                             args...)
    return bootstrap(N, t, value.(s), uncertainty.(s); args...)
end

"""
    LombScargle.bootstrap(N::Integer,
                          times::AbstractVector{Real},
                          signal::AbstractVector{Real},
                          errors::AbstractVector{Real}=ones(signal); ...)

Create `N` bootstrap samples, perform the Lomb–Scargle analysis on them, and
store all the highest peaks for each one in a `LombScargle.Bootstrap` object.
All the arguments after `N` are passed around to `lombscargle`, which see.
"""
bootstrap


"""
    fap(b::Bootstrap, power::Real)

Return the false-alarm probability for `power` in the bootstrap sample `b`.

Its inverse is the `fapinv` function.
"""
fap{T<:AbstractFloat}(b::Bootstrap{T}, power::Real) =
    length(find(x -> x >= power, b.p))/length(b.p)

"""
    fapinv(b::Bootstrap, prob::Real)

Return the power value whose false-alarm probability is `prob` in the bootstrap
sample `b`.

It returns `NaN` if the requested probability is too low and the power cannot be
determined with the bootstrap sample `b`.  In this case, you should enlarge your
bootstrap sample so that `N*fap` can be rounded to an integer larger than or
equal to 1.

This is the inverse of `fap` function.
"""
fapinv{T<:AbstractFloat}(b::Bootstrap{T}, prob::Real) =
    get(b.p, round(Int, length(b.p)*prob), NaN)
