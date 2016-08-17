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
### Code:

function add_at!{R1<:Real,R2<:Real}(arr::AbstractVector{R1},
                                    ind::AbstractVector{R2},
                                    val::Real)
    for i in ind
        arr[mod(i, length(arr) - 1) + 1] += val
    end
end

function add_at!{R1<:Real,R2<:Real,R3<:Real}(arr::AbstractVector{R1},
                                             ind::AbstractVector{R2},
                                             vals::AbstractVector{R3})
    for i in eachindex(ind)
        arr[mod(ind[i], length(arr) - 1) + 1] += vals[i]
    end
end

function extirpolate{R1<:Real,R2<:Real}(X::AbstractVector{R1},
                                        Y::AbstractVector{R2},
                                        N::Integer=0, M::Integer=4)
    x, y = collect(X), collect(Y)
    if N <= 0
        N = trunc(Int, maximum(x) + 0.5*M + 1)
    end
    result = zeros(R2, N)
    integers = find(isinteger, x)
    add_at!(result, trunc(Int, x[integers]), y[integers])
    deleteat!(x, integers)
    deleteat!(y, integers)
    ilo = clamp(trunc(Int, x - div(M, 2)), 0, N - M)
    numerator = y .* [prod(x[j] - ilo[j] - (0:M-1)) for j in eachindex(x)]
    denominator = factorial(M - 1)
    for j in 0:(M - 1)
        if j > 0
            denominator *= j/(j - M)
        end
        ind = ilo + (M - 1 - j)
        add_at!(result, ind, numerator ./ (denominator * (x .- ind)))
    end
    return result
end
