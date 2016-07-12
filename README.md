# LombScargle

[![Build Status](https://travis-ci.org/giordano/LombScargle.jl.svg?branch=master)](https://travis-ci.org/giordano/LombScargle.jl) [![Build status](https://ci.appveyor.com/api/projects/status/vv6mho713fuse6qy/branch/master?svg=true)](https://ci.appveyor.com/project/giordano/lombscargle-jl/branch/master) [![Coverage Status](https://coveralls.io/repos/github/giordano/LombScargle.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/LombScargle.jl?branch=master) [![codecov](https://codecov.io/gh/giordano/LombScargle.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/giordano/LombScargle.jl) [![LombScargle](http://pkg.julialang.org/badges/LombScargle_0.4.svg)](http://pkg.julialang.org/?pkg=LombScargle) [![LombScargle](http://pkg.julialang.org/badges/LombScargle_0.5.svg)](http://pkg.julialang.org/?pkg=LombScargle)

Introduction
------------

`LombScargle.jl` is a [Julia](http://julialang.org/) package to estimate the
[frequency spectrum](https://en.wikipedia.org/wiki/Frequency_spectrum) of a
periodic signal with
[the Lomb–Scargle periodogram](https://en.wikipedia.org/wiki/The_Lomb–Scargle_periodogram).

Another Julia package that provides tools to perform spectral analysis of
signals is [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl), but its methods
require that the signal has been sampled at equally spaced times.  Instead, the
Lomb–Scargle periodogram enables you to study unevenly sampled data, which is a
fairly common case in astronomy.

The algorithm used in this package are reported in the following papers:

* Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
  http://dx.doi.org/10.1088/0067-0049/191/2/247, Bibcode:
  http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
* Zechmeister, M., Kürster, M. 2009, A&A, 496, 577 (URL:
  http://dx.doi.org/10.1051/0004-6361:200811296, Bibcode:
  http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)

Installation
------------

`LombScargle.jl` is available for Julia 0.4 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.add("LombScargle")
```

You may need to update your package list with `Pkg.update()` in order to get the
latest version of `LombScargle.jl`.

Usage
-----

After installing the package, you can start using it with

```julia
using LombScargle
```

The module defines a new `LombScargle.Periodogram` data type, which, however, is
not exported because you will most probably not need to manually construct
`LombScargle.Periodogram` objects.  This data type holds both the frequency and
the power vectors of the periodogram.

The main function provided by the package is `lombscargle`:

```julia
lombscargle{T<:Real}(times::AbstractVector{T}, signal::AbstractVector{T},
                     freqs::AbstractVector{T}=autofrequency(times),
                     errors::AbstractVector{T}=ones(signal);
                     center_data::Bool=true, fit_mean::Bool=true)
```

The mandatory arguments are:

* `times`: the vector of observation times
* `signal`: the vector of observations associated with `times`

Optional arguments are:

* `freqs`: the frequency grid at which the periodogram will be computed.  If not
  provided, the grid will be computed with `autofrequency` function (see below)
* `errors`: the uncertainties associated to each `signal` point

Optional keywords arguments are:

* `fit_mean`: if `true`, fit for the mean of the signal using the Generalised
  Lomb-Scargle algorithm.  If this is `false`, the original algorithm by Lomb
  and Scargle will be employed, which does not take into account a non-null mean
  and uncertainties for observations
* `center_data`: if `true`, subtract the mean of `signal` from `signal` itself
  before performing the periodogram.  This is especially important if `fit_mean`
  is `false`

### Access Frequency Grid and Power Spectrum of the Periodogram ###

```julia
power(p::Periodogram)
freq(p::Periodogram)
freqpower(p::Periodogram)
```

`lombscargle` function return a `LombScargle.Periodogram` object, but you most
probably want to use the frequency grid and the power spectrum.  You can access
these vectors with `freq` and `power` functions, just like in `DSP.jl` package.
If you want to get the 2-tuple `(freq(p), power(p))` use the `freqpower`
function.

### Find Frequencies with Highest Power ###

```julia
findmaxfreq(p::Periodogram)
findmaxfreq(p::Periodogram, threshold::Real)
```

Return the frequencies with the highest power in the periodogram `p`.  If a
second argument `threshold` is provided, return the frequencies with power
larger than or equal to `threshold`.

### Determine Frequency Grid ###

```julia
autofrequency(times::AbstractVector{Real};
              samples_per_peak::Int=5,
              nyquist_factor::Integer=5,
              minimum_frequency::Real=NaN,
              maximum_frequency::Real=NaN)
```

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
