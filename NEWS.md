# History of LombScargle.jl

## v1.0.0 (2020-11-22)

### New Features

* Previously, when using the fast method the vector of the signal was padded to
  a length which is a power of 2, but now you can choose the factors of the
  length of the padded vector with the new keyword argument for
  `LombScargle.lombscargle` and `LombScargle.plan` functions:
  `padding_factors::Union{NTuple{N,<:Integer} where {N},Vector{<:Integer}}`.
  This defaults to `[2,3,5,7]`, which is the optimal choice for FFTW and allows
  for smaller vectors compared to powers of 2.  To reproduce the same results as
  with the previous default setting you need to use `padding_factors=[2]`.

## v0.5.1 (2020-05-15)

### Bug Fixes

* Minor reorganisation of the bootstraping code.

## v0.5.0 (2019-12-08)

### Breaking Changes

* Support for Julia 0.7 was dropped, now the minimum version required is Julia
  v1.0.

## v0.4.0 (2018-08-23)

### Breaking Changes

* Support for Julia 0.6 was dropped.

## v0.3.1 (2018-07-28)

* Minor performance improvements
* New documentation at https://juliaastro.github.io/LombScargle.jl/stable/

## v0.3.0 (2017-04-29)

### New Features

* You can pre-plan a periodogram before actually performing it using
  `LombScargle.plan` function.  This computes some quantities needed afterwards
  and pre-allocate the memory for the actual computation of the periodogram.
  The speed-up is particularly relevant for the fast method.
* Add `flags` and `timelimit` optional keywords to `lombscargle` function, to
  set the FFTW planner flags and the time limit.
* Package license changed to BSD 3-clause "New" or "Revised".

### Breaking Changes

* Support for Julia 0.4 and 0.5 was dropped.
* The `normalization` keyword of `lombscargle` function now must be a `Symbol`
  (instead of `AbstractString`), the default being `:standard` (instead of
  `"standard"`).  The same normalizations as before are supported, the names
  kept the same capitalization.

### Improvements

This version faced several performance improvements, in particular to
`lombscargle` and `LombScargle.model` functions, besides the pre-planning
feature.

* The fast method of `lombscarlge` now is faster, the larger the size of input
  array, the larger the improvement.  In addition, the fast Fourier transform
  computed internally with FFTW library can take advantage of multi-threading
  (call `FFTW.set_num_threads(n)` to use `n` threads) in order to speed-up
  computation.  However, please note that the running time will not scale as `n`
  because computation of the FFT is only a part of the algorithm.  The memory
  footprint of this function is also considerably lower.  To give you an idea of
  the improvement, for an input of 100000 datapoints, a pre-planned periodogram
  is 70% faster than a (non-planned) periodogram in previous version and
  requires almost 90% less memory.  With 4 FFTW threads the speed-up is of 80%.
* The two non-fast methods now are about 20%-30% faster, thanks to the use of
  `sincos` function from math library.  These methods now also support Julia’s
  native
  [multi-threading](http://docs.julialang.org/en/stable/manual/parallel-computing/#multi-threading-experimental).
  Run Julia with `n` threads (e.g., `JULIA_NUM_THREADS=4 julia` for 4 threads)
  in order to gain an `n`-fold scaling.  These functions also eat considerably
  less memory: if the periodogram is pre-planned, all operations are then
  performed in-place, so memory usage of the periodogram only is independent of
  input size.
* The `LombScargle.model` function is now a bit faster and less memory-greedy
  than before.

### Bug Fixes

* PSD normalization with heteroskedastic errors has been fixed.

## v0.2.0 (2016-12-07)

### Breaking Changes

* The fast method is now used when the frequency grid is evenly spaced (a
  `Range` object), no matter what the `times` vector is.  The previous behavior
  was due to a wrong interpretation of the applicability of the method.

### Bug Fixes

* `Periodogram` type now has 4 parameters, one for the type of each field.  Now
  `power`, `freq`, and `times` fields need not to have all the same floating
  point type.
* In the non-fast variant of the Generalised Lomb–Scargle method, when
  `fit_mean` and/or `center_data` are `true`, pre-center the data by subtracting
  from the signal the weighted average of the signal itself, instead of the
  arithmetic mean.

## v0.1.2 (2016-10-17)

### New Features

* New function for performing bootstrap resampling: `LombScargle.bootstrap`.
  The `fap` and `fapinv` functions have now new methods to estimate false-alarm
  probability and its inverse from a bootstrap sample.
* New utilities: `period`, `periodpower`, `findmaxperiod`.

### Bug Fixes

* Fix power in the standard (i.e., `fast = false` variant) generalised
  Lomb–Scargle algorithm with `fit_mean = true`.  You will find different
  results than before, but for the better, previous results were slightly wrong.

## v0.1.1 (2016-08-20)

### New Features

* New function: `LombScargl.model`.  It gives the best fitting Lomb–Scargle
  model at a given frequency.
* `findmaxfreq` function now can take an optional argument (`interval`) to limit
  the search for the maximum frequency to a certain frequency range.

## v0.1.0 (2016-08-18)

### New Features

* The fast, but approximate, method by Press & Rybicki (1989) has been
  implemented.  This has complexity O[N log(N)], to be compared with the O[N^2]
  of the true Lomb–Scargle periodogram.  The implementation in this package is
  based on the one in Astropy.  Related to this method, three new keywords for
  `lombscargle` function has been added: `fast`, `oversampling`, `Mfft`.
* The generalised Lomb–Scargle algorithm by Zechmeister & Kürster is used also
  with `fit_mean=false`, when the user provided the uncertainties.

## v0.0.2 (2016-08-05)

### New Features

* New functions: `findmaxpower`, `prob`, `probinv`, `fap`, `fapinv`.
* New optional keyword for `lombscargle` function: `noise_level`.

## v0.0.1 (2016-07-16)

Initial release.
