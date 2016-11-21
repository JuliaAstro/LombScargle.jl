History of LombScargle.jl
=========================

v0.2.0 (201?-??-??)
-------------------

### Breaking Changes

* The fast method is now used when the frequency grid is evenly spaced (a
  `Range` object), no matter what the `times` vector is.  The previous behavior
  was due to a wrong interpretation of the applicability of the method.

### Bug Fixes

* `Periodogram` type now has 4 parameters, one for the type of each field.  Now
  `power`, `freq`, and `times` fields need need not to have all the same
  floating point type.
* In the non-fast variant of the Generalised Lomb–Scargle method, when
  `fit_mean` and/or `center_data` are `true`, pre-center the data by subtracting
  to the signal the weighted average of the signal itself, instead of the
  arithmetic mean.

v0.1.2 (2016-10-17)
-------------------

### New Features

* New function for performing bootstrap resampling: `LombScargle.bootstrap`.
  The `fap` and `fapinv` functions have now new methods to estimate false-alarm
  probability and its inverse from a bootstrap sample.
* New utilities: `period`, `periodpower`, `findmaxperiod`.

### Bug Fixes

* Fix power in the standard (i.e., `fast = false` variant) generalised
  Lomb–Scargle algorithm with `fit_mean = true`.  You will find different
  results than before, but for the better, previous results were slightly wrong.

v0.1.1 (2016-08-20)
-------------------

### New Features

* New function: `LombScargl.model`.  It gives the best fitting Lomb–Scargle
  model at a given frequency.
* `findmaxfreq` function now can take an optional argument (`interval`) to limit
  the search for the maximum frequency to a certain frequency range.

v0.1.0 (2016-08-18)
-------------------

### New Features

* The fast, but approximate, method by Press & Rybicki (1989) has been
  implemented.  This has complexity O[N log(N)], to be compared with the O[N^2]
  of the true Lomb–Scargle periodogram.  The implementation in this package is
  based on the one in Astropy.  Related to this method, three new keywords for
  `lombscargle` function has been added: `fast`, `oversampling`, `Mfft`.
* The generalised Lomb–Scargle algorithm by Zechmeister & Kürster is used also
  with `fit_mean=false`, when the user provided the uncertainties.

v0.0.2 (2016-08-05)
-------------------

### New Features

* New functions: `findmaxpower`, `prob`, `probinv`, `fap`, `fapinv`.
* New optional keyword for `lombscargle` function: `noise_level`.

v0.0.1 (2016-07-16)
-------------------

Initial release.
