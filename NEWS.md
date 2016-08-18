History of LombScargle.jl
=========================

v0.1.0 (2016-08-18)
-------------------

### New Features ###

* The fast, but approximate, method by Press & Rybicki (1989) has been
  implemented.  This has complexity O[N log(N)], to be compared with the O[N^2]
  of the true Lomb-Scargle periodogram.  The implementation in this package is
  based on the one in Astropy.  Related to this method, three new keywords for
  `lombscargle` function has been added: `fast`, `oversampling`, `Mfft`.
* The generalised Lomb–Scargle algorithm by Zechmeister & Kürster is used also
  with `fit_mean=false`, when the user provided the uncertainties.

v0.0.2 (2016-08-05)
-------------------

### New Features ###

* New functions: `findmaxpower`, `prob`, `probinv`, `fap`, `fapinv`.
* New optional keyword for `lombscargle` function: `noise_level`.

v0.0.1 (2016-07-16)
-------------------

Initial release.
